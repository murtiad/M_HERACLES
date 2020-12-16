function [ClustersOut] = detectRoof(ptCloud,h)
% DETECTROOF
%
% Funtion to detect roof vertices from an aerial point cloud. 
%
% Inputs: 
% - ptCloud: point cloud data 
% - h: mean curvature of each point. Compute using Beksi (2014). If this
% value is not provided, the function will compute it (may take a while)
%
% Outputs:
% - ClustersOut: a struct containing the planar surfaces of the roof.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

format long g

normals=ptCloud.Normal;

%% COMPUTE CURVATURE IF NOT PROVIDED
if nargin<2
    tree = KDTreeSearcher(ptCloud.Location);
    radius = 1;
    f=waitbar(0,'Computing curvatures...');
    c1=zeros(1,ptCloud.Count);
    c2=zeros(1,ptCloud.Count);
    h=zeros(1,ptCloud.Count);
    for i=1:ptCloud.Count
        query = [ptCloud.Location(i,1) ptCloud.Location(i,2) ...
            ptCloud.Location(i,3)];
        [c1(i), c2(i)] = estimateCurvatures(normals, tree, query, radius);
        h(i) = (c1(i)+c2(i))/2; %Mean curvature; if = 0 means point is plane 
        waitbar(i/ptCloud.Count,f)
    end
    close(f)
end
% if needed, save the results
%savefile = strcat('curvature_detectRoof.mat');
%save(savefile, 'h');
%% NORMALS-BASED REGION GROWING
% detect the roof patches and store it in a struct
% see also: regiongrowingnormalsOct.m
[RegionsOct]=regiongrowingnormalsOct(ptCloud,normals,h);

%% DISTANCE-BASED REGION GROWING 
% this part is used to clean up the resulting patches from noises

% extract patch names
nbRegions=numel(fieldnames(RegionsOct));
f=fieldnames(RegionsOct);

% extract number of patches
RegionCount=nbRegions;
for i=1:nbRegions
    RegionName = f{i};
    labels=pcsegdist(RegionsOct.(RegionName),0.2);
    uniqueLabels=unique(labels);
    [nbUniqueLabels,~]=size(uniqueLabels);
    for j=1:nbUniqueLabels
        [~, boolInd]=ismember(labels,uniqueLabels(j,1));
        indexInd=find(boolInd);
        newRegion=select(RegionsOct.(RegionName),indexInd);
        
        % if the resulting patch consists of less than 1000 points, delete
        % it and do not use for further processing
        if newRegion.Count>1000
            RegionCount=RegionCount+1;
            newRegionName=strcat('Region',string(RegionCount));
            RegionsOct.(newRegionName)=newRegion;
        else
            continue
        end
    end
    RegionsOct=rmfield(RegionsOct,RegionName);
end

%% EXTRACT THE PLANE MESHES
f=fieldnames(RegionsOct);
nbRegions=numel(fieldnames(RegionsOct));
objectList=strings(nbRegions,1);
for i=1:nbRegions
    RegionName = f{i};
    thisPtCloud=RegionsOct.(RegionName);

    inPtCloud = thisPtCloud;
    
    % perform RANSAC to extract the plane
    [model,inlierIndices,~] = pcfitplane(inPtCloud,0.2);
    plane = select(inPtCloud,inlierIndices);

    %convert the Matlab surface model from pcfitplane to geom3d format
    a=model.Parameters(1,1);
    b=model.Parameters(1,2);
    c=model.Parameters(1,3);
    d=model.Parameters(1,4);

    P1 = [plane.Location(1,1) plane.Location(1,2) 0];

    P1(1,3) = (-d-(a*P1(1,1))-(b*P1(1,2)))/c;

    planeGeom3d = createPlane(P1, model.Normal); %plane in geom3d format
    
    %Perform 3D Delaunay Triangulation on the point cloud
    TR = delaunayTriangulation(double(plane.Location(:,1)),...
        double(plane.Location(:,2)),double(plane.Location(:,3)));
    
    %extract the surface of the Delaunay mesh TR (otherwise stored in
    %tetrahedral form)
    [~,xf] = freeBoundary(TR);
    [size_xf,~] = size(xf);
    
    %project the points to the mathematical surface model to get 1 surface
    %in the form of the RANSAC plane
    prjtdPts = zeros(size_xf,3);
    
    for j=1:size_xf
        %create a 3D ray for each foint in the mesh surface
        point = [xf(j,1), xf(j,2), xf(j,3)];
        
        %intersect the 3D ray with the mathematical surface
        prjtdPt = projPointOnPlane(point, planeGeom3d);
        
        %store the coordinates of the intersection
        prjtdPts(j,1) = prjtdPt(1,1);
        prjtdPts(j,2) = prjtdPt(1,2);
        prjtdPts(j,3) = prjtdPt(1,3);
    end
    
    %perform 3D Delaunay Triangulation on the projected points
    TR2 = delaunayTriangulation(double(prjtdPts(:,1)),...
        double(prjtdPts(:,2)),double(prjtdPts(:,3)));
    
    %extract the surface of the Delaunay mesh of the projected points
    [~,xf2] = freeBoundary(TR2);
    
    %further smoothing here
    [xf3,F2] = meshParallelSimp(xf2, 0.1);
    
    %store the mesh properties in a struct
    Mesh.Vertices = xf3;
    Mesh.Faces = F2;
    
    %create a Matlab alphashape
    shp = alphaShape(xf2(:,1),xf2(:,2),xf2(:,3),inf);

    %store all of this stuff in the another struct called Object
    objectList(i,1) = strcat('Object',num2str(i));
    Object.PtCloud = plane;
    Object.AlphaShp = shp;
    Object.Mesh = Mesh;
    Object.PlaneGeom3D = planeGeom3d;
    Object.PlaneMatlab = model;
    
    Clusters.(objectList{i}) = Object;
end

%% PERFORM 'SNAPPING'
SnappedCluster = meshSnap(Clusters,5);
%% REORDER THE VERTICE LIST TO COMPLY WITH CITYGML
% CityGML vertices must be sorted clockwise from barycentre

figure('name','Detected Roof Vertices - Patch/Polygon')
for i=1:nbRegions
    
    % name of current object
    ObjectName = objectList{i};
    
    % retrieve the object
    thisObject=SnappedCluster.(ObjectName);
    
    % retrieve the vertices list
    thisVertexList=thisObject.Mesh.Vertices;
    
    % initiate a list of bearing angles
    [nbVertices,~]=size(thisVertexList);
    listBearings=zeros(nbVertices,2);
    listBearings2=zeros(nbVertices,2);
    
    % compute the Principal Component Analysis
    coeffs = pca(thisVertexList);

    % transform the point cloud to PC-aligned system i.e. planar surface
    TransformPCA = thisVertexList*coeffs(:,1:3);
    
    % compute the barycenter of the vertices
    barycenter=[mean(TransformPCA(:,1)),mean(TransformPCA(:,2)), ...
        mean(TransformPCA(:,3))];
    
    % take only the XY (assume the points are coplanar)
    Xa=barycenter(1);
    Ya=barycenter(2);
    
    for j=1:nbVertices
        %indexing
        listBearings(j,1)=j;
        
        % retrieve the vertex's XY 
        Xb=TransformPCA(j,1);
        Yb=TransformPCA(j,2);
        
        % compute the bearing from the barycenter to each vertex
        [alpha, ~] = bearing_surv(Xa,Ya,Xb,Yb);
        
        listBearings(j,2)=alpha;
    end
    
    % sort the vertice list according to this bearing
    [sorted,index]=sort(listBearings(:,2),'ascend');
    listBearings2(:,1)=index;
    listBearings2(:,2)=sorted;
    
    % retrieve the original points, but sorted
    sortedVertices=zeros(nbVertices,3);
    for j=1:nbVertices
        sortedVertices(j,:)=thisVertexList(listBearings2(j,1),:);
    end
    
    % update the result with the new sorted values
     SnappedCluster.(ObjectName).Mesh.Vertices=sortedVertices;
    
    patch(sortedVertices(:,1),sortedVertices(:,2),sortedVertices(:,3),'g')
    view(3)
    hold on
    pcshow(SnappedCluster.(ObjectName).PtCloud)
    hold on
end

% copy the results to the output struct
ClustersOut=SnappedCluster;
end

