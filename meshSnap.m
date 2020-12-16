function [NewClusterStruct] = meshSnap(inputCluster,tol)
% MESHSNAP
%
% Function to merge nearby points within a set of meshes located within a
% distance of TOL.
%
% Inputs: 
% - inputCluster: a struct containing Objects, containing the point cloud
% and mesh information. The Objects are separate entities; the function
% will however snap the mesh vertices by taking into account the
% global/merged state of all Objects within the inputCluster
% - tol: distance threshold between a point with its neighbours to be
% considered as the same point
%
% Outputs:
% - vertices2: list of nx3 coordinates of the simplified mesh
% - faces2: triangular faces of vertices2 resulting froma a Delaunay
% triangulation
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

NewClusterStruct=inputCluster;

%% MERGE THE MESHES
nbRegions=numel(fieldnames(inputCluster));
objectList=strings(nbRegions,1);

% create a list containing Object names (useful later)
for j=1:nbRegions
    objectList(j,1)=strcat('Object',num2str(j));
    [nbVertices,~]=size(inputCluster.(objectList(j,1)).Mesh.Vertices);
    
    % create a new field with index to the original object
    v = zeros(nbVertices, 1);
    v(:) = j;
    inputCluster.(objectList(j,1)).Mesh.ObjectIndex = v;
end

% initialise variables
MergedVertices=inputCluster.Object1.Mesh.Vertices;
MergedFaces=inputCluster.Object1.Mesh.Faces;
MergedIndex=inputCluster.Object1.Mesh.ObjectIndex;
i=1;

% start the loop
while i<nbRegions
    objectName1=objectList(i);
    objectName2=objectList(i+1);
    
    % stack the vertex lists
    MergedVertices = vertcat(MergedVertices, inputCluster.(objectName2).Mesh.Vertices);

    % update the face list
    faceUpdate = inputCluster.(objectName2).Mesh.Faces+max(max(MergedFaces));
    MergedFaces = vertcat(MergedFaces,faceUpdate);
    
    % update the object index list
    MergedIndex = vertcat(MergedIndex, inputCluster.(objectName2).Mesh.ObjectIndex);
    
    i=i+1;
end

%% COMPUTE "SNAPPING"
%initialise the variables
ClustersCopy=inputCluster; % "dynamic" copy of the original Cluster struct
VertexList=MergedVertices; % "dynamic" list of vertices
MergedIndex2=MergedIndex; %"dynamic" list of indices

% start the loop, continue as long as the list of vertices is not empty
while ~isempty(VertexList)
    
    % take the first seed
    vertex1 = VertexList(1,:);
    
    % change VertexList to Matlab point cloud struct, necessary for the
    % next function
    VertexPointCloud=pointCloud(VertexList);
    
    % find the seed's 10 nearest neighbours (this value can be changed, but
    % I doubt it would be necessary)
    [indices,dists] = findNearestNeighbors(VertexPointCloud,vertex1,10);
    
    % for now, delete the seed from the list of nearest neighbours
    indices(1,:)=[];
    dists(1,:)=[];
    
    % number of neighbours inside the tolerance radius
    nb = length(indices);
    
    % initialise a dynamic list of neighbours and their pointer towards the
    % original Object
    NeighborList = []; % empty list of neighbours inside the radius
    OriginList = []; % empty list of neighbours' origins
  
    % start a loop to check each neighbour
    for i=1:nb
        % if the neighbour point's distance is less than the tolerance,
        % consider it as the same point i.e. "snapped"
        if dists(i,1)<tol
            a=MergedIndex2;
            NeighborList(i,1)=indices(i,1);
            OriginList(i,1)= a(indices(i,1));
        else
            continue
        end
    end
    
    % if the seed has neighbours...
    if ~isempty(NeighborList)
        % return the seed to both lists, necessary to compute the
        % coordinates of the new "snapped" point
        NeighborList=[NeighborList;1];
        OriginList=[OriginList;MergedIndex2(1)];
        
        % number of items in the list
        szList=length(NeighborList);

        % how many unique Objects?
        UniqueOrigins=unique(OriginList);
        szUnique=length(UniqueOrigins);

        % create a recapitulative table
        tableToDelete=zeros(szList,2);
        for x=1:szList
            originNm=strcat('Object',string(OriginList(x)));
    
            % retrieve the vector to be deleted from reference list of 
            % vertices
            VectToDelete =  VertexList(NeighborList(x),:);
    
            % find its index in the original Object vertices
            [~, boolDel]=ismember(ClustersCopy.(originNm).Mesh.Vertices,...
                VectToDelete);
            indexDel=find(boolDel,1);
    
            % update the recapitulative table
            tableToDelete(x,1) = OriginList(x);
            tableToDelete(x,2) = indexDel;
        end
        
        % delete the points whose indices are listed in NeighborList in 
        % respective original Objects
        for x=1:szUnique
            originNm=strcat('Object',string(UniqueOrigins(x)));
    
            % find the indices of This Object to be deleted from the 
            % recapitulative table
            ind1 = tableToDelete(:,1) == UniqueOrigins(x);
            thisOriginDel = tableToDelete(ind1,2);
    
            % delete the said point in the dynamic/copied cluster
            ClustersCopy.(originNm).Mesh.Vertices(thisOriginDel,:)=[];
    
            % keep only unique points and delete (possible) redundants
            ClustersCopy.(originNm).Mesh.Vertices=...
            unique(ClustersCopy.(originNm).Mesh.Vertices,'rows');
        end
        
        % initialise XYZ vectors of all the neighbours and the seed 
        x_vect=zeros(szList,1);
        y_vect=zeros(szList,1);
        z_vect=zeros(szList,1);
        
        % retrieve the neighbours and seed's XYZ
        for j=1:szList
            x_vect(j)=VertexList(NeighborList(j),1);
            y_vect(j)=VertexList(NeighborList(j),2);
            z_vect(j)=VertexList(NeighborList(j),3);          
        end
    
        % compute the median to generate a new point
        % why median? more robust than mean/average
        x_new = median(x_vect);
        y_new = median(y_vect);
        z_new = median(z_vect);
    
        new = [x_new y_new z_new];
        
        % add planarity constraint
        listPlanes=zeros(szUnique,9);
        listPlanesParameters=zeros(szUnique,4);
        
        % if the neighbours belong to more than 2 planes...
        if szUnique>2
            for k=1:szUnique
                originNm = strcat('Object',string(UniqueOrigins(k)));
                listPlanesParameters(k,:)=...
                    inputCluster.(originNm).PlaneMatlab.Parameters;
            end
            
            % use Geom3D's three planes intersection function
            %  plane1 = listPlanes(1,:);
            %  plane2 = listPlanes(2,:);
            %  plane3 = listPlanes(3,:);
            %  new2 = intersectThreePlanes(plane1, plane2, plane3);

            % use my least squares multiple planes function
            [new2,~]=intersectMultPlanes(listPlanesParameters);
        
        % if the neighbours belong to two planes...
        elseif szUnique==2
            for k=1:2
                originNm = strcat('Object',string(UniqueOrigins(k)));
                listPlanes(k,:)=inputCluster.(originNm).PlaneGeom3D;
            end
            plane1 = listPlanes(1,:);
            plane2 = listPlanes(2,:);
            
            % intersect the two planes to get a 3D line
            line = intersectPlanes(plane1, plane2);
            
            % project the computed median point into this 3D line
            new2 = projPointOnLine3d(new, line);
        else
            new2=new;
        end
        
        % stact the new point to the un-merged objects
        if szUnique>1
            for x=1:szUnique
                originNm = strcat('Object',string(UniqueOrigins(x)));
                ClustersCopy.(originNm).Mesh.Vertices= ...
                    [ClustersCopy.(originNm).Mesh.Vertices;new2];
            end
        else
            originNm = strcat('Object',string(UniqueOrigins(1)));
            ClustersCopy.(originNm).Mesh.Vertices=...
                [ClustersCopy.(originNm).Mesh.Vertices;new2];
        end
        
        % delete the snapped points from the two lists
        indNeighbors=NeighborList';
        VertexList(indNeighbors,:)=[];
        MergedIndex2(indNeighbors,:)=[];
    else        
        % delete the seed from the two lists if they don't have neighbours   
        VertexList(1,:)=[];
        MergedIndex2(1,:)=[];
    end
end

%% VISUALISE THE 'SNAPPED' RESULT AND STORE IN NEW STRUCT
figure('name','Detected Roof Vertices - Delaunay Mesh')
for x=1:nbRegions
    originNm=strcat('Object',string(x));
    
    % Delaunay triangulation
    TR = delaunayTriangulation...
        (double(ClustersCopy.(originNm).Mesh.Vertices(:,1)),...
        double(ClustersCopy.(originNm).Mesh.Vertices(:,2)),...
        double(ClustersCopy.(originNm).Mesh.Vertices(:,3)));
    
    % retrieve the free boundary mesh
    [F,xf] = freeBoundary(TR);
    [V2, F2] = ensureManifoldMesh(xf, F);
    
    % store the new vertices and faces in a new struct
    NewClusterStruct.(originNm).Mesh.Vertices=V2;
    NewClusterStruct.(originNm).Mesh.Faces=F2;
    
    % plot the mesh
    pcshow(ClustersCopy.(originNm).PtCloud)
    hold on
    drawMesh(V2,F2,'b');
    hold on
end

end

