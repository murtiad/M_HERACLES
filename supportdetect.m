function [Support,STable,remainPc] = supportdetect(Body_ptCloud)
% SUPPORTDETECT
%
% Function to detect building supports (columns and piers)
%
% Inputs: 
% - Body_ptCloud : the point cloud of the building's body without its
% attics. See also ATTICSEGMENT
%
% Outputs:
% - Support : a struct with the the individual support clusters as fields.
% Each field includes the point cloud and its classification, either a
% column or a pier.
% - STable : a recapitulative table with the main results tabled.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

%% Slice the point cloud and take the toppest 
pas = 0.2; %set the slicing interval
thc = 0.05; %set slice thickness
tol_ground=0.2; %set tolerance for min Z
[slices2,slices_list] = slices(Body_ptCloud,pas,thc,tol_ground);
[frontierSliceNb,~] = size(slices_list);

%% Use Euclidean nearest neighbor to detect "islands"
frontierSliceNm = strcat('SLICE',num2str(frontierSliceNb));
frontierSlicePtCld = slices2.(frontierSliceNm);

[nRow,~] = size(frontierSlicePtCld.Location);
    
%check Z max an Z min for this slice
Zmax = max(frontierSlicePtCld.Location(:,3));
Zmin = min(frontierSlicePtCld.Location(:,3));
    
%set an average Z to project the points to
Z=(Zmax+Zmin)/2;
    
%project the points to a plane Z
ProjectedPts = zeros(nRow,3); %initialise matrix
ProjectedPts(:,1)=frontierSlicePtCld.Location(:,1); %x
ProjectedPts(:,2)=frontierSlicePtCld.Location(:,2); %y
ProjectedPts(:,3)=Z; %projected z
    
%create a convexhull around the projected points (required further below)
[U,S]=svd(   bsxfun(@minus,ProjectedPts,mean(ProjectedPts)),   0);  
[~,SliceArea]=convhull(U*S(:,1:2)); %this is required below

%perform the Euclidean segmentation to detect the islands
distThreshold=0.05;
[labels,~] = pcsegdist(frontierSlicePtCld,distThreshold);
[Clusters,listClusters] = clustering(labels,frontierSlicePtCld,10);

%% determine if the clusters are wall or support (column/pier)
nbClusters = numel(fieldnames(Clusters));
%create a lookup table
STable = table("type",0,"type",0,"type");
STable.Properties.VariableNames={'ClusterNo' 'Area_m2' 'Type' 'Circularity' 'SupportType'};
[~, id] = lastwarn;
warning('off', id)

for i=1:nbClusters
    ClusterNm = listClusters{i};
    ThisCluster = Clusters.(ClusterNm);%name of this cluster
    [nRow,~] = size(ThisCluster.Location);

    %set Z as 0 (projected to Z=0)
    %doesn't matter, we only want to compute the surfaces
    Z = 0;
    
    %project the points to a plane Z=0
    ProjectedPts = zeros(nRow,3); %initialise matrix
    ProjectedPts(:,1)=ThisCluster.Location(:,1); %x
    ProjectedPts(:,2)=ThisCluster.Location(:,2); %y
    ProjectedPts(:,3)=Z; %projected z
    
    %create a convexhull around the projected points
    [U,S]=svd(   bsxfun(@minus,ProjectedPts,mean(ProjectedPts)),   0);  
    [k,area]=convhull(U*S(:,1:2));
    
    [ksize,~] = size(k);
    chull_pts = zeros(ksize,3);
    
    for j=1:ksize
        chull_pts(j,1)=ProjectedPts(k(j,1),1);
        chull_pts(j,2)=ProjectedPts(k(j,1),2);
        chull_pts(j,3)=ProjectedPts(k(j,1),3);
    end
        
    %stock the computed area in the table' 1st column, its Z in 2nd column 
    STable{i,1} = convertCharsToStrings(ClusterNm);
    STable{i,2} = area;
    
    % THIS IS THE IMPORTANT PART
    if area < 0.1*SliceArea %if cluster area is lower than 0.1*slice area, it's probably a support at this level
        STable{i,3}="support"; %support
    else
        STable{i,3}="wall"; %wall
    end
    
end

% extract the "support" table
rows = STable.Type=="support";
vars = {'ClusterNo','Area_m2'};
TbSupport = STable(rows,vars);

% extract the "wall" table
% maybe useful in the future
% rows = Table.Type=="wall";
% vars = {'ClusterNo','Area_m2'};
% TbWall = Table(rows,vars);

%% Determine the shape of a support
%%take the support table and detect the geometrical shape
%%if rectangular, it's a pier. If circular, it's a column
nbSupport = height(TbSupport);
disp(strcat('[DING!]',num2str(nbSupport),' supports were detected!'));
bufferSz=0.5; %set the buffer size for your cookie cutter. This parameter must be REFINED manually :(
f = waitbar(0,'Detecting structural support...');
for i=1:nbSupport
    ClusterNm = TbSupport.ClusterNo(i);
    ThisCluster = Clusters.(ClusterNm);
    [nRow,~] = size(ThisCluster.Location);
    
    %set Z as the max Z
    Z=max(ThisCluster.Location(:,3));
    
    %project the points to a plane Z
    ProjectedPts = zeros(nRow,3); %initialise matrix
    ProjectedPts(:,1)=ThisCluster.Location(:,1); %x
    ProjectedPts(:,2)=ThisCluster.Location(:,2); %y
    ProjectedPts(:,3)=Z; %projected z
    
    %create a convexhull around the projected points
    [U,S]=svd(   bsxfun(@minus,ProjectedPts,mean(ProjectedPts)),   0);  
    [k,area]=convhull(U*S(:,1:2));
    
    [ksize,~] = size(k);
    chull_pts = zeros(ksize,3);
    
    for j=1:ksize
        chull_pts(j,1)=ProjectedPts(k(j,1),1);
        chull_pts(j,2)=ProjectedPts(k(j,1),2);
        chull_pts(j,3)=ProjectedPts(k(j,1),3);
    end
    
    %lets compute the circularity and put it on the fourth column of the lookup table
    chull_poly=polyshape(chull_pts(:,1),chull_pts(:,2));
    [~, id] = lastwarn;
    warning('off', id)
    chull_perimeter=perimeter(chull_poly);
    c = (chull_perimeter.^2)./(4*pi*area);
    STable{i,4}=c;
    
    %Fill the lookup table with the classes, put it in the 5th column
    if c<1.12
        STable{i,5}="column";
    else
        STable{i,5}="pier";
    end

    %Let's get cooking!(pardon the pun)
    %Buffering of the chull, and make it a cookiecutter
    polyout1 = polybuffer(chull_poly,bufferSz,'JointType','miter','MiterLimit',2);
    
    %Check if there are points inside the polygon (boolean)
    TFin = isinterior(polyout1,Body_ptCloud.Location(:,1),Body_ptCloud.Location(:,2));
    
    %find row numbers for points which are inside the polygon
    rows = find(TFin(:,1)==1);
    
    %extract the points within our cookiecutter
    cookie = select(Body_ptCloud,rows);

    %find row numbers for points which are outside the polygon
    rowsd = find(TFin(:,1)==0);
    
    %extract the points within our cookiecutter
    dough = select(Body_ptCloud,rowsd);
    
    %Cleaning our cookie
    %Detect the ground plane and store the points in inliers
    maxDistance = 0.05;
    referenceVector = [0,0,1];
    [~,inliers,outliers] = pcfitplane(cookie,maxDistance,referenceVector); %only need the outlier
    
    %Cluster the points, ignoring the ground plane points. Why? to clean the cookie of course
    ptCloudWithoutGround = select(cookie,outliers);
    ptCloudGround = select(cookie,inliers);
    distThreshold = 0.1; %this must also be REFINED

    [labels,~] = pcsegdist(ptCloudWithoutGround,distThreshold);

    %Look for unique labels
    a = unique(labels);
    
    %...create a table with name of label and frequency
    labeltable = [a,histc(labels(:),a)];
    
    %look for the point cloud with the most points after the segmentation
    [aa,~]=sort(labeltable(:,2),'descend');
    
    %find the index of this largest point cloud in table
    rowTable = find(labeltable(:,2)==aa(1,1));
    indexTable = labeltable(rowTable(1,1),1);
    
    %find the index of points with this index in the original point cloud
    rows2 = find(labels(:,1)==indexTable);
    
    %extract the cleaned point cloud
    inPc2 = select(ptCloudWithoutGround,rows2);
    
    %Retrieve the outliers from the pcsegdist process
    rows3 = find(labels(:,1)~=indexTable);
    
    %extract the outlier point cloud
    outPc = select(ptCloudWithoutGround,rows3);
    
    %safeguard if the cluster with most point IS NOT the support
    %add another criterion based on height
    ZminInPc=min(cookie.Location(:,3)); %compute min Z of original cookie
    ZmaxInPc=max(cookie.Location(:,3)); %compute max Z of original cookie
    
    [r,~]=size(aa);
    
    %iterate until we get the cluster which meets our criteria
    for m=1:r
    ZmaxInPc2=max(inPc2.Location(:,3)); %compute max Z of cleaned cookie
    if ZmaxInPc2<ZminInPc+(0.8*(ZmaxInPc-ZminInPc)) %max Z must be lower than 0.8 of original cookie's height
        rowTable = find(labeltable(:,2)==aa(m+1,1));
        indexTable = labeltable(rowTable(1,1),1);
        rows2 = find(labels(:,1)==indexTable);
        inPc2 = select(ptCloudWithoutGround,rows2);
        outPc = select(ptCloudWithoutGround,rows3);
    end
    end
    
    %no harm in more cleaning up
    inPc2 = pcdenoise(inPc2);
    
    %merge the cookie crumbs (outlier point cloud)
    outPc2 = pcmerge(dough,ptCloudGround,0.01);
    outPc3 = pcmerge(outPc2,outPc,0.01);
    
    %store the point cloud in a structure called Object
    Object.PtCloud = inPc2;

    %Plot the results
    pointscolor=uint8(zeros(inPc2.Count,3));
    
    if STable.SupportType(i)=="column"    
        pointscolor(:,1)=255; %red
        Object.Type='Column'; %store the object type
    elseif STable.SupportType(i)=="pier"     
        pointscolor(:,3)=255; %blue
        Object.Type='Pier'; %store the object type
    end
    inPc2.Color=pointscolor;
    fig=figure(1);
    fig.Name='Support Detection';
    axis equal
    title('Support Detection. Red=Columns; Blue=Piers');
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
    pcshow(inPc2)
    hold on

    %store the object in a structure called Support
    Support.(ClusterNm) = Object;
    
    Body_ptCloud = outPc3;
    
    waitbar((i/nbSupport),f)
end
remainPc = Body_ptCloud;
close(f);
