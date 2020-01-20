function [Beams] = beamSeg (ptCloud,normals2,h,beam_width)
% BEAMSEG
%
% Function to perform greedy region growing on a point cloud, based on the 
% normal and curvature values. See http://pointclouds.org/documentation/tutorials/region_growing_segmentation.php  
%
% Inputs: 
% - ptCloud: point cloud data 
% - normals2: normals of each point. Compute with pcnormals
% - h: mean curvature of each point. Compute using Beksi (2014)
% - beam_width: width of the beam... this is the only manual information
% that needs to be provided by the user
%
% Outputs:
% - Beams: a struct containing the segmented regions for each beam and
% their cuboid RANSAC parameters
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

tic
format long g

%% OCTREE-BASED REGION GROWING
[RegionsOct]=regiongrowingnormalsOct(ptCloud,normals2,h,20);

%% FILTER REGIONS, DETECT IF THERE IS L OR Y-SHAPED FACES
threshold_pts = 1000;

nbRegions=numel(fieldnames(RegionsOct));
j=1;

for i=1:nbRegions
    %name of this region
    thisRegionName=strcat('Region',int2str(i));
    
    %point cloud of this region
    thisRegion=RegionsOct.(thisRegionName);
    
    %if the number of points is less than the threshold, ignore
    if thisRegion.Count<threshold_pts
        continue
    end
    
    %if the number of points is more than the threshold, perform axis
    %detection
    [Clusters,Axes,~]=axedetect(thisRegion,beam_width);
    
    %if no axes are detected, use entire point cloud as region
    if isempty(Axes)
        NewRegions(j).ptCloud=thisRegion;
        j=j+1;
    %otherwise, create new regions from the axis-constrained regions
    else
        for k=1:length(Clusters)
            NewRegions(j).ptCloud=Clusters(k).ptCloud;
            j=j+1;
        end
    end
end

% display the new regions (with multiple-axis regions separated using
% axedetect.m)
figure
for i=1:numel(NewRegions)
    thisPtCloud=NewRegions(i).ptCloud;
    labels=zeros(thisPtCloud.Count,1);
    labels(:,1)=1;
    for l=1:thisPtCloud.Count
            labels(l,1)=i+1;
    end
    pcshow(thisPtCloud.Location,labels)
    colormap(hsv(i+1))
    hold on
end

%% CREATE ADJACENCY MATRIX (CONSTRAINT 1: NEIGHBORHOOD)
% that is, determine each region's neighbors
clear adjM;

%threshold of point to plane distance
neighbor_thres = 0.2;
nbRegions=numel(NewRegions);

%initialise the adjacency matrix
adjM=zeros(nbRegions,nbRegions);
f = waitbar(0,'Creating adjacency matrix...','Name','beamSeg.m');

for i=1:nbRegions
    
    % set the reference region
    RefRegion=NewRegions(i).ptCloud;
    
    for j=i+1:nbRegions
        %set the check region
        CheckRegion=NewRegions(j).ptCloud;
        
        %create an octree of the check region to fasten the computation
        OT = OcTree(CheckRegion.Location,'binCapacity',50);  
        [binsIsi, ~, ~] = unique(OT.PointBins);
        [r,~]=size(binsIsi);
        
        %compute the octree bins' centroids coordinates
        binCentroid=zeros(r,3);
        for k=1:r
            x_bin=(OT.BinBoundaries(binsIsi(k),4)+ ...
                OT.BinBoundaries(binsIsi(k),1))/2;
            y_bin=(OT.BinBoundaries(binsIsi(k),5)+ ...
                OT.BinBoundaries(binsIsi(k),2))/2;
            z_bin=(OT.BinBoundaries(binsIsi(k),6)+ ...
                OT.BinBoundaries(binsIsi(k),3))/2;
            
            binCentroid(k,1)=x_bin;
            binCentroid(k,2)=y_bin;
            binCentroid(k,3)=z_bin;
        end
        
        %compute C2C distance between the check Octree bin centroids and
        %the reference point cloud
        dist=euclideanDistanceTwoPointClouds(RefRegion.Location,...
            binCentroid);
        
        %identify the closes distance between a check's octree centroid and
        %the a point from the reference point cloud
        min_dist=min(dist);
        
        %if the minimal distance is closer than the neighbor thresold, than
        %they are neighbors!
        if min_dist < neighbor_thres 
                adjM(i,j)=1;
                adjM(j,i)=1;
        end
    end
    waitbar(i/nbRegions,f)
end
close(f);

% draw a graph to visualise the neighborhood
a=2*pi/nbRegions;
t = 0.05:a:2*pi;
x1 = cos(t);
y1 = sin(t);
xy = [x1',y1'];
if length(xy)>nbRegions
    xy(nbRegions+1,:)=[];
end
figure
gplot(adjM,xy,'.-')
for i=1:nbRegions
    thisRegion=strcat('NewRegions',num2str(i));
    text(x1(i),y1(i),thisRegion,'HorizontalAlignment','center')
end
axis([-1 2 -1 2],'off')

%% ADD THE PARALLELENCY CONSTRAINT (CONSTRAINT 2: PARALLELISM)
%using the same principle as region growing to create the clusters

%create a (dummy) list of available regions
listRegions=1:1:nbRegions;
i=1;

%create a duplicate of the adjacency matrix
adjM_work=adjM;

%initialise the index of the beams
BeamIndex=zeros(nbRegions,2);
BeamIndex(:,1)=listRegions;
seedID=1; %FIRST seed ID

while ~all(isnan(listRegions))
    %immediately, put the current seed out of the function
    listRegions(seedID)=NaN;
    
    %...and put its index in the beam index matrix
    BeamIndex(seedID,2)=i;
    
    %set the seed point cloud
    seedPtCloud = NewRegions(seedID).ptCloud;
    
    %compute its PCA
    seedPCA=pca(seedPtCloud.Location);
    
    %...and extract the 1st order 
    seedVec=[seedPCA(1,1);seedPCA(2,1);seedPCA(3,1)];
    
    %put our seed into the seed list
    listSeed=find(adjM_work(seedID,:)==1);
    
    %put our seed as NaN in the adjacency matrix duplicate
    adjM_work(:,seedID)=NaN;
    adjM_work(seedID,:)=NaN;
    
    %keep on looping until no seeds are left
    while ~isempty(listSeed)
        k=length(listSeed);
        
        %set a secondary dummy seed list
        listSeed2=listSeed;
        
        for j=1:k
            %set the new seed ID (actually the first seed's neighborand 
            %point cloud, and compute also its 1st order PCA parameters
            newSeedID=listSeed(j);
            newSeedPtCloud = NewRegions(newSeedID).ptCloud;
            newSeedPCA = pca(newSeedPtCloud.Location);
            newSeedVec = [newSeedPCA(1,1);newSeedPCA(2,1);newSeedPCA(3,1)];
            
            %definition of parallelism: their axes' vector follows the
            %following relation: v2 = k*v1
            %since the PCA vectors' magnitude is 1 (normalised vectors),
            %the disparity between the vectors should be zero if they are
            %parallel
            disparity=abs(newSeedVec-seedVec);
            magDisparity=sqrt(disparity(1)^2+disparity(2)^2+disparity(3)^2);
            
            %set a tolerance, if the disparity between vectors is less
            %than...
            if magDisparity < 0.1
                BeamIndex(newSeedID,2)=i;
                
                %...means than they are parallel. Now let's check if they
                %are neighbors
                neighbors=find(adjM_work(newSeedID,:)==1);
                
                %if they are neighbors, add the neighbors to the list of
                %seed. Otherwise don't add them
                if ~isempty(neighbors)
                    listSeed2=[listSeed,neighbors];
                end
                
                %'kill' this particular neighbor in the list of regions and
                %the adjacency matrix
                listRegions(newSeedID)=NaN;
                
                adjM_work(:,newSeedID)=NaN;
                adjM_work(newSeedID,:)=NaN;
            end
            
            %'kill' this particular neighbor in the dummy seed list
            idx = listSeed2==newSeedID;
            listSeed2(idx) = [];
            
        end
        %update the seed list
        listSeed=listSeed2;
        
    end
    %find the index of a non-NaN (unkilled) element in the adjacency matrix
    indexNonNan=find(~isnan(adjM_work));
    
    %take the smallest value of indexing
    indexNonNan=min(indexNonNan);
    
    %this can be improved. I used rem because find returns a matricial
    %index instead of row/column index, while I only want the row...
    %maybe use listRegion instead... for the moment too lazy to do it
    indexNonNan=rem(indexNonNan,nbRegions);
    if indexNonNan==0
        indexNonNan=nbRegions;
    end
    
    %renew the seed ID with the newly found index
    seedID=indexNonNan;
    
    %carry on....
    i=i+1;    
end

%so how many beams do we have?
nbBeams=max(BeamIndex(:,2));

%wait a minute! if only one face is present, it cannot constitute a beam
%(needs at least two beam faces)
for i=1:nbBeams
    nbFacesofBeams = sum(BeamIndex(:,2) == i);
    if nbFacesofBeams<2
        idx=BeamIndex(:,2)==i;
        BeamIndex(idx,:)=[];
    end
end

%so here we have the final number of beams
nbBeams=max(BeamIndex(:,2));
disp(strcat('Oh happy day! Found',32,num2str(nbBeams),32,'beam(s)!'));

%% FINAL STEPS: MERGE THE BEAM FACES PT CLOUD TO CREATE THE BEAM PT CLOUD
for i=1:nbBeams
    %find the index of the beam face
    idx=BeamIndex(:,2)==i;
    thisBeamIndex=BeamIndex(idx,:);
    
    %compute the number of faces for this particular beam
    [nbFaces,~]=size(thisBeamIndex);
    
    %initialise the beam point cloud with the first face
    mergedPtCloud=NewRegions(thisBeamIndex(1,1)).ptCloud;
    
    %merge the reste progressively to the first face point cloud
    for j=2:nbFaces
        nextFacePtCloud=NewRegions(thisBeamIndex(j,1)).ptCloud;    
        mergedPtCloud=pcmerge(mergedPtCloud,nextFacePtCloud,0.001);
    end
    
    %store the result in a struct
    Beams(i).ptCloud=mergedPtCloud;
end

%visualise the result
figure
for i=1:numel(Beams)
    thisPtCloud=Beams(i).ptCloud;
    labels=zeros(thisPtCloud.Count,1);
    labels(:,1)=1;
    for l=1:thisPtCloud.Count
            labels(l,1)=i+1;
    end
pcshow(thisPtCloud.Location,labels)
colormap(hsv(i+1))
hold on
end

%% EXTRA STEP: CUBOID RANSAC
%use this part to generate geometric primitives (cuboids) from the
%now-segmented beams
figure
pcshow(ptCloud)
for i=1:length(Beams)
    ptCloud=Beams(i).ptCloud;
    coords=double(ptCloud.Location);
    [model, CuboidParameters, inlierIndices, outlierIndices] = ...
        CuboidRANSAC( coords );
    
    %store the cuboid parameters in the struct
    Beams(i).CuboidModel=model;
    Beams(i).CuboidParameters=CuboidParameters;
    Beams(i).CuboidInlierIdx=inlierIndices;
    Beams(i).CuboidOutlierIdx=outlierIndices;
    
    %visualise the cuboids
    DisplayModel(numel(inlierIndices), model);
    hold on
end
toc