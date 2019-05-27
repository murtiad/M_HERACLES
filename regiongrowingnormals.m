function [Regions] = regiongrowingnormals (ptCloud,normals2,h,AngleThres,CurvThres,neighbor_radius)
% REGIONGROWINGNORMALS
%
% Function to perform greedy region growing on a point cloud, based on the 
% normal and curvature values. See http://pointclouds.org/documentation/tutorials/region_growing_segmentation.php  
%
% Inputs: 
% - ptCloud: point cloud data 
% - normals2: normals of each point. Compute with pcnormals
% - h: mean curvature of each point. Compute using Beksi (2014)
% - AngleThres: threshold of normal angle to be considered as the same
% region
% - CurvThres: threshold of curvature to be considered as potential seed
% - neighbor_radius: search radius of nearest neighbor points
%
% Outputs:
% - Regions: a struct containing the segmented regions.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)
tic
A=ptCloud.Location; %available points
labels=zeros(ptCloud.Count,1);
j=1;
while ~all(isnan(A(:)))
    % name of this region
    thisRegionName=strcat('Region',int2str(j));
    
    % find the index of point having h nearest to 0
    seedID=find(h == min(h));
    seedGeom=A(seedID,:);
    seedNorm=normals2(seedID,:);
    
    % initialize the current region and current seed
    thisSeed=seedID;
    thisRegion=seedID;
    
    % loop for seed here
    while ~isempty(thisSeed)
    % find this seed's nearest neighbors
    
    [neighborIndices,~] = findNeighborsInRadius(ptCloud,seedGeom,neighbor_radius);
    
    % side quest: the above function also gives the seed ID as 'neighbor'.
    % must remove it from the matrix!
    id_dummy=neighborIndices==seedID;
    neighborIndices(id_dummy,:)=[];
    
    % initialize the computation for the neighbors
    [nbNeighbors,~]=size(neighborIndices);
    neighborGeom=zeros(nbNeighbors,3);
    neighborNorm=zeros(nbNeighbors,3);
    neighborCurv=zeros(1,nbNeighbors);
    neighborAngles=zeros(nbNeighbors,1);
    
    % start loop to compute neighbor's parameters
    for k=1:nbNeighbors
        neighborGeom(k,:) = A(neighborIndices(k),:);
        if ~isnan(neighborGeom(k,:))
            neighborNorm(k,:) = normals2(neighborIndices(k),:);
            neighborCurv(k) = h(neighborIndices(k));
        else
            continue
        end
        %compute normal angle between neighbors and seed (in degrees)
        neighborAngles(k,:) = rad2deg(atan2(norm(cross(seedNorm,neighborNorm(k,:))), seedNorm*neighborNorm(k,:)'));
        %if the angle is lower than the thresold, put into the current region
        if neighborAngles(k,:)<AngleThres
            thisRegion=[thisRegion;neighborIndices(k,1)];
            thisRegion=sort(thisRegion);
            %if the curvature is lower than the threshold, put in the current seed
            if neighborCurv(k)<CurvThres
                thisSeed=[thisSeed;neighborIndices(k,1)];
                thisSeed=sort(thisSeed);
            end
            %remove neighbor points from original ptCloud
            for l=1:size(thisRegion)
                A(thisRegion(l),:)=NaN;
                h(thisRegion(l))=NaN;
            end
        else
            continue
        end
               
    end
    
    %remove seed point from original ptCloud
    A(seedID,:)=NaN;
    h(seedID)=NaN;
    
    %remove current seed point from seed list
    id_dummy=find(thisSeed==seedID);
    thisSeed(id_dummy,:)=[];
    
    if ~isempty(thisSeed)
    seedID=thisSeed(1);
    seedGeom=ptCloud.Location(seedID,:);
    seedNorm=normals2(seedID,:);
    else
        continue
    end
    
    end

    ptCloudOut=select(ptCloud,thisRegion);
    if ptCloudOut.Count>50
        Regions.(thisRegionName)=ptCloudOut;
        for l=1:ptCloudOut.Count
            labels(thisRegion(l,1),1)=j+1;
        end
    else
        
        for l=1:ptCloudOut.Count
            labels(thisRegion(l,1),1)=1;
        end
        continue
    end
    
    j=j+1;
end

figure('name','Segmented Point Cloud')
pcshow(ptCloud.Location,labels,'MarkerSize',20)
colormap(hsv(j))
toc