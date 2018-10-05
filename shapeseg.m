function [Clusters,remainpc] = shapeseg(datapc,shpfilepath,bufferSz,distThreshold,gridSz)
% SHAPESEG
%
% Function to use (polygonal) shapefiles to delimit 3D point cloud area to
% be segmented using the pcsegdist function. 
%
% Inputs: 
% - datapc: point cloud data, preferably without ground. But it should also
% work with ground present, provided the shapefile is accurate enough
% - shpfilepath: file path to the shapefile with polygons and attributes
% - bufferSz: set buffer size for shapefile polygon
% - distThreshold: set max distance for pcsegdist
% - gridSz: set grid size for merging remaining point clouds. Use the
% resolution as a rule of thumb
%
% Outputs:
% - Clusters: a struct with the individual segmented objects as fields.
% Each field consists of another struct ("Object") which keeps the point
% cloud data as well as any eventual attributes as described in the
% shapefile
% - remainpc: the remaining point cloud after the objects are segmented
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

clear Clusters;
clear Object;

%Read shapefile data and convert it to native struct
[shapes,objectList,attList]=shpload(shpfilepath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters that may be modified (for test purposes)
% dthres = 1.0; %set denoising max distance (optional). default is ONE(1.0) sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[numObject,~] =size (objectList);
[numAtt,~] = size (attList);
f = waitbar(0,'Creating individual point cloud clusters...');

for i=1:numObject
    %Name of the object 
    objectNm = objectList(i,1);
    
    %Create buffer around the polyshape, to take into account digitising
    %inprecisions
    polyout1 = polybuffer(shapes.(objectNm).polyshape,bufferSz,'JointType','miter','MiterLimit',2);
    
    %Check if there are points inside the polygon (boolean)
    TFin = isinterior(polyout1,datapc.Location(:,1),datapc.Location(:,2));
    
    %find row numbers for points which are inside the polygon
    rows = find(TFin(:,1)==1);
    [nbInlier,~] = size (rows);
    
    %initialise array for segmented point cloud
    pcdummy=zeros(nbInlier,3);
    
    %transfer the inlier point cloud data into the dummy
    for j=1:nbInlier
        pcdummy(j,:)=datapc.Location(rows(j,1),:);
    end
    
    %transfer to point cloud class
    inPc = pointCloud(pcdummy);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Retrieve the outliers outside the shapefile polygon
    rows4 = find(TFin(:,1)==0);
    [szOut2,~] = size(rows4);
    pcdummy4=zeros(szOut2,3);
    for n=1:szOut2
        pcdummy4(n,:)=datapc.Location(rows4(n,1),:);
    end
    outPc = pointCloud(pcdummy4);
%     figure (100+i)
%     pcshow(outPcseg)
%     title(strcat('Remaining points from',{' '},objectNm));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %optional: denoising of the inlier point cloud
    % inPc = pcdenoise(inPc, 'Threshold', dthres);
    
%     %plot the segmented point clouds
%     figure(i)
%     pcshow(inPc);
%     title(objectNm);
    
    %perform Euclidean segmentation to clean the segmented point cloud
    [labels,~] = pcsegdist(inPc,distThreshold);
    
    %Look for unique labels
    a = unique(labels);
    
    %...create a table with name of label and frequency
    table = [a,histc(labels(:),a)];
    
    %look for the point cloud with the most points after the segmentation
    inSeg = max(table(:,2));
    
    %find the index of this largest point cloud in table
    rowTable = find(table(:,2)==inSeg);
    indexTable = table(rowTable(1,1),1);
    
    %find the index of points with this index in the original point cloud
    rows2 = find(labels(:,1)==indexTable);
    
    %initialize dummy array
    pcdummy2=zeros(inSeg,3);
    
    %transfer the final result into the dummy
    for k=1:inSeg
        pcdummy2(k,:)=inPc.Location(rows2(k,1),:);
    end
    
    %transfer the dummy to point cloud class
    inPcseg = pointCloud(pcdummy2);
    
    %Export segmented point cloud to PLY format
    pcwrite(inPcseg,objectNm,'PLYFormat','binary');
    
    %optional: denoising of the inlier point cloud
%     inPcseg = pcdenoise(inPcseg, 'Threshold', dthres);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Retrieve the outliers from the pcsegdist process
    rows3 = find(labels(:,1)~=indexTable);
    [szOut,~] = size(rows3);
    pcdummy3=zeros(szOut,3);
    for m=1:szOut
        pcdummy3(m,:)=inPc.Location(rows3(m,1),:);
    end
    outPcseg = pointCloud(pcdummy3);
%     figure (100+i)
%     pcshow(outPcseg)
%     title(strcat('Remaining points from',{' '},objectNm));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Plot the clusters individually
    figure('name',objectNm)
    pcshow(inPcseg);
    title(strcat(objectNm,{' '},'after Euclidean segmentation'));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %optional: denoising of the inlier point cloud
%     inPcseg = pcdenoise(inPcseg, 'Threshold', dthres);
%     
%     %Plot the clusters individually
%     figure(20+i)
%     pcshow(inPcseg);
%     title(strcat(objectNm,{' '},'after denoising'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Create the struct containing cluster info
    Object.PtCloud = inPcseg; %cluster's point cloud data
    Object.type = shapes.(objectNm).type; %cluster's type (i.e. shapefile name)
    %fill the struct with fields from the shapefile's attributes (name,
    %length, etc.)
    for l=1:numAtt
        Object.(attList{l}) = shapes.(objectNm).(attList{l});
    end
    
    %Create the struct containing inlier clusters
    Clusters.(objectList{i}) = Object;
    
    %Remerge the remaining unsegmented pointclouds
    datapc=pcmerge(outPc,outPcseg,gridSz);
    
    
    waitbar(i/numObject,f)
end
figure('name','RemainPc')
pcshow(datapc)
title('Remaining point cloud')
close(f);
remainpc=datapc;
