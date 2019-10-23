function [Regions2] = regiongrowingnormalsOct (ptCloud,normals2,h,binCapacity)
% REGIONGROWINGNORMALSOCT
%
% Function to perform region growing on a point cloud, based on computed 
% normal and curvature values. This function is based on octree bins and is
% therefore much faster than REGIONGROWINGNORMALS. The octree division
% operation used the OcTree Matlab toolbox from Sven Holcombe (2013).
%
% See also: 
% - REGIONGROWINGNORMALS
%
% Inputs: 
% - ptCloud: point cloud data 
% - normals2: normals of each point. Compute with pcnormals
% - h: mean curvature of each point. Compute using Beksi (2014)
% - binCapacity: number of points in each octree bin; default is 20
%
% Outputs:
% - Regions2: a struct containing the segmented regions.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

if nargin==0, help(mfilename), return, end

if nargin<4, binCapacity=20; end

labels=zeros(ptCloud.Count,1);
labels(:,1)=1;
OT = OcTree(ptCloud.Location,'binCapacity',binCapacity);  

% % OPTIONAL: plot the octree bins created from the previous operation
% figure ('name','Octree Bins')
% pcshow(ptCloud,'MarkerSize',20)
% hold on
% boxH = OT.plot;
% cols = lines(OT.BinCount);
% dim = [.2 .5 .3 .3];
% for i = 1:OT.BinCount
%      set(boxH(i),'Color',cols(i,:),'LineWidth', 1+OT.BinDepths(i))
% end
% doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});

% create a list of octree bins with points inside
[binsIsi, ~, indexBin] = unique(OT.PointBins);
[r,~]=size(binsIsi);

% initialise variables
% NAvg is a lookup table with column order: binID, median of Nx, median of
% Ny, median of Nz, median of mean curvature (h), and standard deviation of
% mean curvature (h)
NAvg=zeros(r,6);
% binCentroid is the XYZ coordinates of the bins' center
binCentroid=zeros(r,3);

% start loop for each bin
for i=1:r
    % get the current bin's centroid coordinates
    % note: in BinBoundaries, the order of the columns is Xmin, Ymin, Zmin,
    % Xmax, Ymax, Zmax
    binCentroid(i,1)=(OT.BinBoundaries(binsIsi(i),4)+ ...
        OT.BinBoundaries(binsIsi(i),1))/2;
    binCentroid(i,2)=(OT.BinBoundaries(binsIsi(i),5)+ ...
        OT.BinBoundaries(binsIsi(i),2))/2;
    binCentroid(i,3)=(OT.BinBoundaries(binsIsi(i),6)+ ...
        OT.BinBoundaries(binsIsi(i),3))/2;
    
    % get the current bin's ID and put it in NAvg's first column
    NAvg(i,1)=binsIsi(i);
    
    % find the indices of the points located inside the current bin
    indexBin2=find(binsIsi==NAvg(i,1));
    isiBin=find(indexBin==indexBin2);
    
    % loop to get the normals of the points inside the current bin
    [r2,~]=size(isiBin);
    for j=1:r2
        nx(j)=normals2(isiBin(j),1);
        ny(j)=normals2(isiBin(j),2);
        nz(j)=normals2(isiBin(j),3);
        
        % ...and their curvature also, while you're at it
        ho(j)=h(isiBin(j));        
    end
    
    % compute the "normal" of the bin, from the median of the
    % normals of the points inside the bin (get it?)
    nx_median=median(nx);
    ny_median=median(ny);
    nz_median=median(nz);
    
    % ...and the median and standar deviation of the curvature also
    ho_median=median(ho);
    ho_std=std(ho);
    
    % put these respective values in our lookup table
    NAvg(i,2)=nx_median;
    NAvg(i,3)=ny_median;
    NAvg(i,4)=nz_median;
    
    NAvg(i,5)=ho_median;
    NAvg(i,6)=ho_std;
    
end

% create a (sparser) point cloud from the bins' centroids 
binPtCloud=pointCloud(binCentroid);

% create a matrix for the bins' normals too. We consider this as
% binPtCloud's normals
binNormals=NAvg(:,2:4);

% OPTIONAL: draw the bin point cloud and its normal directions
figure ('name','Sparse Bin Point Cloud Normals')
pcshow(binPtCloud,'MarkerSize',20)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on
quiver3(binPtCloud.Location(:,1),binPtCloud.Location(:,2), ...
    binPtCloud.Location(:,3),binNormals(:,1),binNormals(:,2), ...
    binNormals(:,3));

% create a matrix for the bins' curvature too. We consider this as
% binPtCloud's curvature
binCurv=transpose(NAvg(:,5));

% compute the distance between bin centroids. We will use this later to
% determine the nearest neighbour threshold
% here I used the index "i", but it can be anything between 1 and the last
% i value from the last loop (the centroid-to-centroid distance should be
% more or less constant)
binRadius=sqrt((OT.BinBoundaries(binsIsi(i),6)- ...
    OT.BinBoundaries(binsIsi(i),3))^2+(OT.BinBoundaries(binsIsi(i),5)- ...
    OT.BinBoundaries(binsIsi(i),2))^2+(OT.BinBoundaries(binsIsi(i),4)- ...
    OT.BinBoundaries(binsIsi(i),1))^2);

% use the regiongrowingnormals function on our sparse bin point cloud
% notice the default angle of 10 degrees (bigger than the default because
% the point cloud is sparser)
% notice also the min point threshold of 27 (assuming all the bin's
% neighbours are considered as the same region)
Regions=regiongrowingnormals(binPtCloud,binNormals,binCurv,10,1, ...
    binRadius*2,27);

% check the number of regions created from the previous operation
nbRegions=numel(fieldnames(Regions));

% start loop to retrieve the other points
for i=1:nbRegions
    
    % name of current region ("Region1","Region2",etc.)
    thisRegionName=strcat('Region',int2str(i));
    
    % initialise a list of this region's member points id
    thisRegion=[];
    
    % loop to retrive the points in the current region
    for j=1:Regions.(thisRegionName).Count
        % find the id of the current bin in the original octree structure
        [~, index1a]=ismember(binCentroid, ...
            Regions.(thisRegionName).Location(j,:),'rows');
        index1=find(index1a,1);
        index2=binsIsi(index1);
        indexRegion=find(OT.PointBins==index2);
        
        % add the point indices located in the current bin into the region
        % point id list
        thisRegion=[thisRegion;indexRegion];    
    end
    
    % create a point cloud from the points' indices in the current region
    ptCloudOut=select(ptCloud,thisRegion);
    
    % create a label for the region (for visualisation purposes)
    for l=1:ptCloudOut.Count
            labels(thisRegion(l,1),1)=i+1;
    end
    
    % put the results in a struct called Regions2
    Regions2.(thisRegionName)=ptCloudOut;
end

% plot the result
figure('name','Octree Segmented Point Cloud')
pcshow(ptCloud.Location,labels,'MarkerSize',20)
colormap(hsv(i+1))
toc