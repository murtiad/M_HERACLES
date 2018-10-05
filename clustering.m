function [Clusters] = clustering(labels,pointcloud,inlierNb)
% CLUSTERING
%
% Function to create individual point clouds in Matlab          
% structure from an original point cloud segmented using the    
% pcsegdist function 
%
% Input : 
% - labels array (output from pcsegdist)
% - the original pre-segmented point cloud
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

clear Clusters

tic
%Look for unique labels
a = unique(labels);
[numClusters,~] =size (a); 
%...create a table with name of label and frequency
table = [a,histc(labels(:),a)];

%filtering for outliers
for i=1:numClusters
    numPts = table (i,2);
    if numPts < inlierNb
        %If it's outlier, fill with zeros
        table(i,:)=[0,0];
    end
end

%Delete the zeros rows
table( ~any(table,2), : ) = [];
[numrowTable,~] =size (table); 
disp('Lookup table created successfully');
toc

tic
f = waitbar(0,'Creating individual point cloud clusters...');
%create dummy for struct field names
stdummy=strings(numrowTable,1);
for i=1:numrowTable
    waitbar(i/numrowTable,f)
    %create dummy for point cloud data
    pcdummy=zeros(table(i,2),3);
    %find row indexes for label in table(i,1)
    rows = find(labels(:,1)==table(i,1));
    %transfer the point cloud coordinates to the dummy
    [numPts2,~]=size(pcdummy);
    for j=1:numPts2
        pcdummy(j,:)=pointcloud.Location(rows(j,1),:);
    end
    %convert pcdummy to point cloud class
    PtCloud=pointCloud(pcdummy);
    %Fill the string dummy with struct field names
    stdummy(i,1) = strcat('PtCloud',num2str(table(i,1)));
    %Plot the cluster individually
    figure(i)
    pcshow(PtCloud);
    title(strcat('Cluster',{' '}, num2str(table(i,1))));
    %Create the struct containing inlier clusters
    Clusters.(stdummy{i}) = PtCloud;
end
close(f);
disp('Clusters created successfully');
toc
