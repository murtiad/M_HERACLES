function [Building,frontierSliceNb,Slices] = atticsegment (ptCloud)
% ATTICSEGMENT
%
% Function to detect a building's attic (space under the roof). 
%
% Inputs: 
% - ptCloud: point cloud data of a building
%
% Outputs:
% - Building: a struct containing the attic and body point clouds.
% - frontierSliceNb: slice number of the last slice before the attic. 
% - Slices: a struct with the individual slices as fields.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

%%
pas = 0.2; %set the slicing interval
thc = 0.05; %set slice thickness
tol_ground=0.2; %set tolerance for min Z
[Slices,slices_list] = slices(ptCloud,pas,thc,tol_ground);
f = waitbar(0,'Detecting the attic...');
[nbSlice,~] = size(slices_list);
Areas=zeros(nbSlice,2); %set dummy matrix to keep surface areas

%begin iteration
for i=1:nbSlice
    ThisSlice = Slices.(slices_list{i});%name of this slice
    [nRow,~] = size(ThisSlice.Location);
    
    %check Z max an Z min for this slice
    Zmax = max(ThisSlice.Location(:,3));
    Zmin = min(ThisSlice.Location(:,3));
    
    %set an average Z to project the points to
    Z=(Zmax+Zmin)/2;
    
    %project the points to a plane Z
    ProjectedPts = zeros(nRow,3); %initialise matrix
    ProjectedPts(:,1)=ThisSlice.Location(:,1); %x
    ProjectedPts(:,2)=ThisSlice.Location(:,2); %y
    ProjectedPts(:,3)=Z; %projected z
    
    %create a convexhull around the projected points
    [U,S]=svd(   bsxfun(@minus,ProjectedPts,mean(ProjectedPts)),   0);  
    [~,area]=convhull(U*S(:,1:2));
    
    %create a matrix with convex hull points
    %not using it at the moment, but maybe in the future?
    %     [ksize,~] = size(k);
    %     chull_pts = zeros(ksize+1,3);
    %     
    %     for j=1:ksize
    %         chull_pts(j,1)=ProjectedPts(k(j,1),1);
    %         chull_pts(j,2)=ProjectedPts(k(j,1),2);
    %         chull_pts(j,3)=Z;
    %     end
    %     chull_pts(ksize+1,:)=chull_pts(1,:);
    
    %stock the computed area in Areas' 1st column, its Z in 2nd column
    Areas(i,1) = area;
    Areas(i,2) = Z;
    waitbar((i/nbSlice)*0.75,f)
end
%create a table of the areas and show it
% AreaTable = table(Areas(:,1),Areas(:,2));
% AreaTable.Properties.VariableNames={'Area_m2' 'Z'};
% AreaTable

%% Check for area differences between consecutive slices
[nIteration,~] = size(Areas);
for i=1:nIteration-1
    ecart = abs(diff([Areas(i,1),Areas(i+1,1)])); %compute absolute area difference
    pctgEcart=(ecart/Areas(i,1))*100; %area difference in percentage
    tolEcart = 10; %set the tolerance of area difference
    if pctgEcart>tolEcart
        Zlimit = (Areas(i,2)+Areas(i+1,2))/2;
        frontierSliceNb = i;%keep the index of the last slice before the roof
        disp(strcat('[DING!]An attic has been detected at Z=',num2str(Zlimit),'. Hooray!'));
        break %if slice n+1 is different to slice n more than the tolerance, stop the iteration here
    elseif (i==nIteration-1) && (pctgEcart<tolEcart)
        Zlimit = (Areas(nIteration,2));
        frontierSliceNb = nIteration;
        disp('[DING?]No attic detected :/');
    end
end

%%
ptindex_body = find(ptCloud.Location(:,3)<=Zlimit); %find point cloud indexes for building body
ptindex_attic = find(ptCloud.Location(:,3)>Zlimit); %find point cloud indexes for roofs

ptCloud_body = select(ptCloud,ptindex_body); %segment the body point cloud
ptCloud_attic = select(ptCloud,ptindex_attic); %segment the roof point cloud

Building.Body = ptCloud_body;
Building.Attic = ptCloud_attic;

figure('name','Attic Segmentation')
fig=pcshowpair(ptCloud_body,ptCloud_attic); %show the two points clouds
title('Attic Segmentation');
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
saveas(fig,strcat('.\03_Output\99_Figs\attic_segmentation.jpg'));

close(f);