function [Clusters,slice_list] = slices(ptCloud,pas,thc,tol_ground)
% SLICES
%
% Function to create multiple vertical slices of a given point cloud. 
%
% Inputs: 
% - ptCloud: point cloud data
% - pas: set the slicing interval
% - thc: set the thickness of each slice
% - tol_ground: set a "buffer" tolerance for minimum ground Z. Set to 0 if
% you are sure with your floor/ground position
%
% Outputs:
% - Clusters: a struct with the individual slices as fields.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

f = waitbar(0,'Creating slices...');
Zmax = max(ptCloud.Location(:,3)); %maximum altitude of point cloud
Zmin = min(ptCloud.Location(:,3)); %minimum altitude of point cloud

nSlices = floor((Zmax-Zmin)/pas);

%create empty list for object names
slice_list = strings(nSlices,1);

tol = thc/2; %set slice thickness in one vertical direction

Z=Zmin+tol_ground;
for i=1:nSlices
    %create list of object names
    slice_list(i,1) = strcat('SLICE',string(i));
    
    ptindex_desc = find(ptCloud.Location(:,3)<Z+tol);
    ptindex_asc = find(ptCloud.Location(:,3)>Z-tol);
    
    ptindex_int = intersect(ptindex_desc,ptindex_asc); %index of points within the slice
    
    ptCloudSlice = select (ptCloud,ptindex_int);

    %write a ply file of the segmented wall point cloud
    pcwrite(ptCloudSlice,strcat('.\03_Output\98_temp\',slice_list(i),'.ply'));
    
    Slices.(slice_list{i}) = ptCloudSlice;
    Z = Z+pas;
    waitbar((i/nSlices),f)
end

disp(strcat('[DING!]',num2str(nSlices),' vertical slices were created with an interval=',num2str(pas),' and thickness=',num2str(thc),' .'));

% %plot the slices
% figure ('Name','Slices')
% for i=1:nSlices
%     fig=pcshow(Slices.(slice_list(i,1)));
%     hold on
% end
% axis equal
% title(strcat('Slices with interval=',num2str(pas),' and thickness=',num2str(thc)))
% xlabel('X (m)')
% ylabel('Y (m)')
% zlabel('Z (m)')
% saveas(fig,strcat('.\03_Output\99_Figs\attic_segmentation.jpg'));

close(f);

Clusters =Slices;
