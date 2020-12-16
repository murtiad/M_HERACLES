function [recenteredPtCloud] = ptCloudRecenter(ptCloud,translation)
% PTCLOUDRECENTER
%
% Funtion to recentre point clouds previously centered using PTCLOUDCENTER.
% This function returns the centred point cloud to its original coordinate 
% system.
%
% Inputs: 
% - ptCloud: point cloud data 
% - translation: translation vector from PTCLOUDCENTER
%
% Outputs:
% - recenteredPtCloud: point cloud with recentred coordinates.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

ptRGB=ptCloud.Color;
ptNormals=ptCloud.Normal;
ptIntensity=ptCloud.Intensity;

x=ptCloud.Location(:,1)+translation(1);
y=ptCloud.Location(:,2)+translation(2);
z=ptCloud.Location(:,3)+translation(3);

recenteredPtCloud=pointCloud([x y z],'Color',ptRGB,'Normal',ptNormals,...
    'Intensity',ptIntensity);
end

