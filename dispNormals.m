function [] = dispNormals(path,step)
% DISPNORMALS
%
% Function to quickly display normal directions in the form of arrows
%
% Inputs: 
% - ptCloud: point cloud data (.PLY)
% - step: the step between points in the point cloud to show the normals
% (default value is 1 point normal for every 10 XYZ points)
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

if nargin<2
  step = 10;
end

ptCloud=pcread(path);

normals = ptCloud.Normal;

x = ptCloud.Location(1:step:end,1);
y = ptCloud.Location(1:step:end,2);
z = ptCloud.Location(1:step:end,3);
u = normals(1:step:end,1);
v = normals(1:step:end,2);
w = normals(1:step:end,3);

figure ('name','Point Cloud Normals')
title('Point Cloud Normals');
pcshow(ptCloud)
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
hold on
h=quiver3(x,y,z,u,v,w);
set(h,'AutoScale','on','AutoScaleFactor', 2);