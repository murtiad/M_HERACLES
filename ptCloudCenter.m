function [centeredPtCloud,translation] = ptCloudCenter(ptCloud,X0,Y0,Z0)
% PTCLOUDCENTER
%
% Funtion to center point clouds with large coordinates (typically coor-
% dinates from projected systems)
%
% Inputs: 
% - ptCloud: point cloud data 
% - X0,Y0,Z0: translation values. If these values are not available, the
% function will use the barycentre instead to translate the point cloud.
%
% Outputs:
% - recenteredPtCloud: point cloud with recentred coordinates.
% - translation: XYZ vector with the translation values.
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

sensorCenter=[mean(ptCloud.Location(:,1)),mean(ptCloud.Location(:,2)), ...
     mean(ptCloud.Location(:,3))];
ptRGB=ptCloud.Color;
ptNormals=ptCloud.Normal;
ptIntensity=ptCloud.Intensity;

if nargin<2
    translation(1)=sensorCenter(1);
    translation(2)=sensorCenter(2);
    translation(3)=sensorCenter(3);
else
    translation(1)=X0;
    translation(2)=Y0;
    translation(3)=Z0;
end

x=ptCloud.Location(:,1)-translation(1);
y=ptCloud.Location(:,2)-translation(2);
z=ptCloud.Location(:,3)-translation(3);

centeredPtCloud=pointCloud([x y z],'Color',ptRGB,'Normal',ptNormals,...
    'Intensity',ptIntensity);
end

