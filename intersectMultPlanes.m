function [vectorXYZ,residuals] = intersectMultPlanes(Planes)
% INTERSECTMULTPLANES
%
% Use classical non-weighted least squares to compute the most probable
% intersection point of multiple 3D planes.
%
% Inputs: 
% - Planes: an Nx4 matrix containing the each plane's parameters (a,b,c,d).
% N is the number of available planes
%
% Outputs:
% - vectorXYZ: a vector containing the most likely 3D coordinates of the
% intersection
% - residuals: an Nx1 matrix containing the residual distances of vectorXYZ
% to the respective planes
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

A=Planes(:,1:3);
F=Planes(:,4)*-1;

X_LS=((A'*A)^-1)*(A'*F);
vectorXYZ=X_LS';

residuals=(A*X_LS)-F;
end

