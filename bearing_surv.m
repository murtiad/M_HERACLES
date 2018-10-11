function [alpha, distance] = bearing_surv(Xa,Ya,Xb,Yb)
% BEARING_SURV
%
% Function to compute bearing angle in the surveying sense, i.e. horizontal
% angle measured from the PROJECTION north. Note that this bearing does not
% give azimuth, which is calculated differently and refers to MAGNETIC
% north. 
% Edit: the function has been extended to compute also the 2D distance
% between the two input points, therefore resolving the first geodetic
% problem
%
% Inputs:
% - Xa, Ya : Cartesian coordinates of the first (origin) point
% - Xb, Yb : Cartesian coordinates of the second (target) point
%
% Outputs:
% - alpha : bearing angle in degrees from point A to point B
% - distance : 2D distance from point A to point B
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

delta_x = Xb-Xa;
delta_y = Yb-Ya;

if delta_x>0 && delta_y>0
    alpha = rad2deg(atan(delta_x/delta_y));
elseif delta_x>0 && delta_y<0
    alpha = 180-abs(rad2deg(atan(delta_x/delta_y)));
elseif delta_x<0 && delta_y<0
    alpha = 180+abs(rad2deg(atan(delta_x/delta_y)));
elseif delta_x<0 && delta_y>0
    alpha = 360-abs(rad2deg(atan(delta_x/delta_y)));
elseif delta_x==0 && delta_y>0
    alpha = 0;
elseif delta_x==0 && delta_y<0
    alpha = 180;
elseif delta_x>0 && delta_y==0
    alpha = 90;
elseif delta_x<0 && delta_y==0
    alpha = 270;
else
    disp('No bearing! Are you sure with your coordinates?');
end

distance = sqrt(delta_x^2+delta_y^2);