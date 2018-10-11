function [Walls,remainPtCloud] = wallSeg(datapc,inc_tol,thres)
% WALLSEG
%
% Function to extract walls from a point cloud of a buildings. The wall
% detection uses normals to determine if they face the same direction or
% not. The output is a struct WALLS, with the number of detected walls as
% number of fields. Each wall object contains the segmented point cloud of
% the wall, the alphashape, and the mesh of the extracted surface. The
% function also creates .ply files for the results, located in the folder
% "03_Outputs".
%
% Inputs: 
% - datapc: point cloud data of a building
% - inc_tol: tolerance for the inclination angle for the wall normals. The
% hypothesis is that a wall has, as the normal's inclination angle, 90
% degrees +/- inc_tol
% - thres: distance threshold for the plane fitting algorithm
%
% Outputs:
% - Walls: a struct with the individual segmented objects as fields.
% Each field consists of another struct ("Object") which keeps the point
% cloud data as well as the alphashape and mesh parameters of the generated
% plane 3D model
% - remainPtCloud: the remaining point cloud after the objects are segmented
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

format long g;
clear Object;
clear Walls;

%Load the point cloud
ptCloud = datapc;

%Compute cloud normals
normals = pcnormals(ptCloud);
[n_row,~]=size(normals);

%Compute horizontal (bearing) and vertical (inclination) angles from
%normals
r=zeros(n_row,2);
for i=1:n_row
    [phi,~] = bearing_surv(0,0,normals(i,1),normals(i,2)); % angle planimetric
    theta = rad2deg(acos(normals(i,3))); % angle inclination
    r(i,1) = phi;
    r(i,2) = theta;
end

%segment the walls. HYPOTHESIS: walls are perpendicular to nadir
rows = find(r(:,2)<90+inc_tol & r(:,2)>90-inc_tol);

%segmented point cloud of ALL the walls
wallsPc = select(ptCloud,rows);

%re-compute normals and normal direction angles for the wallsPc... IMO this
%can be more optimised in the future
normals2 = pcnormals(wallsPc);
[n_row2,~]=size(normals2);

r2=zeros(n_row2,2);
for i=1:n_row2
    [phi,~] = bearing_surv(0,0,normals2(i,1),normals2(i,2)); % angle planimetric
    theta = rad2deg(acos(normals2(i,3))); % angle inclination
    r2(i,1) = phi;
    r2(i,2) = theta;
end

%Plot the histogram of the bearings
figure('Name','BearingsHistogram')
h=histogram(r2(:,1),'BinMethod','fd'); %fd seems to be the best in this case
title('Histogram of Normal Bearings - Spikes may indicate walls')
xlabel('Bearing angle (degrees)')
ylabel('Number of points')

%find local maximas in the histogram. HYPOTHESIS: local maxima of the
%histogram of bearings implies a single wall direction
minProm = ceil(0.05*wallsPc.Count); %set minimal prominence to 5% of total points
[ilm,~]=islocalmax(h.Values,'FlatSelection','center','MinProminence',minProm);
% polarhistogram(deg2rad(r(:,1)),h.NumBins); %optional: plot polarhistogram
iterations_tol = 1;%add number of walls, just in case
iterations=sum(ilm(1,:)~=0)+iterations_tol;

%create empty list for object names
mur_list = strings(iterations,1);

%begin the iterations
for i=1:iterations
    %create list of object names
    mur_list(i,1) = strcat('WALL',string(i));
    
    %fit a plane to the wall point cloud, and then segment the individual
    %wall point clouds
    [model,inlierIndices,outlierIndices] = pcfitplane(wallsPc,thres);
    plane = select(wallsPc,inlierIndices);
    remainPtCloud = select(wallsPc,outlierIndices);
    wallsPc = remainPtCloud;
    
    %convert the Matlab surface model from pcfitplane to geom3d format
    a=model.Parameters(1,1);
    b=model.Parameters(1,2);
    c=model.Parameters(1,3);
    d=model.Parameters(1,4);

    P1 = [plane.Location(1,1) plane.Location(1,2) 0];

    P1(1,3) = (-d-(a*P1(1,1))-(b*P1(1,2)))/c;

    planeGeom3d = createPlane(P1, model.Normal); %plane in geom3d format
    
    %Perform 3D Delaunay Triangulation on the wall point cloud
    TR = delaunayTriangulation(double(plane.Location(:,1)),double(plane.Location(:,2)),double(plane.Location(:,3)));
    
    %extract the surface of the Delaunay mesh TR (otherwise stored in
    %tetrahedral form)
    [~,xf] = freeBoundary(TR);
    [size_xf,~] = size(xf);
    
    %project the points to the mathematical surface model to get 1 surface
    %in the form of the wall
    prjtdPts = zeros(size_xf,3);
    
    for j=1:size_xf
        %create a 3D ray for each foint in the mesh surface
        line = createLine3d(xf(j,1), xf(j,2), xf(j,3), model.Normal(1,1), model.Normal(1,2), model.Normal(1,3));
        
        %intersect the 3D ray with the mathematical surface
        prjtdPt = intersectLinePlane(line, planeGeom3d);
        
        %store the coordinates of the intersection
        prjtdPts(j,1) = prjtdPt(1,1);
        prjtdPts(j,2) = prjtdPt(1,2);
        prjtdPts(j,3) = prjtdPt(1,3);
    end
    
    %perform 3D Delaunay Triangulation on the projected points
    TR2 = delaunayTriangulation(double(prjtdPts(:,1)),double(prjtdPts(:,2)),double(prjtdPts(:,3)));
    
    %extract the surface of the Delaunay mesh of the projected points
    [F,xf2] = freeBoundary(TR2);
    
    %further smoothing here?
    
    %store the mesh properties in a struct
    Mesh.Vertices = xf2;
    Mesh.Faces = F;
    
    %create a Matlab alphashape
    shp = alphaShape(xf2(:,1),xf2(:,2),xf2(:,3),inf);

    %store all of this stuff in the another struct called Object
    Object.PtCloud = plane;
    Object.AlphaShp = shp;
    Object.Mesh = Mesh;
    
    %store the objects in a struct containing all the detected walls
    Walls.(mur_list{i}) = Object;
    
    %write a ply file of the segmented wall point cloud
    pcwrite(plane,strcat('.\03_Output\01_PtCloud\',mur_list(i),'_step3_PtCloud.ply'));

    %write a ply file of the wall's planemesh
    fname = convertStringsToChars(strcat('.\03_Output\02_Meshes\',mur_list(i),'_step3_planeMesh.ply'));
    writeMesh_ply(fname, Mesh.Vertices, Mesh.Faces);
end

%plot the segmented walls
figure ('Name','DetectedWalls')
axis equal
title('Detected Walls Plane Mesh')
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
for i=1:iterations
    drawMesh(Walls.(mur_list(i,1)).Mesh.Vertices,Walls.(mur_list(i,1)).Mesh.Faces);
    hold on
end
