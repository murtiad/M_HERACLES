function setup()
%SETUP Setup paths, etc., for M_HERACLES

% Get dir of executing (i.e. this!) file.
baseDir=fileparts(mfilename('fullpath'));

% Add selected subdirectories.
addpath(fullfile(baseDir,'02_ThirdParty','euclideanDistanceTwoPointClouds(Cheong2016)'),'-end')
%addpath(fullfile(baseDir,'02_ThirdParty','findPointNormals(Taylor2015)'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','CSF(Zhang2016)','matlab'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','m_map1.4'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','OcTree(Holcombe2013)'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','surf2stl'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','CuboidFittingRanSAC(Mehmood2017)'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','estimate_curvatures(Beksi2014)'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom3d(Legland2005)','geom3D'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom3d(Legland2005)','meshes3D'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom2d(Legland2005)','geom2D'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom2d(Legland2005)','utils'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom2d(Legland2005)','polygons2d'),'-end')
addpath(fullfile(baseDir,'02_ThirdParty','geom2d(Legland2005)','polynomialCurves2d'),'-end')
addpath(baseDir,'-end')

fprintf('You can now use M_HERACLES from everywhere.\n')
