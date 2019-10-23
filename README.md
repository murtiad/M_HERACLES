# M_HERACLES
HERitAge by point CLoud procESsing for Matlab

(c) by Arnadi Murtiyoso                          
Photogrammetry and Geomatics Group, ICube UMR 7357 INSA Strasbourg
Contact: arnadi.murtiyoso@insa-strasbourg.fr
https://github.com/murtiad          

A toolbox with functions for processing point cloud data in the context of cultural heritage documentation.

![GIS-based segmentation](https://github.com/murtiad/M_HERACLES/blob/master/test/GitHub_images/complextobuildings.JPG)

![Architectural element segmentation1](https://github.com/murtiad/M_HERACLES/blob/master/test/GitHub_images/segmentation_classif.JPG)

![Architectural element segmentation2](https://github.com/murtiad/M_HERACLES/blob/master/test/GitHub_images/gif_classif_structures.gif)

![Region growing](https://github.com/murtiad/M_HERACLES/blob/master/test/GitHub_images/regiongrowing.JPG)

The code was developped with the Matlab Computer Vision Toolbox installed (2018a), as well as third party dependencies:
- Cloth Simulation Filter (CSF) by Wuming Zhang et al. (Zhang W, Qi J, Wan P, Wang H, Xie D, Wang X, Yan G. An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth Simulation. Remote Sensing. 2016; 8(6):501.): https://github.com/jianboqi/CSF
- M_MAP Toolbox by Rich Pawlowicz: https://www.eoas.ubc.ca/~rich/map.html. Note that I used this because I don't have the Matlab Mapping Toolbox installed :(
- geom3d and geom2d from the matGeom toolbox by David Legland: https://github.com/mattools/matGeom. I modified the drawing function slightly to be able to visualise large topographic coordinates (included in the folder 02_ThirdParty)
- William Beksi's function to compute principal curvatures, available on Matlab File Exchange: https://fr.mathworks.com/matlabcentral/fileexchange/46772-estimate-principal-curvatures?s_tid=prof_contriblnk
- Sven Holcombe's OcTree toolbox, available on Matlab File Exchange: https://fr.mathworks.com/matlabcentral/fileexchange/40732-octree-partitioning-3d-points-into-spatial-subvolumes

Available functions (23/10/2019):
- shapeseg.m : function to use (polygonal) ESRI shapefiles (.shp) to delimit ("cookie cutter" style) 3D point cloud area to be then cleaned using the pcsegdist function. This function generates separate segmented point clouds for each object, with their attributes as per the description in the shapefile.
- clustering.m : function to create individual point clouds in Matlab structure from an original point cloud segmented using the pcsegdist function 
- shpload.m : Loads a polygonal ESRI .shp file and convert it into a struct called 'Shape'. 'Shape' will contain as fields the individual objects in the file. Each field will contain a struct with information available in the .dbf file attached to the .shp file, as well as object type and geometry, both in the form of a list of vertex coordinates and in native Matlab polyshape object type. 
- wallSeg.m : automatic 3D modelling function to detect walls in a building point cloud, then segment them, fit a 3D plane, and generate the (therefore simplified) 3D surface of the walls. In development
- slices.m : function to create multiple vertical slices of a given point cloud. 
- atticsegment.m : function to separate a building's attic (space under the roof) from its body. 
- supportdetect.m : function to detect structural supports (columns or piers) from a building body.
- regiongrowingnormals.m : an implementation of PCL's greedy region growing method, based on normal angles and curvatures. See also theoretical primer as explained in http://pointclouds.org/documentation/tutorials/region_growing_segmentation.php
- regiongrowingnormalsOct.m : region growing based on smoothness constraints (normals and angles) that work much faster by implementing octree sub-divisions
