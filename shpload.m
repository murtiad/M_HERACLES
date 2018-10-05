function [Shape,objectList, fieldNm] = shpload(filepath)
% SHPLOAD
%
% Loads a polygonal ESRI .shp file and convert it into a struct called 'Shape'.          
% 'Shape' will contain as fields the individual objects in the file.    
% Each field will contain a struct with information available in the .dbf
% file attached to the .shp file, as well as object type and geometry, both
% in the form of a list of vertex coordinates and in native Matlab
% polyshape object type.
%
% This function uses the M_MAP toolbox developed by Rich Pawlowicz (Univ.
% of British Columbia)
% https://www.eoas.ubc.ca/~rich/map.html
%
% Input : 
% - filepath of the shapefile
%
% Output :
% - Shape : struct with individual objects as fields
% - objectList : array with a list of objects/records
% - fieldNm : array with a list of attributes in the shapefile
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)


clear Shape
clear Object

%read the shapefile
M=m_shaperead(filepath);

%extract the file name, i.e. object type (buildings,trees, etc)
[~,fname,~] = fileparts(filepath);

%determine the number of objects and the number of fields in the dbf
[numObject,~] =size (M.ncst);
[~,numField] =size (M.fieldnames);

f = waitbar(0,'Converting shapefile to native structure....');

%create dummy for struct field names
stdummy=strings(numObject,1);

%extract field names
fieldNm=transpose(M.fieldnames);

for i=1:numObject
    waitbar(i/numObject,f)
    
    %object type is its filename (e.g. buildings, trees, etc)
    Object.type=fname;
    
    %create dynamic field for the struct Object, from the list of fields
    %present in the .dbf file
    for j=1:numField
        Object.(fieldNm{j}) = M.dbfdata(i,j);
    end
    
    %extract the geometry
    Object.geometry = cell2mat(M.ncst(i,1));
    
    %create a polyshape object and add it to the struct
    Object.polyshape = polyshape(Object.geometry(:,1),Object.geometry(:,2));

    %Fill the string dummy with Shape struct field names
    stdummy(i,1) = strcat(fname,num2str(i));
   
    %Plot the objects individually
%     figure(i)
%     plot(Object.polyshape);
%     title(strcat(fname,{' '}, num2str(i)));
%     pbaspect([1 1 1]);
%     daspect([1 1 1]);
    
    %Create the struct containing the objects
    Shape.(stdummy{i}) = Object;
        
end

close(f);

%create a list of fields for Shape, ergo list of objects. Might be useful
%if to be used in loops
objectList = stdummy;
