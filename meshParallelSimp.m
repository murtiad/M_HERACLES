function [vertices2,faces2] = meshParallelSimp(vertices1, tol)
% MESHPARALLELSIMP
%
% Function to simplify a surface mesh. The function will detect parallels 
% in the shape, and take only the vertices of each of the longest vector in
% the parallel (i.e. the parallel's extremities) as output
%
% Inputs: 
% - vertices1: list of nx3 coordinates of mesh free boundary vertices. The
% mesh should already be a planar surface (cf. roof detection)
% - tol: tolerance for parallelism. The closer the value to 0, the better 
% the original mesh shape is preserved (but the lesser the size reduction)
%
% Outputs:
% - vertices2: list of nx3 coordinates of the simplified mesh
% - faces2: triangular faces of vertices2 resulting froma a Delaunay
% triangulation
%
% (c) Arnadi Murtiyoso (INSA Strasbourg - ICube-TRIO UMR 7357)

% input points
pts=vertices1;
[numPts,~]=size(pts);
indexList=zeros(numPts,1);

for i=1:numPts 
    indexList(i,1)=i;
end

% list all possible vector combinations
VectorList=nchoosek(indexList,2);
[nbCombinations,~]=size(VectorList);

% create a list: columns 1 and 2 are the vector points; colums 3 to 5 the
% xyz non-normalised vector
for i=1:nbCombinations
    index1=VectorList(i,1);
    index2=VectorList(i,2);
    VectorList(i,3)=pts(index1,1)-pts(index2,1);
    VectorList(i,4)=pts(index1,2)-pts(index2,2);
    VectorList(i,5)=pts(index1,3)-pts(index2,3);
end
%% Check for parallel vectors
i=1;
VectorListRef=VectorList;

% keep looping until the VectorList is empty of 'seeds'
while ~isempty(VectorList)
    parallelName = strcat('Parallel',string(i));
    
    % 'seed' vector
    v1 = [VectorList(1,3) VectorList(1,4) VectorList(1,5)];
    [~, index]=ismember(VectorListRef(:,3:5),v1,'rows');
    indexv1=find(index,1);
    
    % normalise the seed vector (v1)
    v1 = normalizeVector3d(v1);
    
    % number of pairs to check the parallelism
    [pairsToCheck,~]=size(VectorList);
    
    % add this seed to the parallel's list of vectors
    listParallel = indexv1;
    ListToDelete=[];
    
    for j=2:pairsToCheck 
        % second vector (to be checked against v1)
        v2 = [VectorList(j,3) VectorList(j,4) VectorList(j,5)];
        [~, index2]=ismember(VectorListRef(:,3:5),v2,'rows');
        indexv2=find(index2,1);
        
        % normalise v2
        v2 = normalizeVector3d(v2);
        
        % check the disparity between the two vectors.
        disparity=abs(v1-v2);
        
        % magDisparity should be 0 if they are parallel
        magDisparity=sqrt(disparity(1)^2+disparity(2)^2+disparity(3)^2);
        
        % set parallelism tolerance here
        if magDisparity < tol
           listParallel=vertcat(listParallel,indexv2);
           ListToDelete = vertcat(ListToDelete,j);
        end
    end
    
    % sanity check. don't delete anything if there is nothing to delete
    if ~isempty(ListToDelete)
        VectorList(ListToDelete,:) = [];
    end
    VectorList(1,:)=[];
    
    % store the vector indices belonging to the same parallel to a struct
    Struct.(parallelName).VectorIndices=listParallel;
    
    % moving on...
    i=i+1;
end
%% compute parallel vector distances
nbFields=length(fieldnames(Struct));
verticesSimp = [];
for i=1:nbFields
    parallelName = strcat('Parallel',string(i));
    [nbVectors,~]=size(Struct.(parallelName).VectorIndices);
    for j=1:nbVectors
        thisVectorIndex=Struct.(parallelName).VectorIndices(j,1);
        v=[VectorListRef(thisVectorIndex,3) VectorListRef(thisVectorIndex,4) VectorListRef(thisVectorIndex,5)];
        distance = sqrt(v(1)^2+v(2)^2+v(3)^2);
       Struct.(parallelName).Distances(j,1)=distance;
    end
    % determine the longest vector in a parallel
    [M,I] = max(Struct.(parallelName).Distances);
    Struct.(parallelName).LongestVecMag = M;
    Struct.(parallelName).LongestVecInd = Struct.(parallelName).VectorIndices(I,1);
    VectorIndex = Struct.(parallelName).VectorIndices(I,1);
    
    % extract the corresponding vertices
    point1Index = VectorListRef(VectorIndex,1);
    point2Index = VectorListRef(VectorIndex,2);
    
    coordPoint1 = pts(point1Index,:);
    coordPoint2 = pts(point2Index,:);
    
    Struct.(parallelName).ExtremePts = vertcat(coordPoint1,coordPoint2);
    dummy = Struct.(parallelName).ExtremePts;
    
    verticesSimp = vertcat(verticesSimp,dummy);
end

% delete doubles
verticesSimp = unique (verticesSimp,'rows');

%% delaunay triangulation
TR = delaunayTriangulation(verticesSimp(:,1),verticesSimp(:,2),verticesSimp(:,3));

[faces2,vertices2] = freeBoundary(TR);
[numPts2,~] = size(vertices2);

gain=(abs(numPts2-numPts)/numPts)*100;

%disp(strcat('Heyo! Vertices reduced by:',num2str(gain),'%'));
end

