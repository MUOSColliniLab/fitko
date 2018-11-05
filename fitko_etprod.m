function C = fitko_etprod(indexC, A, indexA, B, indexB)
% Einstein tensor product between n-dimensional arrays. Only one index can
% be multiplied.
%
% Example: C = fitko_etprod('ijk', A, 'ijl', B, 'lk');

dimA = numel(indexA);
dimB = numel(indexB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A

% verifiy the matching of the indexes between A and C
indexAC = zeros(size(indexA));
for i=1:numel(indexA)
    for j=1:numel(indexC)
        if indexA(i) == indexC(j)
            indexAC(i) = j;
        end
    end
end

% sort the indexes and permute A
indexAC(indexAC==0)=numel(indexC)+1;
[indexAC,sort_indexAC] = sort(indexAC);
if numel(sort_indexAC) > 1
    A = permute(A, sort_indexAC);
end

% size of A, without the multiplied index
sizeA = size(A);
if numel(sizeA) < dimA
    sizeA(dimA) = 1;
end
sizeA(dimA) = [];
indexAC(dimA) = [];

% reshape A
ncolA = 1;
for i=1:(dimA-1)
    ncolA = ncolA*size(A,i);
end
A = reshape(A, [ncolA size(A,dimA)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% B

% verifiy the matching of the indexes between B and C
indexBC = zeros(size(indexB));
for i=1:numel(indexB)
    for j=1:numel(indexC)
        if indexB(i) == indexC(j)
            indexBC(i) = j;
        end
    end
end

% sort the indexes and permute B
indexBC(indexBC==0)=numel(indexC)+1;
[indexBC,sort_indexBC] = sort(indexBC);
if numel(sort_indexBC) > 1
    B = permute(B, sort_indexBC);
end

% size of B, without the multiplied index
sizeB = size(B);
if numel(sizeB) < dimB
    sizeB(dimB) = 1;
end
sizeB(dimB) = [];
indexBC(dimB) = [];

% reshape B
ncolB = 1;
for i=1:(dimB-1)
    ncolB = ncolB*size(B,i);
end
B = reshape(B, [ncolB size(B,dimB)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% C 

% matrix multiplication
C = A * B.';

% reconstruct matrix C
C = C(:);
sizeC = [sizeA sizeB];
C = reshape(C, sizeC);

% permute C with the final ordering
[~,sort_indexAC_indexBC] = sort([indexAC indexBC]);
if numel(sort_indexAC_indexBC) > 1
    C = permute(C, sort_indexAC_indexBC);
end


end