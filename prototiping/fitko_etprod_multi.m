function C = fitko_etprod_multi(indexC, A, indexA, B, indexB)
% Einstein tensor product between n-dimensional arrays. Only one index can
% be multiplied.
%
% Example: C = fitko_etprod('ijk', A, 'ijl', B, 'lk');

dimA = numel(indexA);
dimB = numel(indexB);
dimC = numel(indexC);

% verifiy the matching of the indexes 
indexAC = compare_index(indexA,indexC);
indexBC = compare_index(indexB,indexC);
indexAB = compare_index(indexA,indexB);

if numel(indexAC==0) ~= numel(indexAB==0)
    error('Multiplied indexes do not match.')
end

% number of multiplied indexes
N = sum(indexAC==0);

% assign index
indexAC(indexAC==0) = (1:sum(indexAB~=0)) + dimC;
indexBC(indexAB(indexAB~=0)) = indexAC(indexAB~=0);


% cosa fare:
% mettere la dimensione più grande per ultima
%
% calcolarsi tutte le combinazioni degli indici secondo, terzo, ecc...
% [ca, cb, cc] = ndgrid(A, B, C);
% combs = [ca(:), cb(:), cc(:)]
%
% finire in un solo ciclo for accumulando matrici C, utilizzare squeeze per
% rimuovere le dimensioni singleton













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A

% sort the indexes and permute A
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

function indexXC = compare_index(indexX,indexC)
    % verifiy the matching of the indexes between X and C
    indexXC = zeros(size(indexX));
    for i=1:numel(indexX)
        for j=1:numel(indexC)
            if indexX(i) == indexC(j)
                indexXC(i) = j;
            end
        end
    end
end