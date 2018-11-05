indexC = 'abce';
indexA = 'aibc';
indexB = 'ei';
rng(0)
A = rand(2,3,4,5) + 1i * rand(2,3,4,5);
B = rand(3) + 1i * rand(3);

%%%
tic
for i = 1:1000
%      C = etprod(indexC,A,indexA,B,indexB);
end
toc
tic
for i = 1:1000
    C_f = fitko_etprod(indexC,A,indexA,B,indexB);
end
toc
check=C-C_f;
sum(check(:))

%%

indexC = 'j';
indexA = 'ji';
indexB = 'i';
rng(0)
A = rand(1,3) + 1i * rand(1,3);
B = rand(3,1) + 1i * rand(3,1);

%%%
tic
for i = 1:1000
    C = etprod(indexC,A,indexA,B,indexB);
end
toc
tic
for i = 1:1000
    C_f = fitko_etprod(indexC,A,indexA,B,indexB);
end
toc
check=C-C_f;
sum(check(:))

%%
indexC = 'abce';
indexA = 'aibcj';
indexB = 'jei';
rng(0)
A = rand(2,3,4,5,7) + 1i * rand(2,3,4,5,7);
B = rand(7,3,3) + 1i * rand(7,3,3);

C_f = fitko_etprod_multi(indexC,A,indexA,B,indexB);
