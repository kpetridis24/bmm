n = 1.5e6;
d = 0.4; % approximate number of true elements per row

% generate random sparse matrices
F = sprand( n, n, d/n ) > 0;
A = sprand( n, n, d/n ) > 0;
B = sprand( n, n, d/n ) > 0;

% add to coo matrices 
[rowF, colF] = find(F);
[rowA, colA] = find(A);
[rowB, colB] = find(B);
cooF = [colF, rowF];
cooA = [colA, rowA];
cooB = [colB, rowB];

% sort coo
cooF = sortrows(cooF,[2 1]);
cooA = sortrows(cooA,[2 1]);
cooB = sortrows(cooB,[2 1]);

% count nnz
nnzF = length(rowF);
nnzA = length(rowA);
nnzB = length(rowB);

% maksed-BMM
tic;
C = F.*(A*B) > 0; 
toc

[rowC, colC] = find(C);
cooC = [colC, rowC];
cooC = sortrows(cooC,[2 1]);
nnzC = length(rowC);

% write data to files
dlmwrite('F.mtx', [n n nnzF] , 'delimiter', ' ', 'precision',  10);
dlmwrite('F.mtx', cooF , '-append', 'delimiter', ' ', 'roffset', 0, 'precision',  10);
dlmwrite('A.mtx', [n n nnzA] , 'delimiter', ' ', 'precision',  10);
dlmwrite('A.mtx', cooA , '-append', 'delimiter', ' ', 'roffset', 0, 'precision',  10);
dlmwrite('B.mtx', [n n nnzB] , 'delimiter', ' ', 'precision',  10);
dlmwrite('B.mtx', cooB , '-append', 'delimiter', ' ', 'roffset', 0, 'precision',  10);
dlmwrite('C.mtx', [n n nnzC] , 'delimiter', ' ', 'precision',  10);
dlmwrite('C.mtx', cooC , '-append', 'delimiter', ' ', 'roffset', 0, 'precision',  10);
