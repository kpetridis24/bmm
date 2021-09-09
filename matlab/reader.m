clear;
format long g;

%% select matrix
% fid = fopen('com-Youtube.mtx');
% n = 1134890;
% nnz = 2987624;

% fid = fopen('belgium_osm.mtx');
% n = 1441295;
% nnz = 1549970;

% fid = fopen('s6.mtx');
% nnz = 12;
% n = 6;

fid = fopen('s12.mtx');
nnz = 19;
n = 12;

%% read
e = textscan(fid,'%f %f', nnz);
data = cell2mat(e);
fclose(fid);

S = sparse(data(:, 2), data(:, 1), 1, n, n, nnz);


%% masked bmm
tic;
C = (S.*(S*S)) > 0;
toc
[row, col] = find(C);
cooC = [col, row];
cooC = sortrows(cooC,[2 1]);
% dlmwrite('bmm_res_belgium_osm.mtx', cooC , 'delimiter', ' ', 'precision',  10);


