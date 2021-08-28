clear;
format long g;

%% com-Youtube
fid = fopen('com-Youtube.mtx');
e = textscan(fid,'%f %f', 2987624);
data = cell2mat(e);
fclose(fid);

S = sparse(data(:, 2), data(:, 1), 1, 1134890, 1134890, 2987624);

%% s6
% fid = fopen('s6.mtx');
% e = textscan(fid,'%f %f', 12);
% data = cell2mat(e);
% fclose(fid);
% 
% S = sparse(data(:, 2), data(:, 1), 1, 6, 6, 12);

%% s12
% fid = fopen('s12.mtx');
% e = textscan(fid,'%f %f', 19);
% data = cell2mat(e);
% fclose(fid);
% 
% S = sparse(data(:, 2), data(:, 1), 1, 12, 12, 19);

%% COO
C = S.*(S*S) > 0;
[row, col] = find(C);
cooC = [col, row];
cooC = sortrows(cooC,[2 1]);
dlmwrite('bmm_res_com-Youtube.mtx', cooC , 'delimiter', ' ', 'precision',  10);


