%% development test on gsolver
% comparision is conducted with matlab ga in unitest_comparewith_gasolver
clearvars;
workdir = pwd;
% idcs = strfind(workdir, '\');
% upperfolder = workdir(1: idcs(end)-1);
% problem_folder = strcat(upperfolder,'\problems\SMD');
% addpath(problem_folder);

pop.F = [3, 3; 
    2, 5; 
    1, 2;
    1, 3; 
    2, 2; 
    2, 1; 
    1, 1; 
    2, 3;
    4, 5;
    7, 1;
    3, 1];
pop.C = [
    -1, -2;
    -3, -4;
    -5,-1;
    -1, -3;
    -4, -1;
    0, 0;
    1, 1;
    2, -1;
    -1, 1;
    3, 4;
    5, -2];

pop = pop_sort(pop);
    