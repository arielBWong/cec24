clear all;
clc;

relative_path = fullfile(pwd, 'problems');
addpath(relative_path);
relative_path = fullfile(pwd,'methods', 'globalsolver');
addpath(relative_path);
relative_path = fullfile(pwd,'methods', 'globalsolver', 'ND_Sort');
addpath(relative_path);
relative_path = fullfile(pwd,'methods', 'causal_relation');
addpath(relative_path);

nested_EA_forBLOP();