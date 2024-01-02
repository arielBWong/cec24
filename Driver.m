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


% setpref('Internet','SMTP_Server','smtp-mail.outlook.com');
% setpref('Internet','E_mail','bing.wang@adfa.edu.au');
% sendmail('ariel.bingwang@outlook.com','Hello From MATLAB!', ...
%     'Thanks for using sendmail.');

problems = {'DS5(5,4)'};
np = length(problems);
seedmax = 21;
choices = {true, false};
nc = length(choices);

strc_id = 1;
for i = 1 : np
    for j = 1 : seedmax
        for k = 1 : nc
            onerun_parameters(strc_id). problem_str = prob_str{i};
            onerun_parameters(strc_id). seed = j;
            onerun_parameters(strc_id). strategy = choices{k}; % proposed  method
            strc_id = strc_id + 1;
        end
    end
end

nrun = length(onerun_parameters);
parfor ii = 1:nrun
    nested_EA_forBLOP(onerun_parameters(ii).seed,...
        onerun_parameters(ii).problem_str, ...
        onerun_parameters(ii).strategy);
end
