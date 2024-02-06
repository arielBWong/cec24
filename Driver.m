clear all;
clc;

mail = 'ariel.bingwang@outlook.com';
password = 'Carol1984123'; 
server = 'smtp-mail.outlook.com';
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port', '587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');
setpref('Internet','E_mail', mail);
setpref('Internet','SMTP_Server', server);
setpref('Internet','SMTP_Username', mail);
setpref('Internet','SMTP_Password', password);

relative_path = fullfile(pwd, 'problems');
addpath(relative_path);
relative_path = fullfile(pwd, 'methods', 'globalsolver');
addpath(relative_path);
relative_path = fullfile(pwd, 'methods', 'globalsolver', 'ND_Sort');
addpath(relative_path);   
relative_path = fullfile(pwd, 'methods', 'causal_relation');
addpath(relative_path);
relative_path = fullfile(pwd, 'postprocess');
addpath(relative_path);

relative_path = fullfile(pwd, 'methods');
addpath(relative_path);

nested_EA_forBLOP(1, 'DS4D(3,2)', 3);

problems = {'DS4(5,4)', 'DS5(5,4)', 'DS4D(5, 4)', 'DS5D(5, 4)', 'DS4(4, 3)', 'DS5(4, 3)', 'DS4D(4, 3)', 'DS5D(4, 3)', 'DS4(3,2)', 'DS5(3,2)', 'DS4D(3, 2)', 'DS5D(3, 2)',};
np = length(problems);
seedmax = 21;
choices = {0, 3};
nc = length(choices);

strc_id = 1;
for i = 1 : np
    for j = 1 : seedmax
        for k = 1 : nc
            onerun_parameters(strc_id). problem_str = problems{i};
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





setpref('Internet','E_mail','ariel.bingwang@outlook.com');
sendmail('bing.wang@adfa.edu.au','Your code has finished running!', ...
    'As above');