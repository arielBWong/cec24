function [result] = ranksum_sigtest_smallerBetter(metric_base, metric_new)
% argument1: baseline
% argument2: new method
% check which is smaller
% result -1: baseline is better
% result 1: new method is better
% result 0: no difference
% 1 means  the former is smaller,  second is bigger,

[p1,h1,stats1] = ranksum( metric_base, metric_new,  'alpha', 0.05, 'tail', 'left');
[p2,h2,stats2] = ranksum( metric_new, metric_base,  'alpha', 0.05, 'tail', 'left');
if h1 == 1 && h2 == 0
    result = -1;
elseif h2 == 1 && h1 == 0
    result = 1;
else
    result = 0;
end

end