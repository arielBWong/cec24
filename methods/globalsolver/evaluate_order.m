function [pop, archive, front_idx]= evaluate_order(pop, archive, funh_obj, funh_con, cx, gen, param, varargin)
% This function evaluate given x population, add to pop (N to 2N)
% and sort pop with nd and feasibility consideration
% input
%           pop : population of previous generation
%           archive : record of all evolutionary process
%           funh_obj : function handle of objective function
%           funh_con: function handle of constraints
%           cx: child population generated from last step
%           gen: current generation
%           param: evolution parameter (gen popsize)
% output
%           pop: extended and sorted population
%           archive: extended archive
%-------------------------
child.X=cx;
% This is NSGA-II
output = funh_obj(child.X);

if isstruct(output)
    child.F = output.f;
    child.A = output.addon;
    child.Mdl = output.mdl;
    child.trgdata = output.trgdata;
else
    child.F = output;
    child.A = [];
    child.Mdl = {};
    child.trgdata = {};
end


child.C = funh_con(child.X);

archive.sols=[archive.sols;[repmat(gen,size(child.X, 1),1), child.X, child.F, child.C, child.A]];

% Appending X F C to pop
pop.X = [pop.X; child.X];
pop.F = [pop.F; child.F];
pop.C = [pop.C; child.C];
pop.A = [pop.A; child.A];
pop.Mdl = [pop.Mdl, child.Mdl]; % deliminator is ,
pop.trgdata = [pop.trgdata, child.trgdata];

[pop, front_idx] = pop_sort(pop);




    

end

