function [pop,archive] = initialize_pop(funh_obj, funh_con, num_xvar, xl, xu, initmatrix, param)
% funh_obj, funh_con, num_xvar, ops, evals, xu, xl
% this method generate the initial population of global solver
%
%   input
%       funh_obj                                              : function handle for objective function
%       funh_con                                              : function handle to contraints
%       num_xvar                                              : number of design variables
%       xu                                                    : 1d row vector, upper bound of design variables
%       xl                                                    : 1d row vector, lower bound of design variables
%       initmatrix                                            :  partial population inserted to evolution, put in [] if no partial
%                                                                          population is given
%       param                                                 : evolutionary parameters (popsize, gen)
%   output
%       pop                                                   : population structure with fields: X(design variable),
%                                                                                    F(objective), C(constraint)
%       archive                                               : record  all solutions encountered over evolution
%                                                                                      * usage under development
%---------------------------------------------------------------------------------
N = param.popsize;

% create initial population,  w.r.t.  given partial population intimatrix
n_init = size(initmatrix, 1);
n_rest = N - n_init;
% X_pop  = repmat(xl, n_rest, 1) + repmat(xu - xl, n_rest, 1) .* lhsdesign(n_rest, num_xvar);
X_pop =  unifrnd(repmat(xl,n_rest,1),repmat(xu,n_rest,1));
X_pop  = [X_pop; initmatrix];

% make sure initialization unique
X_pop = unique(X_pop, 'rows','stable');


% objective and constraints
output  = funh_obj(X_pop);

if isstruct(output) % customized 
    F_pop   = output.f;
    A_pop   = output.addon;
    Mdl_pop = output.mdl;
    trg_pop = output.trgdata;
else
    F_pop   = output;
    A_pop   = [];
    Mdl_pop = {};
    trg_pop = {};
end

C_pop     = funh_con(X_pop);
% [~,ids,~] = nd_sort(F_pop, (1:size(F_pop,1))');
ids = 1: size(X_pop, 1);
% Storing relevant information pop and archive
pop.X = X_pop(ids,:);
pop.F = F_pop(ids,:);

if isempty(A_pop)
    pop.A = [];
else
    pop.A = A_pop(ids, :);
end

%---dictionary list

if isempty(Mdl_pop)
    pop.Mdl = {};
else
    pop.Mdl = Mdl_pop(ids);
end

if isempty(trg_pop)
    pop.trgdata = {};
else
    pop.trgdata = trg_pop(ids);
end


if isempty(C_pop)
    pop.C = [];
else
    pop.C = C_pop(ids, :);
end

% archive does not save dictionary type component
archive.sols=[repmat(0,N,1),  pop.X, pop.F, pop.C, pop.A];

return

