function [X_Child] = generate_child_DE(lb, ub,  pop_old, param)
% DE from platemo


NP = param.popsize;
D = size(pop_old, 2);

[CR,F,proM,disM] = deal(1, 0.5, 1, 20);

Lower = repmat(lb, NP, 1);
Upper = repmat(ub, NP, 1);

Parent1 = pop_old;
Parent2 = pop_old(randperm(size(pop_old,1)),:);
Parent3 = pop_old(randperm(size(pop_old,1)),:);

% Differental evolution
Site = rand(NP,D) < CR;
Offspring       = Parent1;
Offspring(Site) = Offspring(Site) + F*(Parent2(Site)-Parent3(Site));


%% Polynomial mutation
Site  = rand(NP,D) < proM/D;
mu    = rand(NP,D);
temp  = Site & mu <= 0.5;
Offspring       = min(max(Offspring,Lower),Upper);
Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
    (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
temp = Site & mu>0.5;
Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
    (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));

X_Child = Offspring;   
end