function new_pf = mirror_ref_selection(pf_candidate, n)
% 
% 

[refpoint, mirror] = mirror_direction(n);

if n > size(pf_candidate, 1)
    raise('PF front size is smaller than expected selected points');
end

nadir = max(pf_candidate, [], 1);
ideal = min(pf_candidate, [], 1);

pf_norm = (pf_candidate - ideal) ./ (nadir - ideal);
new_pf = [];



for i = 1 : n
    % cla reset;
    vec1 = pf_norm - mirror(i, :);
    vec2 = refpoint(i, :) - mirror(i, :);
       
    theta = sum(vec1 .* vec2, 2) ./ (sqrt(sum(vec1.^2, 2)) .* sqrt(sum(vec2.^2, 2)));
    theta = acos(theta); % theta is between [0, pi]
    
    % plot([mirror(i, 1), refpoint(i, 1)], [mirror(i, 2), refpoint(i, 2)]); hold on;
    % scatter(pf_norm(:, 1), pf_norm(:, 2));
    
    
    [~, idx] = sort(theta);
    % plot(pf_norm(idx(1), 1), pf_norm(idx(1), 2), 'k*');
    new_pf = [new_pf; pf_norm(idx(1), :)];
    pf_norm(idx(1), :) = [];    
end

new_pf = ideal + (nadir - ideal) .* new_pf;

end

function [refpoint, mirror] = mirror_direction(n)
% This function returns n parallel mirror directions
% for creating uniformly distributed PF points

% (1) generate reference point
refpoint_candidate = 0:n-1;
refpoint_candidate = refpoint_candidate ./ (n-1);

refpoint = [];
for i= 1:n
    coordinate = [refpoint_candidate(i), refpoint_candidate(n-i+1)];
    refpoint = [refpoint; coordinate];
end
mirror = refpoint - 1;

end