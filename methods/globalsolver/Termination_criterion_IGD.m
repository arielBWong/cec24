function [termination_flag,termination_vector]  = Termination_criterion_IGD(archive, sliding_window, threshold, gen, termination_vector, verbose)
% Archive is a list of solutions array


termination_flag = false;

% IGD stopping condition is not changing 
% fill in termination vector
vector_id = mod(gen, size(termination_vector, 2)); % [3, 10]
if vector_id == 0  vector_id = 10; end

NDF_t = archive(end).FU; 
NDF_t_1 = archive(end-1).FU;
[deltaZ_ideal, deltaZ_nadir, igd_phi] = IGD_Z_termination(cat(1, NDF_t{:}), cat(1, NDF_t_1{:}));
termination_vector(1, vector_id) = deltaZ_ideal;
termination_vector(2, vector_id) = deltaZ_nadir;
termination_vector(3, vector_id) = igd_phi;


if gen <= sliding_window
    termination_flag = false;
    return
else
    if verbose
        fprintf('Termination condition : %0.3f, %0.3f, %0.3f, compared to %0.3f, gen %d \n',max(termination_vector(1,:)), max(termination_vector(2, :)),  max(termination_vector(3, :)), threshold, gen );
    end
    if max(termination_vector(1, :)) < threshold && max(termination_vector(2, :)) < threshold && max(termination_vector(3, :))< threshold
        termination_flag = true;
    end
end
end



function [z_ideal, z_nadir, phi] = IGD_Z_termination(NDF_t, NDF_t_1)

ideal_t = min(NDF_t, [], 1);
nadir_t = max(NDF_t, [], 1);
denominator = nadir_t - ideal_t;

% deal with special cases (one ND front then no normalization, if  no feasible solutions, treat same as one ND)
if size(NDF_t, 1) == 1
    denominator = 1;
end

% 
ideal_t_1 = min(NDF_t_1, [], 1);
nadir_t_1 = max(NDF_t_1, [], 1);

z_ideal = abs(ideal_t_1 - ideal_t) ./ denominator;
z_nadir = abs(nadir_t_1 - nadir_t) ./ denominator;

z_ideal = max(z_ideal);
z_nadir = max(z_nadir);

% deal with special cases, two adjacent generation both have one soluiton,
% do not consider its suggestion
if size(NDF_t, 1) == 1 &&  size(NDF_t_1, 1) == 1
    z_ideal = 1;
end

% diversity, normalize on t
NDF_t_norm = (NDF_t - ideal_t) ./ denominator;
NDF_t_1_norm = (NDF_t_1 - ideal_t) ./ denominator;

phi = mean(min(pdist2(NDF_t_norm, NDF_t_1_norm),[],2));

end
