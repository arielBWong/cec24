function [termination_flag, termination_vector] = Termination_criterion_HV(archive, sliding_window, threshold, gen, termination_vector, verbose)

termination_flag = false;

% decide which element to fill
vector_id =  mod(gen, size(termination_vector, 2)); 
if vector_id == 0
    vector_id = 10;
end

% decide reference point
if gen >= sliding_window
    considered_generations = gen: -1: (gen-sliding_window+1);
    considered_archives =  archive(considered_generations);
else
    considered_archives = archive;
end
FU = [considered_archives.FU];
FU = cat(1, FU{:});
ref = max(FU, [], 1);
current_FU = archive(end).FUs;
hv = Hypervolume(current_FU, ref);
termination_vector(1, vector_id) = hv;

if gen <=  sliding_window
    termination_flag = false;
else
    HV_max = max(termination_vector(1,:));
    HV_min = min(termination_vector(1,:));
    H =(HV_max - HV_min) / (HV_max + HV_min);
    
    if verbose
        fprintf('Termination condition : %0.4f, compared to %0.3f, generation %d \n',H, threshold, gen);
    end

    if H<= threshold
        termination_flag = true;
    end 
end

end
