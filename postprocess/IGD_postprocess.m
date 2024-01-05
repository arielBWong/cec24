problems = {'DS4(5, 4)', 'DS5(5,4)'};
np = length(problems);

seedmax = 21;
median_number = 11;
record_folderT = fullfile(pwd, 'results');

method_causal_igd = ones(np, seedmax) .* inf;
method_ea_igd = ones(np, seedmax) .* inf;
method_modified_igd = ones(np, seedmax) .* inf;

method_casual_ULFE = ones(np, seedmax);
method_casual_LLFE = ones(np, seedmax);

method_ea_ULFE = ones(np, seedmax);
method_ea_LLFE = ones(np, seedmax);

method_modified_ULFE = ones(np, seedmax);
method_modified_LLFE = ones(np, seedmax);



% write into csv file
output_file = sprintf("median_igd.csv");
output_file = fullfile(pwd, "postprocess", output_file);

for ii = 1:np
    prob = eval(problems{ii});
    for jj = 1:seedmax
        record_file = sprintf("%s_CRchecking_seed_%d.mat", prob.name, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_causal_igd(ii, jj) = records{3}(end);
        method_casual_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_casual_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records

        record_file = sprintf("%s_EA_seed_%d.mat", prob.name, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_ea_igd(ii, jj) = records{3}(end);
        method_ea_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_ea_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records


        record_file = sprintf("%s_CRchecking_modified_seed_%d.mat", prob.name, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_modified_igd(ii, jj) = records{3}(end);
        method_modified_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_modified_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records
    end

end

fp = fopen(output_file, "w");
fprintf(fp, ['problem name, Causal modified IGD, Causal median IGD, EA median IGD,  Causal modified ULFE, causal ULFE,  EA ULFE, ...' ...
    'Causal modified LLFE, Causal LLFE, EA LLFE \n' ]);

for ii = 1:np
    prob = eval(problems{ii});
    fprintf(fp, "%s,", prob.name);

    [median_value, median_id0] = find_median(method_modified_igd(ii, :), median_number);
    fprintf(fp, "%0.4f ,", median_value);

    [median_value, median_id1] = find_median(method_causal_igd(ii, :), median_number);
    fprintf(fp, "%0.4f ,", median_value);
    %fprintf(fp, "%d \t", median_id1);
    
    [median_value, median_id2] = find_median(method_ea_igd(ii, :), median_number);
    fprintf(fp, "%0.4f ,", median_value);
    % fprintf(fp, "%d \t", median_id2);

    plot_medianID(median_id0, median_id1, median_id2, prob);
    
    fprintf(fp, "%d ,", median(method_modified_ULFE(ii, :)) );
    fprintf(fp, "%d ,", median(method_casual_ULFE(ii, :)) );
    fprintf(fp, "%d ,", median(method_ea_ULFE(ii, :)));
    fprintf(fp, "%d ,", median(method_modified_LLFE(ii, :)) );
    fprintf(fp, "%d ,", median(method_casual_LLFE(ii, :)));
    fprintf(fp, "%d ,", median(method_ea_LLFE(ii, :)));
    fprintf(fp, "\n");
end
fclose(fp);

function [median_value, median_id] = find_median(vec, median_number)
[vec, ids] = sort(vec);
median_value = vec(median_number);
median_id = ids(median_number);
end

function plot_medianID(id0, id1, id2, prob)

record_folderT = fullfile(pwd, 'results');

record_file = sprintf("%s_CRchecking_seed_%d.mat", prob.name, id1);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);

IGD_perGen = records{3};
f1 = figure('Position',[100, 100, 600, 600]);

FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
end
semilogy(accumulatedFE_perGen, IGD_perGen, '-o', 'LineWidth',2,'Color', [0 0.4470 0.7410], 'DisplayName', 'Causal Check'); hold on;


clear records;
record_file = sprintf("%s_EA_seed_%d.mat", prob.name, id2);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);
IGD_perGen = records{3};

FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
end
semilogy(accumulatedFE_perGen, IGD_perGen, '-x', 'LineWidth',2,'Color', [0.8500 0.3250 0.0980],'DisplayName', 'Baseline EA'); hold on;


clear records;
record_file = sprintf("%s_CRchecking_modified_seed_%d.mat", prob.name, id0);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);
IGD_perGen = records{3};

FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
end
semilogy(accumulatedFE_perGen, IGD_perGen, '-diamond', 'LineWidth',2,'Color', [0.9290 0.6940 0.1250],'DisplayName', 'Causal modified'); hold on;

xlabel('Generations', 'FontSize',20);
ylabel('IGD', 'FontSize', 20);
title(prob.name, 'FontSize', 20);
legend('FontSize', 20);

ax = gca;
ax.XAxis.FontSize = 15;
ax.YAxis.FontSize = 15;

savename = sprintf("IGD_convergence_plot_%s.fig", prob.name);
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);

close(f1);


end




