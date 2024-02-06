

problems = {'DS4(3,2)', 'DS5(3,2)', 'DS4D(3, 2)', 'DS5D(3, 2)', 'DS4(4, 3)', 'DS5(4, 3)', 'DS4D(4, 3)', 'DS5D(4, 3)', 'DS4(5,4)', 'DS5(5,4)', 'DS4D(5, 4)', 'DS5D(5, 4)'};

% problems = {'DS4(5, 4)', 'DS5(5,4)'};
np = length(problems);

seedmax = 21;
median_number = 11;
record_folderT = fullfile(pwd, 'results');

method_causal_igd = ones(np, seedmax) .* inf;
method_ea_igd = ones(np, seedmax) .* inf;
method_dual_igd = ones(np, seedmax) .* inf;

method_casual_ULFE = ones(np, seedmax);
method_casual_LLFE = ones(np, seedmax);

method_ea_ULFE = ones(np, seedmax);
method_ea_LLFE = ones(np, seedmax);

method_dual_ULFE = ones(np, seedmax);
method_dual_LLFE = ones(np, seedmax);

% write into csv file
output_file = sprintf("median_igd.csv");
output_file = fullfile(pwd, "postprocess", output_file);

for ii = 1:np
    prob = eval(problems{ii});
    nl = length(prob.ll_bu);
    for jj = 1:seedmax
        record_file = sprintf("%s_CRchecking_nl_%d_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_causal_igd(ii, jj) = records{3}(end);
        method_casual_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_casual_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records

        record_file = sprintf("%s_EA_nl_%d_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_ea_igd(ii, jj) = records{3}(end);
        method_ea_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_ea_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records

        record_file = sprintf("%s_dualpop_nl_%d_modified_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_dual_igd(ii, jj) = records{3}(end);
        method_dual_ULFE(ii, jj) = sum(records{4}(:, 1));
        method_dual_LLFE(ii, jj) = sum(records{4}(:, 2));
        clear records
    end
end

fp = fopen(output_file, "w");
fprintf(fp, ['problem name, Causal IGD, Dual IGD,  Significance, EA median IGD,  Significance \n' ]);

for ii = 1:np
    prob = eval(problems{ii});
    fprintf(fp, "%s & ", prob.name);

    [median_value0, median_id0] = find_median(method_causal_igd(ii, :), median_number);
    fprintf(fp, "%0.4f &", median_value0);

    [median_value1, median_id1] = find_median(method_dual_igd(ii, :), median_number);
    fprintf(fp, "%0.4f & ", median_value1);
    %fprintf(fp, "%d \t", median_id1);

    [result] = ranksum_sigtest_smallerBetter(method_dual_igd(ii, :), method_causal_igd(ii, :));
    if result == -1
        fprintf(fp, " $\\uparrow_{1}$ &  ");
    elseif result == 1
        fprintf(fp, "$\\downarrow_{1}$ & ");
    else
        fprintf(fp, '$\\approx$');
    end
    
    [median_value2, median_id2] = find_median(method_ea_igd(ii, :), median_number);
    fprintf(fp, "%0.4f &  ", median_value2);

    [result] = ranksum_sigtest_smallerBetter(method_ea_igd(ii, :), method_causal_igd(ii, :));
    if result == -1
        fprintf(fp, "$\\uparrow_{1}$ &  ");
    elseif result == 1
        fprintf(fp, "$\\downarrow_{1}$ &");
    else
        fprintf(fp, '$\\approx$ &');
    end

    
    plot_medianID(median_id0, median_id1, median_id2, prob);
    % plot_truncated_medianID(median_id0, median_id1, median_id2, prob);

    causal_fe = method_casual_ULFE(ii, :) + method_casual_LLFE(ii, :);
    dual_fe = method_dual_ULFE(ii, :) + method_dual_LLFE(ii, :);
    ea_fe = method_ea_ULFE(ii, :) + method_ea_LLFE(ii, :);

    fprintf(fp, '%d &', median(causal_fe));
    fprintf(fp, '%d &', median(dual_fe));
    fprintf(fp, '%d ', median(ea_fe));



    
%     fprintf(fp, "%d ,", median(method_modified_ULFE(ii, :)) );
%     fprintf(fp, "%d ,", median(method_casual_ULFE(ii, :)) );
%     fprintf(fp, "%d ,", median(method_ea_ULFE(ii, :)));
%     fprintf(fp, "%d ,", median(method_modified_LLFE(ii, :)) );
%     fprintf(fp, "%d ,", median(method_casual_LLFE(ii, :)));
%     fprintf(fp, "%d ,", median(method_ea_LLFE(ii, :)));
    fprintf(fp, "\\\\\n");
end
fclose(fp);



output_file = sprintf("truncated_median_igd.csv");
output_file = fullfile(pwd, "postprocess", output_file);

fp = fopen(output_file, "w");
fprintf(fp, 'problem name, Causal modified IGD, Causal median IGD,  Significance,  EA IGD, Significance\n' );
% truncated FE compare
method_truncated_causal_igd = ones(np, seedmax) .* inf;
method_ea_igd = ones(np, seedmax) .* inf;
method_truncated_dualpop_igd = ones(np, seedmax) .* inf;

method_truncated_causal_fe = ones(np, seedmax) .* inf;
method_ea_fe = ones(np, seedmax) .* inf;
method_truncated_dualpop_fe = ones(np, seedmax) .* inf;

for ii = 1:np
    prob = eval(problems{ii});
    fprintf(fp, '%s & ', prob.name);
    nl = length(prob.ll_bu);
    for jj = 1:seedmax
        % base 
        record_file = sprintf("%s_EA_nl_%d_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        method_ea_igd(ii, jj) = records{3}(end);
        FE_perGen = records{4};
        FE_perGen = sum(FE_perGen, 2);
        max_FE = sum(FE_perGen);
        method_ea_fe(ii, jj) = max_FE;
        
        % causal
        clear records
        record_file = sprintf("%s_CRchecking_nl_%d_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        igd_perGen = records{3};
        FE_perGen = records{4};
        FE_perGen = sum(FE_perGen, 2);

        accumulatedFE_perGen = [];
        accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
        for kk = 2:length(FE_perGen)
            accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(kk)];
            if accumulatedFE_perGen(end)>max_FE
                break;
            end
        end
        accumulatedFE_perGen = accumulatedFE_perGen(1:end-1);
        num = length(accumulatedFE_perGen);
        method_truncated_causal_igd(ii, jj) = records{3}(num);
        method_truncated_causal_fe(ii, jj) = accumulatedFE_perGen(end);


        % modified
        clear records
        record_file = sprintf("%s_dualpop_nl_%d_modified_seed_%d.mat", prob.name, nl, jj);
        record_file = fullfile(record_folderT, prob.name, record_file);
        load(record_file);

        igd_perGen = records{3};
        FE_perGen = records{4};
        FE_perGen = sum(FE_perGen, 2);

        accumulatedFE_perGen = [];
        accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
        for kk = 2:length(FE_perGen)
            accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(kk)];
            if accumulatedFE_perGen(end)>max_FE
                break;
            end
        end
        accumulatedFE_perGen = accumulatedFE_perGen(1:end-1);
        num = length(accumulatedFE_perGen);
        method_truncated_dualpop_igd(ii, jj) = records{3}(num);
        method_truncated_dualpop_fe(ii, jj) = accumulatedFE_perGen(end);

        clear records
    end

    [median_value0, median_id0] = find_median(method_truncated_causal_igd(ii, :), median_number);
    fprintf(fp, "%0.4f & ", median_value0);

    [median_value1, median_id1] = find_median(method_truncated_dualpop_igd(ii, :), median_number);
    fprintf(fp, "%0.4f &", median_value1);
    %fprintf(fp, "%d \t", median_id1);

    [result] = ranksum_sigtest_smallerBetter(method_truncated_causal_igd(ii, :), method_truncated_dualpop_igd(ii, :));
    if result == -1
        fprintf(fp, "$\\downarrow_{1}$ & ");
    elseif result == 1
        fprintf(fp, "$\\uparrow_{1}$ & ");
    else
        fprintf(fp, '$\\approx$ &');
    end
    
    [median_value2, median_id2] = find_median(method_ea_igd(ii, :), median_number);
    fprintf(fp, "%0.4f & ", median_value2);

     [result] = ranksum_sigtest_smallerBetter(method_truncated_causal_igd(ii, :), method_ea_igd(ii, :));
    if result == -1
        fprintf(fp, "$\\downarrow_{1}$ & ");
    elseif result == 1
        fprintf(fp, "$\\uparrow_{1}$  & ");
    else
        fprintf(fp, '$\\approx $ & ');
    end

    fprintf(fp, '%d &', median(method_truncated_causal_fe(ii, :)));
    fprintf(fp, '%d &', median(method_truncated_dualpop_fe(ii, :)));
    fprintf(fp, '%d ', median(method_ea_fe(ii, :)));

    fprintf(fp, '\\\\\n');
end
fclose(fp);




function [median_value, median_id] = find_median(vec, median_number)
[vec, ids] = sort(vec);
median_value = vec(median_number);
median_id = ids(median_number);
end

function plot_medianID(id0, id1, id2, prob)

record_folderT = fullfile(pwd, 'results');
nl = length(prob.ll_bu);
record_file = sprintf("%s_CRchecking_nl_%d_seed_%d.mat", prob.name, nl, id0);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);

IGD_perGen = records{3};
f1 = figure('Position',[100, 100, 800, 600]);

FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
end
first_checking = records{5};
accumulatedFE_perGen = [first_checking(1); accumulatedFE_perGen];
IGD_perGen = [first_checking(2), IGD_perGen];
semilogy(accumulatedFE_perGen, IGD_perGen, '-o', 'LineWidth',3, 'MarkerSize',14,'Color', [0 0.4470 0.7410], 'DisplayName', 'BLEA-VAC'); hold on;


clear records;
record_file = sprintf("%s_EA_nl_%d_seed_%d.mat", prob.name, nl, id2);
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
semilogy(accumulatedFE_perGen, IGD_perGen, '-x', 'LineWidth',3, 'MarkerSize',14,'Color', [0.8500 0.3250 0.0980],'DisplayName', 'BLEA'); hold on;


clear records;
record_file = sprintf("%s_dualpop_nl_%d_modified_seed_%d.mat", prob.name, nl, id1);
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
semilogy(accumulatedFE_perGen, IGD_perGen, '-diamond', 'LineWidth',3, 'MarkerSize',14, 'Color', [0.9290 0.6940 0.1250],'DisplayName', 'BLEA-DUA'); hold on;

xlabel('FEs', 'FontSize', 45);
ylabel('IGD', 'FontSize', 45);
% title(prob.name, 'FontSize', 20);
legend('FontSize', 30);

ax = gca;
ax.XAxis.FontSize = 30;
ax.YAxis.FontSize = 30;
box on;
ax.LineWidth = 2;


savename = sprintf("IGD_convergence_plot_%s_nl_%d.fig", prob.name, nl);
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);

savename = sprintf("IGD_convergence_plot_%s_nl_%d.png", prob.name, nl);
savename = fullfile(pwd, 'postprocess', savename);
saveas(f1, savename);

% pause(1);
close(f1);


end



function plot_truncated_medianID(id0, id1, id2, prob)


record_folderT = fullfile(pwd, 'results');
f1 = figure('Position',[100, 100, 800, 600]);
nl = length(prob.ll_bu);
%----- EA ----
record_file = sprintf("%s_EA_nl_%d_seed_%d.mat", prob.name, nl, id2);
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
semilogy(accumulatedFE_perGen, IGD_perGen, '-x', 'LineWidth',3, 'MarkerSize',12,'Color', [0.8500 0.3250 0.0980],'DisplayName', 'BLEA'); hold on;

% set FE base
FE_max = accumulatedFE_perGen(end);
clear record_file

%----- checking ----
record_file = sprintf("%s_CRchecking_nl_%d_seed_%d.mat", prob.name, nl, id0);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);

IGD_perGen = records{3};
FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
    if accumulatedFE_perGen(end)>FE_max
        break;
    end
end
accumulatedFE_perGen = accumulatedFE_perGen(1:end-1);
num = length(accumulatedFE_perGen);
IGD_perGen = IGD_perGen(1:num);
semilogy(accumulatedFE_perGen, IGD_perGen, '-o','LineWidth',3, 'MarkerSize',14,'Color', [0 0.4470 0.7410], 'DisplayName', 'BLEA-VAC'); hold on;


%--- modified-----
clear records;
record_file = sprintf("%s_dualpop_nl_%d_modified_seed_%d.mat", prob.name, nl, id1);
record_file = fullfile(record_folderT, prob.name, record_file);
load(record_file);
IGD_perGen = records{3};

FE_perGen = records{4};
FE_perGen = sum(FE_perGen, 2);

accumulatedFE_perGen = [];
accumulatedFE_perGen = [accumulatedFE_perGen;   FE_perGen(1)];
for ii = 2:length(FE_perGen)
    accumulatedFE_perGen = [accumulatedFE_perGen; accumulatedFE_perGen(end) +  FE_perGen(ii)];
    if accumulatedFE_perGen(end)>FE_max
        break;
    end
end
accumulatedFE_perGen = accumulatedFE_perGen(1:end-1);
num = length(accumulatedFE_perGen);
IGD_perGen = IGD_perGen(1:num);
semilogy(accumulatedFE_perGen, IGD_perGen, '-diamond', 'LineWidth',3, 'MarkerSize',14,'Color', [0.9290 0.6940 0.1250],'DisplayName', 'BLEA-DUA'); hold on;

xlabel('FEs', 'FontSize',45);
ylabel('IGD', 'FontSize', 45);
% title(prob.name, 'FontSize', 30);
legend('FontSize', 30);

ax = gca;
ax.XAxis.FontSize = 30;
ax.YAxis.FontSize = 30;
box on;
ax.LineWidth = 2;

savename = sprintf("IGD_truncated_convergence_plot_%s_nl_%d.fig", prob.name, nl);
savename = fullfile(pwd, 'postprocess', savename);
savefig(savename);

savename = sprintf("IGD_truncated_convergence_plot_%s_nl_%d.png", prob.name, nl);
savename = fullfile(pwd, 'postprocess', savename);
saveas(f1, savename);

close(f1);


end
