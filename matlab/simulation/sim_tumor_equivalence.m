clear;
close all;

opts = spreadsheetImportOptions("NumVariables", 125);

% Specify sheet and range
opts.Sheet = "Table S5";
opts.DataRange = "A4:DU197";

% Specify column names and types
opts.VariableNames = ["SampleID", "Tissue", "Patient", "LiveCells", "EndothelialCells", "EpithelialCells", "Fibroblasts", "ImmuneCells", "Other", "FibroblastsFAPSMA", "FibroblastsFAPSMA1", "FibroblastsFAPSMA2", "FibroblastsFAPSMA3", "TCells", "NaturalKillerCells", "Granulocytes", "BCells", "PlasmaCells", "plasmacytoidDendriticCells", "MyeloidCells", "Basophils", "TCellsCD4", "TCellsCD8", "TCellsPD1", "TAMsPDL1", "T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "Ep01", "Ep02", "Ep03", "Ep04", "Ep05", "Ep06", "Ep07", "Ep08", "Ep09", "Ep10", "Ep11", "Ep12", "Ep13", "Ep14", "Ep15", "Ep16", "Ep17", "Ep18", "Ep19", "Ep20", "Ep21", "Ep22", "Ep23", "Ep24", "Ep25", "Ep26", "Ep27", "Ep28", "Ep29", "Ep30", "Ep31", "Ep32", "Ep33", "Ep34", "Ep35", "Ep36", "Ep37", "Ep38", "Ep39", "Ep40", "Ep41", "Ep42", "Ep43", "Ep44", "Ep45", "EpithelialCellsKi67", "EpithelialCellsCA9", "EpithelialCellsERa", "EpithelialCellsERadim", "EpithelialCellsERabright", "EpithelialCellsPRB", "EpithelialCellsPRBdim", "EpithelialCellsPRBbright", "EpithelialCellsHER2", "EpithelialCellsHER2dim", "EpithelialCellsHER2bright", "PhenotypicAbnormalityScore", "IndividualityScore", "LargestCluster", "EcosystemClass", "ImmuneClassTIG"];
opts.SelectedVariableNames = ["SampleID", "Tissue", "Patient", "LiveCells", "EndothelialCells", "EpithelialCells", "Fibroblasts", "ImmuneCells", "Other", "FibroblastsFAPSMA", "FibroblastsFAPSMA1", "FibroblastsFAPSMA2", "FibroblastsFAPSMA3", "TCells", "NaturalKillerCells", "Granulocytes", "BCells", "PlasmaCells", "plasmacytoidDendriticCells", "MyeloidCells", "Basophils", "TCellsCD4", "TCellsCD8", "TCellsPD1", "TAMsPDL1", "T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09", "T10", "T11", "T12", "T13", "T14", "T15", "T16", "T17", "T18", "T19", "T20", "M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "Ep01", "Ep02", "Ep03", "Ep04", "Ep05", "Ep06", "Ep07", "Ep08", "Ep09", "Ep10", "Ep11", "Ep12", "Ep13", "Ep14", "Ep15", "Ep16", "Ep17", "Ep18", "Ep19", "Ep20", "Ep21", "Ep22", "Ep23", "Ep24", "Ep25", "Ep26", "Ep27", "Ep28", "Ep29", "Ep30", "Ep31", "Ep32", "Ep33", "Ep34", "Ep35", "Ep36", "Ep37", "Ep38", "Ep39", "Ep40", "Ep41", "Ep42", "Ep43", "Ep44", "Ep45", "EpithelialCellsKi67", "EpithelialCellsCA9", "EpithelialCellsERa", "EpithelialCellsERadim", "EpithelialCellsERabright", "EpithelialCellsPRB", "EpithelialCellsPRBdim", "EpithelialCellsPRBbright", "EpithelialCellsHER2", "EpithelialCellsHER2dim", "EpithelialCellsHER2bright", "PhenotypicAbnormalityScore", "IndividualityScore", "LargestCluster", "EcosystemClass", "ImmuneClassTIG"];
opts.VariableTypes = ["string", "categorical", "string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical"];
opts = setvaropts(opts, [1, 3], "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 2, 3, 124, 125], "EmptyFieldRule", "auto");

% Import the data
data = readtable("1-s2.0-S0092867419302673-mmc5.xlsx", opts, "UseExcel", false);

clear opts

tumor = data.BCells(data.Tissue == "Tumor") / 100;
juxta = data.BCells(data.Tissue == "Juxta-tumoral") / 100;

[R, P] = ttest2(tumor, juxta)
ksdensity(tumor)
hold on
ksdensity(juxta)

tumor_train = tumor;%(1 : 50);
tumor_test = tumor;%(50 + 1 : end);

juxta_train = juxta;%(1 : 15);
juxta_test = juxta;%(15 + 1 : end);

mu1 = mean(tumor_train);
mu0 = mean(juxta_train);

sigma1 = std(tumor_train);
sigma0 = std(juxta_train);

N0 = 1000;
N1 = 1000;

lower = -0.05;
upper = 0.20;

alpha = 0.05;

%res_mat0 = table();
%res_mat1 = table();

for M1 = 12:2:20
    for M0 = 12:2:20

        a1 = ((1 - mu1) / sigma1 / sigma1 - 1 / mu1) * mu1 * mu1;
        b1 = a1 * (1 / mu1 - 1);

        a0 = ((1 - mu0) / sigma0 / sigma0 - 1 / mu0) * mu0 * mu0;
        b0 = a0 * (1 / mu0 - 1);

        a = [a0 a1];
        b = [b0 b1];
        N = [N0 N1];
        M = [M0 M1];

        Ep = a ./ (a + b);
        Vp = a .* b .* (a + b + N) ./ ((N .* (a + b) .^ 2) .* (a + b + 1));

        nu = sum(Vp ./ M) ^ 2 / sum((Vp ./ M) .^ 2 ./ (M - 1));
        Et_lower = (Ep(2) - (Ep(1) + lower)) / sqrt(sum(Vp ./ M));
        Et_upper = - (Ep(2) - (Ep(1) + upper)) / sqrt(sum(Vp ./ M));
        t_star = tinv(1 - alpha, nu);

        % res_mat0(M0, M1) = min([1, normcdf(t_star - Et_lower) + normcdf(t_star - Et_upper)]);
        res_mat0(M0, M1) = min([1, tcdf(t_star - Et_lower, nu) + tcdf(t_star - Et_upper, nu)]);
        
        rng(2020);
        N_ITER = 10000;
        P_lower = zeros(1, N_ITER);
        P_upper = zeros(1, N_ITER);
        
        
        for ii = 1:N_ITER
            sample_mask1 = randsample(length(tumor_test), M1);
            sample_mask0 = randsample(length(juxta_test), M0);
            p1 = binornd(N1, tumor_test(sample_mask1)) / N1;
            p0 = binornd(N0, juxta_test(sample_mask0)) / N0;

            m1 = mean(p1);
            m0 = mean(p0);
            s1 = std(p1) .^ 2;
            s0 = std(p0) .^ 2;

            nu = (s0 / M0 + s1 / M1) .^ 2 ./ ((s0 / M0) .^ 2 / (M0 - 1) + (s1 / M1) .^ 2 / (M1 - 1));

            T_lower = (m1 - (m0 + lower)) ./ sqrt(s1 / M1 + s0 / M0);
            T_upper = -(m1 - (m0 + upper)) ./ sqrt(s1 / M1 + s0 / M0);
            
            P_lower(ii) = tcdf(T_lower, nu, 'upper');
            P_upper(ii) = tcdf(T_upper, nu, 'upper');
        end
        res_mat1(M0, M1) = 1 - sum((P_lower < alpha) & (P_upper < alpha)) / N_ITER;
        
    end
end

