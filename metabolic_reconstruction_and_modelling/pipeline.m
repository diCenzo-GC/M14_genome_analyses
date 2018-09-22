%% Init and load data

% Start cobra toolbox
initCobraToolbox();
changeCobraSolver('gurobi', 'all');
addpath(genpath('Tn-Core-v1.2'));

% Load the ortholog list
orthologs = table2cell(readtable('orthologs.txt', 'Delimiter', '\t', 'ReadVariableNames', false));

% Load the fixed name list for M14
M14names = table2cell(readtable('M14_names.txt', 'Delimiter', '\t', 'ReadVariableNames', false));

% Load the metabolic models
iGD726 = readCbModel('iGD726.xml');
iGD1575 = readCbModel('iGD1575.xml');
model = readCbModel('M14_reAnnotated_model.xml');

%% Update gene names in the M14 model

% Rename the CDS genes in the gene list
genesNew = {};
for n = 1:length(model.genes)
    pos = strmatch(model.genes{n}, M14names(:,1), 'exact');
    if ~isempty(pos)
        genesNew{n,1} = M14names{pos,2};
    else
        genesNew{n,1} = model.genes{n};
    end
end
genes = horzcat(model.genes, genesNew);
genes = sortrows(genes, -1);
model.genes = genesNew;

% Fix the grRules
for n = 1:length(genes)
    for m = 1:length(model.rxns)
        if strfind(model.grRules{m}, genes{n,1})
            model.grRules{m} = strrep(model.grRules{m}, genes{n,1}, genes{n,2});
        end
    end
end

% Remove duplicate genes
model = tncore_duplicates(model);

%% Remove non-orthologous genes from meliloti models

% Rename Rm1021 names with M14 names in iGD726
genesNew = {};
for n = 1:length(iGD726.genes)
    pos = strmatch(iGD726.genes{n}, orthologs(:,2), 'exact');
    if ~isempty(pos)
        genesNew{n,1} = orthologs{pos,1};
    else
        genesNew{n,1} = iGD726.genes{n};
    end
end
genes = horzcat(iGD726.genes, genesNew);
genes = sortrows(genes, -1);
iGD726.genes = genesNew;

% Fix the grRules in iGD726
for n = 1:length(genes)
    for m = 1:length(iGD726.rxns)
        if strfind(iGD726.grRules{m}, genes{n,1})
            iGD726.grRules{m} = strrep(iGD726.grRules{m}, genes{n,1}, genes{n,2});
        end
    end
end

% Remove rxns without orthologs in iGD726
posC = strmatch('smc', iGD726.genes);
posB = strmatch('smb', iGD726.genes);
posA = strmatch('sma', iGD726.genes);
genesToDeleteA = iGD726.genes(posA);
genesToDeleteB = iGD726.genes(posB);
genesToDeleteC = iGD726.genes(posC);
genesToDelete = unique(union(union(genesToDeleteA, genesToDeleteB), genesToDeleteC));
[iGD726_deleted, ~, constrRxnNames] = deleteModelGenes(iGD726, genesToDelete);
iGD726_deleted = removeRxns(iGD726_deleted, constrRxnNames);
iGD726_deleted = tncore_delete(iGD726_deleted);
save('temp.mat');

% Rename Rm1021 names with M14 names in iGD1575
genesNew = {};
for n = 1:length(iGD1575.genes)
    pos = strmatch(iGD1575.genes{n}, orthologs(:,2), 'exact');
    if ~isempty(pos)
        genesNew{n,1} = orthologs{pos,1};
    else
        genesNew{n,1} = iGD1575.genes{n};
    end
end
genes = horzcat(iGD1575.genes, genesNew);
genes = sortrows(genes, -1);
iGD1575.genes = genesNew;

% Fix the grRules in iGD1575
for n = 1:length(genes)
    for m = 1:length(iGD1575.rxns)
        if strfind(iGD1575.grRules{m}, genes{n,1})
            iGD1575.grRules{m} = strrep(iGD1575.grRules{m}, genes{n,1}, genes{n,2});
        end
    end
end

% Remove rxns without orthologs in iGD1575
posC = strmatch('smc', iGD1575.genes);
posB = strmatch('smb', iGD1575.genes);
posA = strmatch('sma', iGD1575.genes);
genesToDeleteA = iGD1575.genes(posA);
genesToDeleteB = iGD1575.genes(posB);
genesToDeleteC = iGD1575.genes(posC);
genesToDelete = unique(union(union(genesToDeleteA, genesToDeleteB), genesToDeleteC));
[iGD1575_deleted, ~, constrRxnNames] = deleteModelGenes(iGD1575, genesToDelete);
iGD1575_deleted = removeRxns(iGD1575_deleted, constrRxnNames);
iGD1575_deleted = tncore_delete(iGD1575_deleted);
save('temp2.mat');

%% Expand the M14 model

% Expand with both models
M14_expanded = tncore_expand(model, iGD726_deleted, 0);
M14_expanded_noDeadends = tncore_expand(M14_expanded, iGD1575_deleted, 1);
M14_expanded_withDeadends = tncore_expand(M14_expanded, iGD1575_deleted, 0);
M14_expanded_noDeadends = rmfield(M14_expanded_noDeadends, 'rxnNotes');
M14_expanded_withDeadends = rmfield(M14_expanded_withDeadends, 'rxnNotes');

% Remove demand reactions
demandRxns = M14_expanded_noDeadends.rxns(strmatch('Demand_', M14_expanded_noDeadends.rxns));
M14_expanded_noDeadends = removeRxns(M14_expanded_noDeadends, demandRxns);
demandRxns = M14_expanded_withDeadends.rxns(strmatch('Demand_', M14_expanded_withDeadends.rxns));
M14_expanded_withDeadends = removeRxns(M14_expanded_withDeadends, demandRxns);

% Remove reactions differing only by a proton
rxnsToRemove = {};
for n = 1:length(M14_expanded_noDeadends.rxns)
    sameName = M14_expanded_noDeadends.rxns(strmatch(M14_expanded_noDeadends.rxns{n}, M14_expanded_noDeadends.rxns));
    if length(sameName) > 1
        for m = 2:length(sameName)
            [metList_1] = parseRxnFormula(cell2mat(printRxnFormula(M14_expanded_noDeadends, sameName{1,1}))); 
            [metList_2] = parseRxnFormula(cell2mat(printRxnFormula(M14_expanded_noDeadends, sameName{m,1}))); 
            metList_union = union(metList_1, metList_2);
            metList_inter = intersect(metList_1, metList_2);
            metlist_diff = setdiff(metList_union, metList_inter);
            if isempty(metlist_diff)
                rxnsToRemove = vertcat(rxnsToRemove, sameName{m});
            elseif length(metlist_diff) == 1
                if strmatch('cpd00067', metlist_diff)
                    rxnsToRemove = vertcat(rxnsToRemove, sameName{m});
                end
            end
        end
    end
end
M14_expanded_noDeadends = removeRxns(M14_expanded_noDeadends, rxnsToRemove);
M14_expanded_noDeadends = tncore_remove(M14_expanded_noDeadends);

rxnsToRemove = {};
for n = 1:length(M14_expanded_withDeadends.rxns)
    sameName = M14_expanded_withDeadends.rxns(strmatch(M14_expanded_withDeadends.rxns{n}, M14_expanded_withDeadends.rxns));
    if length(sameName) > 1
        for m = 2:length(sameName)
            [metList_1] = parseRxnFormula(cell2mat(printRxnFormula(M14_expanded_withDeadends, sameName{1,1}))); 
            [metList_2] = parseRxnFormula(cell2mat(printRxnFormula(M14_expanded_withDeadends, sameName{m,1}))); 
            metList_union = union(metList_1, metList_2);
            metList_inter = intersect(metList_1, metList_2);
            metlist_diff = setdiff(metList_union, metList_inter);
            if isempty(metlist_diff)
                rxnsToRemove = vertcat(rxnsToRemove, sameName{m});
            elseif length(metlist_diff) == 1
                if strmatch('cpd00067', metlist_diff)
                    rxnsToRemove = vertcat(rxnsToRemove, sameName{m});
                end
            end
        end
    end
end
M14_expanded_withDeadends = removeRxns(M14_expanded_withDeadends, rxnsToRemove);
M14_expanded_withDeadends = tncore_remove(M14_expanded_withDeadends);

%% Save

save('allWorkspace.mat');

%% Biolog

% Get exchange reactions
EXlist = M14_expanded_noDeadends.rxns(strmatch('EX_', M14_expanded_noDeadends.rxns));
EXlist2 = M14_expanded_noDeadends.rxns(strmatch('Demand_', M14_expanded_noDeadends.rxns));
EXlist3 = iGD1575.rxns(strmatch('EX_', iGD1575.rxns));

% Prepare base medium
carbonFreeMinimal = {'EX_cpd00001_e0','EX_cpd00007_e0','EX_cpd00099_e0'...
    'EX_cpd00009_e0','EX_cpd00149_e0','EX_cpd00013_e0','EX_cpd00030_e0',...
    'EX_cpd00034_e0','EX_cpd00048_e0','EX_cpd00058_e0','EX_cpd00063_e0',...
    'EX_cpd00067_e0','EX_cpd00099_e0','EX_cpd00205_e0','EX_cpd00254_e0',...
    'EX_cpd00971_e0','EX_cpd10515_e0','EX_cpd10516_e0'};

% Set model boundaries
M14_expanded_noDeadends = changeRxnBounds(M14_expanded_noDeadends, EXlist,0,'l');
M14_expanded_noDeadends = changeRxnBounds(M14_expanded_noDeadends, EXlist2,0,'u');
M14_expanded_noDeadends = changeRxnBounds(M14_expanded_noDeadends, carbonFreeMinimal, -10, 'l');
iGD1575 = changeRxnBounds(iGD1575, EXlist3,0,'l');
iGD1575 = changeRxnBounds(iGD1575, carbonFreeMinimal, -10, 'l');

% Set the objective function of the models
M14_expanded_noDeadends = changeObjective(M14_expanded_noDeadends,'bio1');
iGD1575 = changeObjective(iGD1575,'biomass_bulk_c0');

% Define list of the models
testModels = { M14_expanded_noDeadends; iGD1575};

% Model names
models = { 'M14_expanded_noDeadends'; 'iGD1575'};

% Define list of the biolog test compounds
biolog = {'EX_cpd00020_e0','EX_cpd00023_e0','EX_cpd00027_e0',...
    'EX_cpd00029_e0','EX_cpd00033_e0','EX_cpd00035_e0','EX_cpd00036_e0',...
    'EX_cpd00041_e0','EX_cpd00047_e0','EX_cpd00051_e0','EX_cpd00053_e0',...
    'EX_cpd00054_e0','EX_cpd00064_e0','EX_cpd00069_e0','EX_cpd00076_e0',...
    'EX_cpd00082_e0','EX_cpd00100_e0','EX_cpd00105_e0','EX_cpd00106_e0',...
    'EX_cpd00107_e0','EX_cpd00108_e0','EX_cpd00119_e0','EX_cpd00121_e0',...
    'EX_cpd00129_e0','EX_cpd00130_e0','EX_cpd00138_e0','EX_cpd00142_e0',...
    'EX_cpd00154_e0','EX_cpd00155_e0','EX_cpd00156_e0','EX_cpd00158_e0',...
    'EX_cpd00159_e0','EX_cpd00161_e0','EX_cpd00179_e0','EX_cpd00182_e0',...
    'EX_cpd00208_e0','EX_cpd00224_e0','EX_cpd00246_e0','EX_cpd00281_e0',...
    'EX_cpd00308_e0','EX_cpd00314_e0','EX_cpd00322_e0','EX_cpd00366_e0',...
    'EX_cpd00382_e0','EX_cpd00386_e0','EX_cpd00392_e0','EX_cpd00396_e0',...
    'EX_cpd00417_e0','EX_cpd00567_e0','EX_cpd00588_e0','EX_cpd00794_e0',...
    'EX_cpd00797_e0','EX_cpd00851_e0','EX_cpd01133_e0','EX_cpd01307_e0',...
    'EX_cpd02175_e0','EX_cpd03198_e0','EX_cpd05158_e0','EX_cpd11585_e0',...
    'EX_cpd11588_e0','EX_cpd11589_e0','EX_cpd11592_e0','EX_cpd00024_e0',...
    'EX_cpd00040_e0','EX_cpd00060_e0','EX_cpd00066_e0','EX_cpd00072_e0',...
    'EX_cpd00079_e0','EX_cpd00089_e0','EX_cpd00094_e0','EX_cpd00136_e0',...
    'EX_cpd00137_e0','EX_cpd00139_e0','EX_cpd00141_e0','EX_cpd00157_e0',...
    'EX_cpd00164_e0','EX_cpd00180_e0','EX_cpd00211_e0','EX_cpd00212_e0',...
    'EX_cpd00222_e0','EX_cpd00232_e0','EX_cpd00248_e0','EX_cpd00266_e0',...
    'EX_cpd00280_e0','EX_cpd00306_e0','EX_cpd00320_e0','EX_cpd00339_e0',...
    'EX_cpd00361_e0','EX_cpd00374_e0','EX_cpd00380_e0','EX_cpd00432_e0',...
    'EX_cpd00438_e0','EX_cpd00453_e0','EX_cpd00477_e0','EX_cpd00489_e0',...
    'EX_cpd00550_e0','EX_cpd00599_e0','EX_cpd00607_e0','EX_cpd00609_e0',...
    'EX_cpd00611_e0','EX_cpd00666_e0','EX_cpd00728_e0','EX_cpd00750_e0',...
    'EX_cpd01055_e0','EX_cpd01067_e0','EX_cpd01107_e0','EX_cpd01113_e0',...
    'EX_cpd01246_e0','EX_cpd01363_e0','EX_cpd01502_e0','EX_cpd01799_e0',...
    'EX_cpd02143_e0','EX_cpd02351_e0','EX_cpd02599_e0','EX_cpd03161_e0',...
    'EX_cpd03561_e0','EX_cpd03696_e0','EX_cpd03734_e0','EX_cpd03737_e0',...
    'EX_cpd05161_e0','EX_cpd05192_e0','EX_cpd05240_e0','EX_cpd09457_e0',...
    'EX_cpd10719_e0','EX_cpd11602_e0','EX_cpd11685_e0','EX_cpd11717_e0',...
    'EX_cpd11879_e0','EX_cpd13391_e0','EX_cpd13392_e0','EX_cpd16821_e0',...
    'EX_cpd00039_e0','EX_cpd00122_e0','EX_cpd00132_e0','EX_cpd00162_e0',...
    'EX_cpd00185_e0','EX_cpd00249_e0','EX_cpd00276_e0','EX_cpd00492_e0',...
    'EX_cpd00589_e0','EX_cpd00737_e0','EX_cpd00751_e0','EX_cpd00832_e0',...
    'EX_cpd01030_e0','EX_cpd01171_e0','EX_cpd01200_e0','EX_cpd01262_e0',...
    'EX_cpd02274_e0','EX_cpd03884_e0','EX_cpd04349_e0','EX_cpd11594_e0',...
    'EX_cpd11601_e0','EX_cpd11603_e0','EX_cpd11748_e0','EX_cpd00080_e0',...
    'EX_cpd00117_e0','EX_cpd00118_e0','EX_cpd00184_e0','EX_cpd00227_e0',...
    'EX_cpd01242_e0','EX_cpd01293_e0','EX_cpd00197_e0','EX_cpd01524_e0'};
%biolog = intersect(biolog, EXlist);

% Prepare output matrix
results = cell(length(biolog)+1,length(testModels)+3);

% Call metabolites the way they are called in the actual biolog experiments
iGD1575.metNames = lower(iGD1575.metNames);
iGD1575.metNames = strrep(iGD1575.metNames,'_e0','');
iGD1575.metNames = strrep(iGD1575.metNames,'_','-');
iGD1575.metNames = regexprep(iGD1575.metNames,'(.*)ate$','$1ic acid');
iGD1575.metNames = regexprep(iGD1575.metNames,'hydroxy ([^\d]*)','hydroxy$1');
iGD1575.metNames = strrep(iGD1575.metNames,'tartric','tartaric');
iGD1575.metNames = strrep(iGD1575.metNames,'gly-pro-l','gly-pro');
iGD1575.metNames = strrep(iGD1575.metNames,'gly-asp-l','gly-asp');
iGD1575.metNames = strrep(iGD1575.metNames,'gly-glu-l','gly-glu');
iGD1575.metNames = strrep(iGD1575.metNames,'phosphic acid','phosphate');

% Find out the real names of the compounds
metIDs = strrep(biolog,'EX_','');
metIDs = strrep(metIDs,'_e0','[e]');
metNames = cell(length(metIDs),3);
metNames(:,2) = strrep(metIDs,'_e0','');
metNames(:,3) = strrep(metIDs,'_e0','');
for n = 1:length(metIDs)
    metNames{n,1} = iGD1575.metNames(strmatch(metIDs(n),iGD1575.mets,'exact'));
end

% Add headers to output
results{1,1} = 'exchangeRxn';
results{1,2} = 'cpdNumber';
results{1,3} = 'metaboliteName';

% Do phenotype microarray and save output
for n = 1:length(testModels)
    results{1,n+3} = models{n};
    for m = 1:length(biolog)
        tmpTest = changeRxnBounds(testModels{n}, biolog{m}, -1, 'l');
        growth = optimizeCbModel(tmpTest, 'max');
        results{m+1,n+3} = round(growth.f, 4);
        results{m+1,2} = metNames{m,2};
        results{m+1,3} = metNames{m,1};
        results{m+1,1} = biolog(m);
    end
end

% Transfer results to table
outputBIOLOG = cell2table(results);

save('Biolog_output.mat');
writetable(outputBIOLOG, 'biologResults.xlsx', 'Sheet', 1, ...
    'WriteVariableNames', false);
clear;









