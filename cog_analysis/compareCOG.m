%% Change working directory

cd COG_categories/

%% M14 versus AA7

% Load the data
clear;
chromosome = table2cell(readtable('Sinorhizobium_sp_M14.counted.txt', ...
    'ReadVariableNames', false));
chromid = table2cell(readtable('Ensifer_adhaerens_Casida_A.counted.txt', ...
    'ReadVariableNames', false));

% Calculate Fisher statistics
output = cell(size(chromosome, 1) + 1, 9);
output{1,1} = 'COG';
output{1,2} = 'M14_Count';
output{1,3} = 'M14_Total';
output{1,4} = 'M14_Percent';
output{1,5} = 'CasidaA_Count';
output{1,6} = 'CasidaA_Total';
output{1,7} = 'CasidaA_Percent';
output{1,8} = 'Unadjusted_P_Value';
output{1,9} = 'Adjusted_P_Value';
for n = 1:size(chromosome, 1)
    table = table([chromosome{n,2}; chromid{n,2}], [chromosome{n,3}; chromid{n,3}],...
        'VariableNames', {'COG','Total'}, 'RowNames', {'Chromosome','Chromid'});
    [~, pValue] = fishertest(table);
    adj_pValue = pValue * 21;
    output{n+1,1} = chromosome{n,1};
    output{n+1,2} = chromosome{n,2};
    output{n+1,3} = chromosome{n,3};
    output{n+1,4} = round(chromosome{n,4}, 3);
    output{n+1,5} = chromid{n,2};
    output{n+1,6} = chromid{n,3};
    output{n+1,7} = round(chromid{n,4}, 3);
    output{n+1,8} = round(pValue, 3);
    output{n+1,9} = round(adj_pValue, 3);
    clearvars table
end

% Export the data
output = cell2table(output);
writetable(output, 'M14_versus_CasidaA_compared.txt',...
    'WriteVariableNames', false, 'Delimiter', '\t');

%% M14 versus KH-21-134

% Load the data
clear;
chromosome = table2cell(readtable('Sinorhizobium_sp_M14.counted.txt', ...
    'ReadVariableNames', false));
chromid = table2cell(readtable('Ensifer_adhaerens_OV14.counted.txt', ...
    'ReadVariableNames', false));

% Calculate Fisher statistics
output = cell(size(chromosome, 1) + 1, 9);
output{1,1} = 'COG';
output{1,2} = 'M14_Count';
output{1,3} = 'M14_Total';
output{1,4} = 'M14_Percent';
output{1,5} = 'OV14_Count';
output{1,6} = 'OV14_Total';
output{1,7} = 'OV14_Percent';
output{1,8} = 'Unadjusted_P_Value';
output{1,9} = 'Adjusted_P_Value';
for n = 1:size(chromosome, 1)
    table = table([chromosome{n,2}; chromid{n,2}], [chromosome{n,3}; chromid{n,3}],...
        'VariableNames', {'COG','Total'}, 'RowNames', {'Chromosome','Chromid'});
    [~, pValue] = fishertest(table);
    adj_pValue = pValue * 21;
    output{n+1,1} = chromosome{n,1};
    output{n+1,2} = chromosome{n,2};
    output{n+1,3} = chromosome{n,3};
    output{n+1,4} = round(chromosome{n,4}, 3);
    output{n+1,5} = chromid{n,2};
    output{n+1,6} = chromid{n,3};
    output{n+1,7} = round(chromid{n,4}, 3);
    output{n+1,8} = round(pValue, 3);
    output{n+1,9} = round(adj_pValue, 3);
    clearvars table
end

% Export the data
output = cell2table(output);
writetable(output, 'M14_versus_OV14_compared.txt',...
    'WriteVariableNames', false, 'Delimiter', '\t');

%% M14 versus KT2440

% Load the data
clear;
chromosome = table2cell(readtable('Sinorhizobium_sp_M14.counted.txt', ...
    'ReadVariableNames', false));
chromid = table2cell(readtable('Sinorhizobium_sp_A49.counted.txt', ...
    'ReadVariableNames', false));

% Calculate Fisher statistics
output = cell(size(chromosome, 1) + 1, 9);
output{1,1} = 'COG';
output{1,2} = 'M14_Count';
output{1,3} = 'M14_Total';
output{1,4} = 'M14_Percent';
output{1,5} = 'A49_Count';
output{1,6} = 'A49_Total';
output{1,7} = 'A49_Percent';
output{1,8} = 'Unadjusted_P_Value';
output{1,9} = 'Adjusted_P_Value';
for n = 1:size(chromosome, 1)
    table = table([chromosome{n,2}; chromid{n,2}], [chromosome{n,3}; chromid{n,3}],...
        'VariableNames', {'COG','Total'}, 'RowNames', {'Chromosome','Chromid'});
    [~, pValue] = fishertest(table);
    adj_pValue = pValue * 21;
    output{n+1,1} = chromosome{n,1};
    output{n+1,2} = chromosome{n,2};
    output{n+1,3} = chromosome{n,3};
    output{n+1,4} = round(chromosome{n,4}, 3);
    output{n+1,5} = chromid{n,2};
    output{n+1,6} = chromid{n,3};
    output{n+1,7} = round(chromid{n,4}, 3);
    output{n+1,8} = round(pValue, 3);
    output{n+1,9} = round(adj_pValue, 3);
    clearvars table
end

% Export the data
output = cell2table(output);
writetable(output, 'M14_versus_A49_compared.txt',...
    'WriteVariableNames', false, 'Delimiter', '\t');




quit

