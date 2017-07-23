% Get stats for all symbiont introduction years.


%fls = dir( fullfile( 'D:/GoogleDrive/Coral_Model_Steve/SymbiontGenetics_V9_DualSymbiont/', 'longMortYears*.mat' ) );
fls = dir( fullfile( 'D:/GoogleDrive/Coral_Model_Steve/SymbiontGenetics_V9_DualSymbiont/', 'longMortYears_rcp*.mat' ) );
quants = [0.05 0.25 0.5 0.75 0.95];
fprintf('Symbiont Introduction Based on 5 Years of Mortality\n');
fprintf('Quantiles\n');
fprintf('%7.2f  ',quants);
fprintf('\n');
for ii = 1:numel(fls)
    load(fls(ii).name, 'longMortYears');
    q = quantile(longMortYears, quants);
    fprintf('%7.0f  ', q);
    fprintf('%s\n', fls(ii).name);
end
fprintf('\nMeans\n');
fprintf('Reefs with no   Mean for reefs\n');
fprintf('addition        with addition\n');
for ii = 1:numel(fls)
    load(fls(ii).name, 'longMortYears');
    noAdd = sum(longMortYears ==0);
    longMortYears(longMortYears == 0) = NaN;
    fprintf('%13.0f  %13.0f  ', noAdd, mean(longMortYears, 'omitnan'));
    fprintf('%s\n', fls(ii).name);
end


fls = dir( fullfile( 'D:/GoogleDrive/Coral_Model_Steve/SymbiontGenetics_V9_DualSymbiont/', 'firstBleachYears_rcp*.mat' ) );
quants = [0.05 0.25 0.5 0.75 0.95];
fprintf('\nSymbiont Introduction Based on the First year of bleaching\n');

fprintf('Quantiles\n');
fprintf('%7.2f  ',quants);
fprintf('\n');
for ii = 1:numel(fls)
    load(fls(ii).name, 'firstBleachYears');
    q = quantile(firstBleachYears, quants);
    fprintf('%7.0f  ', q);
    fprintf('%s\n', fls(ii).name);
end
fprintf('\nMeans\n');
fprintf('Reefs with no   Mean for reefs\n');
fprintf('addition        with addition\n');
for ii = 1:numel(fls)
    load(fls(ii).name, 'firstBleachYears');
    noAdd = sum(firstBleachYears ==0);
    firstBleachYears(firstBleachYears == 0) = NaN;
    fprintf('%13.0f  %13.0f  ', noAdd, mean(firstBleachYears, 'omitnan'));
    fprintf('%s\n', fls(ii).name);
end