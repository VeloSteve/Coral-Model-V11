% This produces tables of bleaching statistics in selected years.  It also
% saves one set (total mortality, bleaching, and frequent bleaching) in a
% .mat file for use in cross-run plotting.
%
% Inputs:
% bleachEvents, bleachState, mortState - organized (reef, year, coral type)
%     and containing logical flags for each state or event.
% lastYearAlive - the last year each reef is alive (year)
% lastBleachEvent - last bleaching event for each year and coral type
%     (year, coral type)
% thisRun - list of reef numbers used in this run, which can serve to index
% the other arrays
% allLatLon - longitude and latitude for all reefs
% outputPath - where to write output
% RCP, E, OA - constants used for labeling
% bleachParams is just passed on to an out .mat file as a record of the
%     setup of this run.
%   
function [permBleached, percentMortality] = ...
        Stats_Tables(bleachEvents, bleachState, mortState, lastYearAlive, ...
        lastBleachEvent, frequentBleaching, thisRun, allLatLon, outputPath, startYear, RCP, E, OA, ...
        bleachParams)
    % Subset all of the input arrays which list all reefs to just those
    % which are active in this run.  We don't care about reef IDs, just the
    % number and their latitude.
    if length(thisRun) < length(allLatLon)
        bleachEvents = bleachEvents(thisRun, :, :);
        bleachState = bleachState(thisRun, :, :);
        mortState = mortState(thisRun, :, :);
        lastYearAlive = lastYearAlive(thisRun);
        lastBleachEvent = lastBleachEvent(thisRun, :);
        frequentBleaching = frequentBleaching(thisRun, :, :);
        % And here we only need latitude. Discard longitude.
        latitude = allLatLon(thisRun, 2);
    else
        latitude = allLatLon(:, 2);
    end
    
    % The number of coral types is one less than the number of the last
    % index of bleachState.
    nc = size(bleachState, 3) - 1;
    ny = size(bleachState, 1);
    nr = size(thisRun);
    
    % Divide the world's reefs into 3 equal parts by latitude.
    % For everyx = 1 an equal split into 3 parts by latitude would be 642.
    % 7, 15 gives the closest match: 627, 654, 644 
    eqLim = 7;
    loLim = 15;
    lower = [0 eqLim loLim 0];
    upper = [eqLim loLim 90 90];
        
    % Warning: the psw2 optimization code assumes that 1950 is in column 2
    % of the output array, so it needs to be first here.
    % Points for graphing and estimating mortality years
    %years = [1950 1965 1980 1990 2000 2010 2016 2020 2030 2033 2040 2050 2060 2070 2075 2085 2095 2100];
    % Extended version for 400 years
    %years = [1860 1875 1880 1885 1890 1895 1900 1925 1950 1980 2000 2010 2016 2020 2030 2033 2040 2050 2060 2070 2075 2085 2095 2100 2125 2150 2175 2200 2205 2210 2215 2220 2225 2250 2260 ];
    % Minimal points for faster optimization runs:
    % years = [1950 2100];
    % Default, for a readable table.
    years = [1950 2000 2016 2050 2075 2100];

    % All tables start out with the same prealocation.
    permBleached = zeros(5,length(years));
    percentMortality = permBleached;
    highFrequency = permBleached;
    unrecovered = permBleached;
    frequentLive = permBleached; % For percent of still-living reefs seeing frequent bleaching.
    allStress = permBleached;
     
    % Use region indexes to select from the various input arrays.
    indEq = find(abs(latitude) <= eqLim);
    indLo = find((abs(latitude) > eqLim) & (abs(latitude) <= loLim));
    indHi = find(abs(latitude) > loLim);
    indexes = {indEq, indLo, indHi, 1:length(thisRun)};

    % Get divisors for each region based on the indexes
    numEq = length(indEq);
    numLo = length(indLo);
    numHi = length(indHi);
    numTotal = numEq + numLo + numHi;
    latCounts = [numEq numLo numHi numTotal];    
    assert(length(latitude) == numEq + numLo + numHi, ...
        'Reefs in subsets should match total reefs in this run. eq/lo/hi = %d/%d/%d, total = %d', ...
        numEq, numLo, numHi, length(latitude));
    
    % This loop builds the interior of all six tables at once.  They are
    % organized by year (columns) and by latitude range (rows).
    for lat = 1:4
        ind = indexes{lat};
        % Get subsets for the current latitude range.
        % XXX Index nc+1 is all reefs, but the old code used just massive!
        % in an initial 97-reef test nc+1 and 1 give the same results here.
        % 2100 values:
        % old code   7/25/2017 1PM
        % 96        56
        % 84.38     34.38
        % 60        55
        % 77.32     48.45
        % Why so different?  The old code counted reefs which actually
        % recovered!
        lastBleachLat = lastBleachEvent(ind, 1);
        lastYearAliveLat = lastYearAlive(ind);
        frequentBleachingLat = frequentBleaching(ind, :, :);
        fbCombo = frequentBleachingLat(:, :, 1) ...
            | frequentBleachingLat(:, :, 2); % either coral counts
        mortLat = mortState(ind, :, :);
        mortCombo = mortLat(:, :, 1) | mortLat(:, :, 2);
        bleachLat = bleachState(ind, :, :);
        bleachCombo = bleachLat(:, :, 1) | bleachLat(:, :, 2);
        bleachOrMort = mortCombo | bleachCombo;
        stressCombo = bleachOrMort | fbCombo;
        
        
        for n = 1:length(years)
            yr = years(n);
            yIndex = yr - startYear + 1;
            permBleached(1, n) = yr;
            permBleached(lat+1, n) = 100*length(lastBleachLat(lastBleachLat <= yr)) / latCounts(lat);

            % Permanent mortality
            percentMortality(1, n) = yr;
            mortalityCount = length(lastYearAliveLat(lastYearAliveLat <= yr));
            percentMortality(lat+1, n) = 100* mortalityCount / latCounts(lat);

            fb = fbCombo(:, yIndex);
            highFrequency(1, n) = yr;
            bc = nnz(fb);
            highFrequency(lat+1, n) = 100 * bc  / latCounts(lat);

            % To consider only living reefs for frequent bleaching, just
            % change the divisor.
            frequentLive(1, n) = yr;
            live = latCounts(lat) - mortalityCount;
            if live > 0
                frequentLive(lat+1, n) = 100 * bc / live;
            else
                frequentLive(lat+1, n) = nan;
            end

            % This is now defined to mean that either coral type is 
            % bleached or dead.
            unrecovered(1, n) = yr;
            bmc = nnz(bleachOrMort(:, yIndex));
            unrecovered(lat+1, n) = 100 * bmc / latCounts(lat);
            
            % The number we plot later - any negative state including
            % frequent bleaching, mortality, or current bleaching.
            allStress(1, n) = yr;
            sc = nnz(stressCombo(:, yIndex));
            allStress(lat+1, n) = 100 * sc / latCounts(lat);
           
        end
    end
       
    % All tables have the same label columns.
    labels = cell(5, 3);
    labels{1,1} = 'Year      ';
    labels{2,1} = 'Equatorial';
    labels{3,1} = 'Low       ';
    labels{4,1} = 'High      ';
    labels{5,1} = 'All Reefs ';
    labels{1, 2} = 'Total Reefs';
    labels{2, 2} = numEq;
    labels{3, 2} = numLo;
    labels{4, 2} = numHi;
    labels{5, 2} = numTotal;
    labels{1, 3} = 'Max Latitude';
    labels{2, 3} = eqLim;
    labels{3, 3} = loLim;
    labels{4, 3} = max(latitude(:));


    logTwo('Permanently bleached reefs as of the date given:\n');
    printTable(labels, permBleached, length(years));

    logTwo('\nPermanent mortality as of the date given:\n');
    printTable(labels, percentMortality, length(years));

    logTwo('\nPercentage of reefs with more than one massive coral bleaching event in the previous 10 years.\n');
    printTable(labels, highFrequency, length(years));

    logTwo('\nPercentage of LIVING reefs with more than one bleaching event in the previous 10 years.\n');
    printTable(labels, frequentLive, length(years));

    logTwo('\nPercentage of reefs with at least one coral type in an unrecovered state.\n');
    printTable(labels, unrecovered, length(years));

    logTwo('\nPercentage of reefs with any form of current bleaching or mortality.\n');
    printTable(labels, allStress, length(years));
    
    
    % Data for plotting bleaching histories.
    % Instead of creating the plot here, save the data for use in an
    % outside program which can compare different runs.
    xForPlot = years;
    yForPlot = allStress(5, :);
    yEq = allStress(2, :);
    yLo = allStress(3, :);
    yHi = allStress(4, :);
    save(strcat(outputPath, 'bleaching/BleachingHistory', RCP, 'E=', num2str(E), 'OA=', num2str(OA), '.mat'), 'xForPlot', 'yForPlot', 'yEq', 'yLo', 'yHi', 'bleachParams');
end


function printTable(labels, num, len)

    % The number of elements in years changes the required format spec.
    format1 = strcat('%s ', repmat(' %6d ', 1, len),  ' %12s %12s\n');
    format2 = strcat('%s ', repmat(' %6.2f ', 1, len),  ' %12d %12.1f\n');
    logTwo(format1, labels{1, 1}, num(1, :), labels{1, 2:3});
    for i = 2:5
        logTwo(format2, labels{i, 1}, num(i, :), labels{i, 2}, labels{i, 3});
    end
end