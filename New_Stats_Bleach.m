function [percentBleached, percentMortality] = New_Stats_Bleach(mEvents, thisRun, allLatLon, ...
    outputPath, bleachParams, RCP, E, OA) 
    % Careful - mEvents here is bEvents in the calling code.
    % mEvents is the list of all bleaching and mortality-related events for
    % this run.
    % thisRun is the list of reefs simulated in this run. 
    % allLatLon is the full list of locations for all 1925 reefs.
    % For everyx = 1:
    % an equal split into 3 parts by latitude would be 642
    % 7, 15 gives the closest match: 627, 654, 644 
    eqLim = 7;
    loLim = 15;
    lower = [0 eqLim loLim 0];
    upper = [eqLim loLim 90 90];
    % Look for when massive corals die for the last time.  We really only
    % want to count when BOTH corals die...  TODO
%TODO
%HERE - changing MASS to BRAN increase the curve steepness, which is good.
% XXX did this contribute to making it hard to reach the 2% target?
%next, try looking at both types (but must avoid double counting)
    lastReefBleach = mEvents([mEvents.last] == 1 & strcmp({mEvents.coral}, 'MASS') & strcmp({mEvents.event}, 'BLEACH'));
    lastReefMort = mEvents([mEvents.last] == 1 & strcmp({mEvents.coral}, 'REEF') & strcmp({mEvents.event}, 'MORTALITY'));
    anyMassBleach = mEvents(strcmp({mEvents.coral}, 'MASS') & strcmp({mEvents.event}, 'BLEACH'));
    anyBranBleach = mEvents(strcmp({mEvents.coral}, 'BRAN') & strcmp({mEvents.event}, 'BLEACH'));
    
    % For anyBleach, we only want one event per year per reef, even if both
    % branching and massive bleached.
    anyBleach = onePerYearPerReef(anyMassBleach, anyBranBleach);
    
    % Now look for a state rather than event - reef has at least one coral
    % type with a bleaching more recent than the last REVIVE or RECOVER, so
    % get all of those.
    % Just the recoveries and revivals:
    bleachStateMass = mEvents(strcmp({mEvents.coral}, 'MASS') ...
        & (strcmp({mEvents.event}, 'RECOVER') | strcmp({mEvents.event}, 'REVIVE') ...
        | strcmp({mEvents.event}, 'MORTALITY')) );
    bleachStateBran = mEvents(strcmp({mEvents.coral}, 'BRAN') ...
        & (strcmp({mEvents.event}, 'RECOVER') | strcmp({mEvents.event}, 'REVIVE') ...
        | strcmp({mEvents.event}, 'MORTALITY')) );
    % Add in the other relevant events from above
    bleachStateMass = [bleachStateMass anyMassBleach];
    bleachStateBran = [bleachStateBran anyBranBleach];
    

    
    % Warning: the psw2 optimization code assumes that 1950 is in column 2
    % of the output array, so it needs to be first here.
    % Warning: expanding this array works automatically below EXCEPT that
    % the fprintf lines will need extra copies of the repeated format
    % elements.
    % Enough points for comparison to Logan et al. 2013:
    % years = [1950 1980 2000 2010 2016 2033 2050 2060 2075 2100];
    % and even more for estimating mortality years
    %years = [1950 1965 1980 1990 2000 2010 2016 2020 2030 2033 2040 2050 2060 2070 2075 2085 2095 2100];

    % Extended version for 400 years
    %years = [1860 1875 1880 1885 1890 1895 1900 1925 1950 1980 2000 2010 2016 2020 2030 2033 2040 2050 2060 2070 2075 2085 2095 2100 2125 2150 2175 2200 2205 2210 2215 2220 2225 2250 2260 ];
    % Minimal points for faster optimization runs:
    % years = [1950 2100];
    % For control400 only, end at 2260.
    % Default, for a readable table.
    %years = [1950 2000 2016 2050 2100];
    years = [1950 2000 2016 2050 2075 2100];

    % All tables start out with the same size.
    percentBleached = zeros(5,length(years));
    percentMortality = percentBleached;
    highFrequency = percentBleached;
    unrecovered = percentBleached;
    frequentLive = percentBleached; % For percent of still-living reefs seeing frequent bleaching.
   
    % Get divisors for each region
    latLon = allLatLon(thisRun, 2);
    numEq = nnz(abs(latLon) <= eqLim);
    numLo = nnz(abs(latLon) > eqLim & abs(latLon) <=loLim);
    numHi = nnz(abs(latLon) > loLim);
    numTotal = numEq + numLo + numHi;
    latCounts = [numEq numLo numHi numTotal];
    
    assert(length(latLon) == numEq + numLo + numHi, 'Reefs in subsets should match total reefs in this run.');

    % Separate passes for the latitude rows. eq, lo, hi, and all
    for lat = 1:4
        subLast = lastReefBleach(abs([lastReefBleach.lat]) <= upper(lat) & abs([lastReefBleach.lat]) >= lower(lat));
        subMort = lastReefMort(abs([lastReefMort.lat]) <= upper(lat) & abs([lastReefMort.lat]) >= lower(lat));
        subHigh = anyBleach(abs([anyBleach.lat]) <= upper(lat) & abs([anyBleach.lat]) >= lower(lat));
        subMass = bleachStateMass(abs([bleachStateMass.lat]) <= upper(lat) & abs([bleachStateMass.lat]) >= lower(lat));
        subBran = bleachStateBran(abs([bleachStateBran.lat]) <= upper(lat) & abs([bleachStateBran.lat]) >= lower(lat));
        
        ct = 1; % Column in the output table.
        for n = 1:length(years)
            yr = years(n);
            percentBleached(1, ct) = yr;
            percentBleached(lat+1, ct) = 100* length(subLast([subLast.year] <= yr)) / latCounts(lat);

            percentMortality(1, ct) = yr;
            mortalityCount = length(subMort([subMort.year] <= yr));
            percentMortality(lat+1, ct) = 100* mortalityCount / latCounts(lat);
   
            highFrequency(1, ct) = yr;
            % Subset is just the last 10 years.
            % The sums count how many times a k value appears more than once
            bc = countHF(subHigh, yr);
            highFrequency(lat+1, ct) = 100 * bc  / latCounts(lat);
            
            frequentLive(1, ct) = yr;
            live = latCounts(lat) - mortalityCount;
            if live
                frequentLive(lat+1, ct) = 100 * bc / live;
            else
                frequentLive(lat+1, ct) = nan;
            end

            unrecovered(1, ct) = yr;
            % subset down to what matters
            subsetM = subMass([subMass.year] <= yr);
            subsetB = subBran([subBran.year] <= yr);
            % countBleachedR is expensive, and is called once for each
            % table cell. Since only the last entries for each reef are
            % used in any one call, it could be worth doing that filtering
            % between loops, possibly putting the years loop outside of the
            % latitudes.
            bleachedCt = countBleachedR(subsetM, subsetB);

            unrecovered(lat+1, ct) = 100 * bleachedCt / latCounts(lat);

            ct = ct + 1;
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
    labels{4, 3} = max(latLon(:));

    % Combine the last 3 arrays to indicate any form of bleaching or
    % mortality.  Summing percentages is frowned upon, but the divisors
    % are the same in each case.
    allBleaching = unrecovered;
    for i = 2:5
        for j = 1:ct-1
            allBleaching(i,j) = allBleaching(i,j) + highFrequency(i,j) ...
                + percentMortality(i,j);
        end
    end
    
    logTwo('Permanently bleached reefs as of the date given:\n');
    printTable(labels, percentBleached, length(years));

    logTwo('\nPermanent mortality as of the date given:\n');
    printTable(labels, percentMortality, length(years));

    logTwo('\nPercentage of reefs with more than one bleaching event in the previous 10 years.\n');
    printTable(labels, highFrequency, length(years));
    
    logTwo('\nPercentage of LIVING reefs with more than one bleaching event in the previous 10 years.\n');
    printTable(labels, frequentLive, length(years));

    logTwo('\nPercentage of reefs with at least one coral type in an unrecovered state.\n');
    printTable(labels, unrecovered, length(years));

    logTwo('\nPercentage of reefs with any form of current bleaching or mortality.\n');
    printTable(labels, allBleaching, length(years));
    
    % Data for plotting bleaching histories.
    % Instead of creating the plot here, save the data for use in an
    % outside program which can compare different runs.
    xForPlot = years;
    yForPlot = allBleaching(5, :);
    yEq = allBleaching(2, :);
    yLo = allBleaching(3, :);
    yHi = allBleaching(4, :);
    save(strcat(outputPath, 'bleaching/BleachingHistory', RCP, 'E=', num2str(E), 'OA=', num2str(OA), '.mat'), 'xForPlot', 'yForPlot', 'yEq', 'yLo', 'yHi', 'bleachParams');
end


function [bc] = countHF(events, yr)
    subset =  events(([events.year] <= yr & [events.year] >= yr-9));
    if isempty(subset)
        bc = 0;
    else
        % accumarray returns an array from 1 to the number of the highest-numbered reef
        % found.
        % bc is the number of reefs with more than one bleaching
        % event in this time window.
        bc = sum(accumarray([subset.k]', 1)>1);
    end
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
