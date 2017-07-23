%%
%   Identify when bleaching and recovery happens for a reef.  The
%   definitions for this purpose are:  
%   Bleaching: The ratio of symbiont density to coral density decreases to
%   10% or less of the value one year before for massive corals, or to 20%
%   or less for branching corals.
%   Recovery: The ratio increases to at least 60% of the pre-bleaching
%   value.
%   Only a bleached reef can recover, and a bleached reef can't bleach again
%   without recovering first.
%
%   Notes:
%   1) Mortality should probably be redefined now that we have a better
%      criterion for bleaching.
%   2) If the symbiont/coral ratio recovers, but the absolute amount of
%      coral is much less, is that still considered recovery?  Should there
%      be some limit?
%   3) Does it make sense to use the optimization script to choose
%      bleaching and recovery thresholds which help match expected values?
%   4) Rebuild the output plots and other summary outputs?
%%
function [bEvents, multiPlot, dCov] = ...
    Get_bleach_freq(C, C_seed, S, S_seed, k, time, latlon, TIME, dt, maxReefs, ...
                    initIndex, startYear, b, multiPlot)
    persistent sBleachSum cBleachSum bBleachSum bMortSum;
    persistent sRecSum massRecSum seedRecSSum seedRecCSum massActSum;
    
    if k == 1
        sBleachSum = 0;
        cBleachSum = 0;
        bBleachSum = 0;
        bMortSum = 0;
        sRecSum = 0;
       	massRecSum = 0;
       	seedRecSSum = 0;
       	seedRecCSum = 0;
        massActSum = 0;
    end

    %% Constants and options
    MASS = 1;
    BRAN = 2;
    typeName = {'MASS', 'BRAN'};
    
    %% Put all adjustable parameters in a structure.
    % Massive before branching    
    % Define bleaching and recovery
    %bleachFraction = [0.12 0.12];
    sBleach = b.sBleach;
    cBleach = b.cBleach;
    seedThresh = C_seed .* b.cSeedThresholdMult;
% out of use    seedThreshRecover = C_seed .* b.cSeedRecoverMult;
    %sRecoverFraction = b.sRecoverFraction;
    % out of use cRecoverFraction = b.cRecoverFraction;
    yearsToMortality = b.yearsToMortality;
    yearsAverage = b.yearsRunningAverage;
    sRecoverySeedMult = b.sRecoverySeedMult;
    cRecoverySeedMult = b.cRecoverySeedMult;

    %requireCoralRecovery = false;

    preBleachSymbiont = [0.0 0.0];
    % out of use preBleachCoral = [0.0 0.0];
    
    % fprintf('Sym bleach: %6.3f %6.3f  Coral bleach: %6.3f %6.3f \n', sBleach, cBleach);
    
    
    %% Record each bleaching or recovery event
    eventCount = 0;  % index into mEvent array
    bleachCount = 0;
    bleachCount8510 = 0;
    % Initialize to a larger than average size so there's not much array
    % growth.
    bEvents(32).last = 0; 
       
    stepsPerYear = 12 / dt;
    assert(length(C)*dt == length(TIME), 'Expected months in TIME to equal timesteps in C times dt.');
    
    % For both S and C we will compare to a trailing average, but it's a
    % per-date average so seasonal variation is preserved.  This started as
    % a 5-year average, but it's tuneable.
    %cCompare = trailingSameDateAverage(C, yearsAverage, stepsPerYear);
    %sCompare = trailingSameDateAverage(S, yearsAverage, stepsPerYear);
    
    % Find the minimum coral cover and symbiont density for each year.
    [crow, ccol] = size(C);
    [srow, scol] = size(S);
    
    %% WARNING: Note that S is collapsed to one type of symbiont for each
    %  coral type.  That's good in this routine since any alternate
    %  symbiont is considered just as good for recovery.  Do NOT pass this
    %  version of S back to other parts of the code where all columns may
    %  be needed.
    if multiPlot.active && any(k == multiPlot.reefs)
        % If S will be plotted, preserve a copy before the change.
        SPlot = S;
    end
    if scol > 2
        S(:, 1) = sum(S(:, 1:2:end), 2);
        S(:, 2) = sum(S(:, 2:2:end), 2);
        scol = 2;
    end

        

    Cmin(crow, ccol) = 0; % preallocate Cmin
    Smin(srow, scol) = 0; % preallocate Smin
    % Reshape each column of C into an array with one year (96 time steps)
    % per row, so the "min" calculation is simple.
    %{
    % initIndex is not an exact multiple of the time step, but needs to be
    % here.
    idx = stepsPerYear*ceil(initIndex/stepsPerYear);
    for i = 1:ccol
        Cmin(:, i) = min(reshape(C(1:idx, i), stepsPerYear, []))';
    end
    averageCMin = mean(Cmin); % to get a single global minimum
        for j = 1:scol
        Smin(:, j) = min(reshape(S(1:idx, j), stepsPerYear, []))';
    end
    averageSMin = mean(Smin);
    %}

    % Make copies of C and S, but with the annual mimimum for each column
    % stored in each row.
    for j = 1:ccol
        for i = 1:stepsPerYear:crow
            Cmin(i:i+stepsPerYear-1, j) = min(C(i:i+stepsPerYear-1, j));
        end
    end
    for j = 1:scol
        for i = 1:stepsPerYear:crow
            Smin(i:i+stepsPerYear-1, j) = min(S(i:i+stepsPerYear-1, j));
        end
    end
    % Average over recent years, to increase capture of gradual drops.
    if yearsAverage > 1
        Cmin = trailingSameDateAverage(Cmin, yearsAverage, stepsPerYear);
        Smin = trailingSameDateAverage(Smin, yearsAverage, stepsPerYear);
    end

    % Just for diagnostics, at least for now:
    if multiPlot.active
        bMark = zeros(length(S), 2);
        rMark = zeros(length(S), 2);
        mMark = zeros(length(S), 2);
    end
    % Flags for when corals are considered bleached and dead.
    bState = zeros(length(S), 4);
    mState = zeros(length(S), 6);
    % Mortality events are plotted from the main program, originally using
    % an array called dCov from the old Get_mort_freq function.
    % Duplicate that here even though it's redundant.
    % By starting with nan's, only values set to 1 will plot.
    %dCov = nan(length(S), 2);
    dCov(length(S), 2) = 0;
    
    % Columns
    % For 1 and 2, use existing MASS and BRAN values.
    BOTH = 3;
    EITHER = 4;
    START = 5;
    END = 6;
       
    sBleachCount = 0;
    cBleachCount = 0;
    bBleachCount = 0; 
    bMortCount = 0;
    sRecCount = 0;
    massRecCount = 0;
    seedRecSCount = 0;
    seedRecCCount = 0;
    massActCount = 0;
    % if massive coral has mortality due to seed threshold (not bleaching),
    % keep track and allow it to recover the same way.
    massiveSeedMort = false;
    for cType = MASS:BRAN
        bleached = false;
        anyBleached = false;
        lastBleach = 0;
        lastBleachEvent = 0;
        mort = false;
        %fprintf('Checking for bleaching, type %f \n', cType);
        % Big change 1/9/17 - just look at annual values, not every time
        % step.
        for i = floor(1.5*stepsPerYear):stepsPerYear:length(S(:,cType))

            if bleached
                % Require 60% of previous S value AND > 5*seed for massive
                % symbionts and coral.
                % Optionally, required coral to recover to 30% of previous
                % value as well (requireCoralRecovery).
                % Note that S_seed is a single constant, C_seed has two
                % values.

                % Note: requireCoralRecovery is off, but implementation
                % below is really "allow coral recovery". Fix naming or
                % logic if enabled.
                %sRec = Smin(i, cType) > preBleachSymbiont(cType) * sRecoverFraction(cType);
                % 1/10/17 TEST XXX sRec = Smin(i, cType) > Smin(i-stepsPerYear, cType) * sRecoverFraction(cType);
                sRec = true;
                %cRec = requireCoralRecovery && Cmin(i, cType) > preBleachCoral(cType) * cRecoverFraction(cType);
             
                massRec = massiveSeedMort && cType == MASS;
                seedRecS = Smin(i, cType) > sRecoverySeedMult(cType)*S_seed(cType);
                % debugC = Cmin(i, cType);
                seedRecC = Cmin(i, cType) > cRecoverySeedMult(cType)*C_seed(cType);
                if seedRecC && ((seedRecS && sRec ) || massRec)

                    sRecCount = sRecCount + sRec;
                    massRecCount = massRecCount + massRec;
                    seedRecSCount = seedRecSCount + seedRecS;
                    seedRecCCount = seedRecCCount + seedRecC;
                    massActCount = massActCount + (~(seedRecS && sRec ) && massRec);

                    bleached = false;
                    massiveSeedMort = false;
                    bState(i:end, cType) = 0;
                    if multiPlot.active; rMark(i, cType) = 1; end;
                    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i), typeName(cType), 'RECOVER', k, latlon);
                    if mort
                        mort = false;
                        [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i), typeName(cType), 'REVIVE', k, latlon);
                        mState(i:end, cType) = 0;
                    end
                else
                    % If there's no recovery, there may also be mortality
                    % (5 years bleached).  Here it is per coral type -
                    %  overall mortality will be deduced in New_Stats_Bleach
                    if ~mort
                        if (time(i) - lastBleach) >= yearsToMortality*365
                            mort = true;
                            mMark(i, cType) = 0.5;
                            if multiPlot.active; mMark(i, cType) = 1; end;
                            
                            [bEvents, eventCount] = ...
                                makeBEvent(bEvents, eventCount, time(i), typeName(cType), 'MORTALITY', k, latlon);
                            dCov(i, cType) = 1;
                            mState(i:end, cType) = 1;
                        end
                    end
                end
            else
                % Not bleached.  Check for bleaching.
                % Line 1: drops in symbionts
                % Line 2: drops in coral cover
                sB = Smin(i, cType) < Smin(i-stepsPerYear, cType)*sBleach(cType);
                cB = Cmin(i, cType) < Cmin(i-stepsPerYear, cType)*cBleach(cType);
                if sB || cB
                    bBleachCount = bBleachCount + (sB && cB);
                    sBleachCount = sBleachCount + (sB && ~cB);
                    cBleachCount = cBleachCount + (cB && ~sB);
                  
                    bleached = true;
                    bState(i:end, cType) = 1;
                    anyBleached = true;
                    lastBleach = time(i);
                    if reefYear(lastBleach) >= 1985 && reefYear(lastBleach) <= 2010
                        bleachCount8510 = bleachCount8510 + 1;
                    end
                    bleachCount = bleachCount + 1;
                    if multiPlot.active; bMark(i, cType) = 1; end;
                    preBleachSymbiont(cType) = Smin(i-stepsPerYear, cType);
                    % out of use preBleachCoral(cType) = Cmin(i-stepsPerYear, cType);
                    %fprintf('Reef %d, S = %0.2e, year = %d, %s bleached\n', ...
                    %    k, S(i, cType), reefYear(time(i)), typeName{cType});
                    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i), typeName(cType), 'BLEACH', k, latlon);
                    lastBleachEvent = eventCount;
                end
            end
            % Most of this "for" is about bleaching, but also check for
            % mortality based on an extreme low value of coral cover.
            % Note that we also consider that type of coral bleached, so we
            % look for recovery, not bleaching.
%HERE: mort interaction with the bleaching logic needs work
%    - for example, if we check for bleaching recovery after coral-based mortality, what's the old symbiont value?  Skip that part?
            if ~mort
                % Massives must maintain 5*Seed, and branching 30*Seed
                if C(i, cType) < seedThresh(cType)
                    %fprintf('Declaring low cover mortality for reef %d at i = %d, in %d based on %s coral.\n', k, i, reefYear(time(i)), typeName{cType});
                    mort = true;
                    bMortCount = bMortCount + 1;
                    if cType == MASS
                        massiveSeedMort = true;
                    end
                    % mort implies bleached.  The if just helps with
                    % debugging.
                    if ~bleached
                        bleached = true;
                    end
                    if multiPlot.active; mMark(i, cType) = 1; end;

                    [bEvents, eventCount] = ...
                        makeBEvent(bEvents, eventCount, time(i), typeName(cType), 'MORTALITY', k, latlon);
                    mState(i:end, cType) = 1;
                    bState(i:end, cType) = 1;
                    dCov(i, cType) = 1;

                    % Set a higher-than actual value, since we want to be
                    % well above the current level to declare recovery.
                    % out of use preBleachCoral(cType) = seedThreshRecover(cType);
                    
                    % Not setting anyBleached since this is not bleaching
                    % in the normal sense.
                    % anyBleached = true;
                end
            end
        end  % end time loop
        
        % The last bleaching event for each coral type should be flagged if
        % there was no subsequent recovery.
        if anyBleached
            bEvents(lastBleachEvent).last = 1;
        end
    end  % end coral type loop
    
    % Store how many times bleaching occurred for direct use in mapping.
    % This is easier than adding events later, and includes maps with zero
    % bleaching.
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'REEF', 'BLEACHCOUNT', k, latlon, bleachCount);
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'REEF', 'BLEACHCOUNT8510', k, latlon, bleachCount8510);

    % Diagnostic for how often coral and symbionts trigger bleaching.
    %fprintf('For reef %4d, bleaching triggered by S/C/Both %2d %2d %2d Or coral mort %2d\n', k, sBleachCount, cBleachCount, bBleachCount, bMortCount);
    sBleachSum = sBleachSum + sBleachCount;
    cBleachSum = cBleachSum + cBleachCount;
    bBleachSum = bBleachSum + bBleachCount; 
    bMortSum = bMortSum + bMortCount;

    sRecSum  = sRecSum + sRecCount;
    massRecSum = massRecSum + massRecCount;
    seedRecSSum = seedRecSSum + seedRecSCount;
    seedRecCSum = seedRecCSum + seedRecCCount;
    massActSum = massActSum + massActCount;
    if k == maxReefs
        fprintf('Summary of bleaching and mortality triggers:\n');
        fprintf('Bleaching was triggered by S/C/Both %4d %4d %4d times and coral mortality %4d times\n', sBleachSum, cBleachSum, bBleachSum, bMortSum);
        triggerSum = (sBleachSum + cBleachSum + bBleachSum + bMortSum)/100;
        fprintf('By S/C/Both percent %5.2f%% %5.2f%% %5.2f%% and coral mortality %5.2f%% \n', sBleachSum/triggerSum, cBleachSum/triggerSum, bBleachSum/triggerSum, bMortSum/triggerSum);
	
        fprintf('Summary of recovery triggers:\n');
        fprintf('Recovery was flagged by symbiont recovery and massive coral recovery %4d and %4d times and by symbiont and coral seed requirements %4d and %4d times\n', sRecSum, massRecSum, seedRecSSum, seedRecCSum);
        fprintf('Recovery was dependent on massive coral recovery %4d times.\n', massActSum);
    end


    % Check for mortality of both coral types
    mState(:, BOTH) = mState(:, MASS) & mState(:, BRAN);
    mState(:, EITHER) = mState(:, MASS) | mState(:, BRAN);

    % Start and end reef mort
    mState(2:end, START) = mState(2:end, BOTH) > mState(1:end-1,BOTH);
    mState(2:end, END) = mState(2:end, BOTH) < mState(1:end-1,BOTH);
    % Save events for full-reef mortality
    mS = find(mState(:, START) > 0);
    mE = find(mState(:, END) > 0);
    lastRevival = 0;
    if ~isempty(mE)
        for i = mE'
            [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i), 'REEF', 'REVIVE', k, latlon);
        end
        lastRevival = i;
    end
    if ~isempty(mS)
        for i = mS'
            [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i), 'REEF', 'MORTALITY', k, latlon);
        end
        % If the last mortality is after the last revival (or there is no revival), mark it as
        % permanent.
        if i > lastRevival
            bEvents(eventCount).last = 1;
        end
    end
    
    %% New special entry type, 2/18/2017
    %  We want to know when a reef is in a state of mortality for 5 years.
    %  For each mortality event, see if mortality ends within 5 years.
    %  If not, note the year and ignore the rest of the array.  If we never
    %  find a matching case, store an event anyway, with the year set to
    %  zero.
    span = stepsPerYear * 5;
    if isempty(mS)
        % No mortality, mark with zero time.
        [bEvents, eventCount] = makeBEvent(bEvents, eventCount, 0, 'REEF', 'LONGMORT', k, latlon);
    else
        found = false;
        for i = mS'
            if ~found
                % Mortality started.  Either it ends within 5 years (keep
                % looking), never ends (save this time+5 years), ends in more
                % than 5 years (save), or the run ends in less than 5 years
                % (mark with zero).
                if isempty(mE) || isempty(find(mE > i, 1))
                    if i+span > length(time)
                        % never in mortality for 5 years.
                        [bEvents, eventCount] = makeBEvent(bEvents, eventCount, 0, 'REEF', 'LONGMORT', k, latlon);
                    else 
                        [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i+span), 'REEF', 'LONGMORT', k, latlon);
                    end
                    found = true;
                else
                    % If there's an end in mE within 5 years, this is not what
                    % we're looking for.  
                    nextEnd = find(mE > i, 1);  % first index in the array of indexes
                    if mE(nextEnd) > i + span
                        [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(i+span), 'REEF', 'LONGMORT', k, latlon);
                        found = true;
                    end
                end
            end
        end
        if ~found
            [bEvents, eventCount] = makeBEvent(bEvents, eventCount, 0, 'REEF', 'LONGMORT', k, latlon);
        end

    end

    
    % Bleaching state
    bState(:, BOTH) = bState(:, MASS) & bState(:, BRAN);
    bState(:, EITHER) = bState(:, MASS) | bState(:, BRAN);
    
    % Answer the question: For how many years was this reef bleached during
    % 1985 to 2010?
    % The bState array has a one or zero for each timestep.  We want to
    % count each year as either in or out.  Or do we?  We could also count
    % the number of flagged steps and round to years...
    i1985 = 1+(1985-startYear)*stepsPerYear; % beginning of 1985
    i2010 = (2011-startYear)*stepsPerYear;   % end of 2010
    bYears(MASS) = sum(bState(i1985:i2010, MASS))/stepsPerYear;
    bYears(BRAN) = sum(bState(i1985:i2010, BRAN))/stepsPerYear;
    bYears(BOTH) = sum(bState(i1985:i2010, BOTH))/stepsPerYear;
    bYears(EITHER) = sum(bState(i1985:i2010, EITHER))/stepsPerYear;
    % Special events which act as summaries
    % Note that these are NOT in sync with the current 85_10 bleaching stat
    % which is targeted at 2%.  That use the number of events, ignoring
    % their duration.
    % OBSOLETE?
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'REEF', 'BLEACH8510', k, latlon, bYears(BOTH));
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'BRAN', 'BLEACH8510', k, latlon, bYears(BRAN));
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'MASS', 'BLEACH8510', k, latlon, bYears(MASS));
    [bEvents, eventCount] = makeBEvent(bEvents, eventCount, time(1), 'EITHER', 'BLEACH8510', k, latlon, bYears(EITHER));
    
    % Delete unused event storage locations before returning.
    bEvents(eventCount+1:end) = [];
    
    % yearAgoFraction = S ./ circshift(S, stepsPerYear, 1);
    
     
    if multiPlot.active && any(k == multiPlot.reefs)
        % Reset to a new figure if "panels" have been used on the plot
        if multiPlot.count == multiPlot.panels
            multiPlot.figure = multiPlot.figure + (~multiPlot.print);
            multiPlot.count = 0;
        end
        
        %cType = dType;    
        startYr = multiPlot.startYr;  % for the plot - overall start of data is at startYear
        endYr = multiPlot.endYr; % was 1920-2050
        startComp = max(1, (startYr-startYear)*stepsPerYear);
        endComp = (endYr-startYear)*stepsPerYear;

        % Make a single figure instead of one per coral type:
        % TODO - this uses the second "div" from above.  Compute a single
        % best value.
        %disp('In Get_bleach_freq plot.');
        % Scale both symbionts to largest coral value
        % Max values for each type:
        maxC = max(C(startComp:endComp, :));
        maxS = max(SPlot(startComp:endComp, :));
        % single value: div = max(maxS) / max(maxC);
        % separate symbiont scales, based on larger of the two coral
        % values
        maxC = max(maxC);
        % 1/30/17 keep all symbionts on the same scale
        maxS = max(maxS);
        div = maxS ./ maxC;
        div = round(1+log2(div));
        div = round(2 .^ div, 1, 'significant');
        SPlot = SPlot ./ div;

            %fprintf('maxC = %0.2e, maxS = %0.2e, %0.2e, div = %d, %d\n', maxC, maxS, div);
        titleText = sprintf('Reef %d, symbiont based', k);
        multiPlot.count = multiPlot.count + 1;
        % Note that only bEvents and multiPlot are returned from this
        % function.  At this point it's okay to modify values for plotting
        % purposes only.
        mState(:, BOTH) = mState(:, MASS) & mState(:, BRAN);
        bMark = bMark - (bMark == 0);
        rMark = rMark - (rMark == 0);
        mMark = mMark - (mMark == 0);

        % Note hardwired "true" as 2nd argument forces output each time.
        % Starting 2/17/17, flag each set of y values with true for a 
        % line on the left axis and false for symbols scaled right.
        graphCompare(multiPlot, true, k, titleText, '', ...
            ((startComp:endComp)/stepsPerYear)+1860, 'index', ...
            C(startComp:endComp, MASS), 'Massive Coral' , true, ...
            C(startComp:endComp, BRAN), 'Branching Coral' , true, ...
            SPlot(startComp:endComp, MASS), strcat('Mass Symb / ', num2str(div)), true, ...
            SPlot(startComp:endComp, BRAN), strcat('Bran Symb / ', num2str(div)), true, ...
            SPlot(startComp:endComp, MASS+2), strcat('Alt Mass Symb / ', num2str(div)), true, ...
            SPlot(startComp:endComp, BRAN+2), strcat('Alt Bran Symb / ', num2str(div)), true, ...
            linspace(5000000, 5000000, endComp-startComp+1), 'Massive seed x5', true, ...
            linspace(500000, 500000, endComp-startComp+1), 'Branching seed x5', true, ...
            bMark(startComp:endComp, MASS), 'Mass Bleaching', false, ...
            rMark(startComp:endComp, MASS), 'Mass Recovery', false, ...
            mMark(startComp:endComp, MASS), 'Mass Mortality', false, ...
            bMark(startComp:endComp, BRAN), 'Bran Bleaching', false, ...
            rMark(startComp:endComp, BRAN), 'Bran Recovery', false, ...
            mMark(startComp:endComp, BRAN), 'Bran Mortality', false ...           
          );
    end

end

%% Convert a timestamp (in days) to a whole number year.  If zero is provide,
%  zero is returned.
function [thisYear] = reefYear(d)
    % Note that these bases are just a starting point for the calculation.
    % They do not need to match the first date of the run, but they must be
    % that date or earlier.
    %baseDate = 679352;  % January 1, 1860
    %baseYear = 1860;
    if isempty(d)
        thisYear = [];
        return;
    end
    thisYear = floor(1860 + (d-679352)/365.25);  % Not perfect over long spans, but okay for the model.
end

%% Add one entry to the bEvents list.
function [bEvents, eventCount] = makeBEvent(bEvents, eventCount, t, type, event, k, latlon, yearCount)
    eventCount = eventCount + 1;
    bEvents(eventCount).year = reefYear(t);
    bEvents(eventCount).coral = char(type);
    bEvents(eventCount).event = event;
    bEvents(eventCount).last = 0;  % Override the last one to 1 outside loop
    bEvents(eventCount).k = k;
    bEvents(eventCount).lat = latlon(2);
    bEvents(eventCount).lon = latlon(1);
    if exist('yearCount','var')
        bEvents(eventCount).count = yearCount;
    else
        bEvents(eventCount).count = 0;
    end
end

%{
Debug info:
i increments by 96 per year.
Matlab date number increments by one unit per day
i       year  Matlab date number
 1000   1871    683530
 2000   1881
13008   1996
14000   2006
14100   2007
14300   2009
14800   2015
15000   2017
15200   2019
15300   2020
16080   2028
23040   2101    767390  (Jan 10)

%}
