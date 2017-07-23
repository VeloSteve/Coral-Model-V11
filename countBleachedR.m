%% Call a recursive function for each of two reefs subsets and count the
%  number of reefs in a bleached state.
function [count] = countBleachedR(subsetM, subsetB)
    a = listBleached(subsetM);
    b = listBleached(subsetB);
    count = length(unique([a b]));
end

%% Find a list of the reefs which are bleached in a given subset of events,
%  typically constrained to some range of years.  If the subset is large
%  (> 15) and contains more than one reef, split it recursively.
%  Bleached means that there has been bleaching NOT followed by a recovery,
%  revival or mortality. Note that EITHER coral type bleached is counted as
%  reef bleaching here.
function [bleachList] = listBleached(subset)
    big = 15;  % values from 15 to 25 give fastest results.
    bleachList = [];
    % Note that Matlab complains about comparing length to zero, but the
    % result of isempty is not the same in this case!
    if isempty(subset) || length(subset) == 0
        return;
    end
    reefs = unique([subset.k]);

    if length(reefs) > big
        maxK = max(reefs);
        minK = min(reefs);
        if maxK > minK
            half = floor((maxK + minK)/2);
                a = listBleached(subset([subset.k] <= half));
                b = listBleached(subset([subset.k] >  half));
                bleachList = unique([a b]);
            return;
        end
    end
    % fprintf('Recursive check of %d reefs starting with k = %d.\n', length(reefs), min(reefs));
    % For each reef look at the various events to see if it's considered
    % bleached.
    for k = reefs
        % Get just one reef's events.  This uses the most time, in 13462 calls with everyx=5 in one example.
        last = subset([subset.k] == k);
        % Latest year with an event.
        maxY = max([last.year]);
        % Subset to only events in that last year.
        last = last([last.year] == maxY);

        if length(last) > 1
            % If all are non-bleaching or all bleaching, just keep one
            % If not, it's a problem
            state = 0; % 0 = start, 1 = bleach, 2 = recover
            for i = 1:length(last)
                e = last(i);
                if strcmp({e.event}, 'BLEACH') || strcmp({e.event}, 'MORTALITY')
                    if state == 0 || state == 1
                        state = 1;
                    elseif state == 2
                        %error('Recovery and bleaching in the same year!');
                        fprintf('ERROR - recovery and bleaching for reef %d\n', k);
                    end
                else
                    if state == 0 || state == 2
                        state = 2;
                    elseif state == 1
                        %error('Recovery and bleaching in the same year.');
                        fprintf('ERROR -- recovery and bleaching for reef %d\n', k);

                    end
                end
            end
            % if we get here, all events in last are effectively the same.
            last = last(1);
        end
        if ~isempty(last) && length(last) > 0
            ev = {last.event};
            if strcmp(ev, 'REVIVE') || strcmp(ev, 'RECOVER') || strcmp(ev, 'MORTALITY')
                % not included
            elseif strcmp(ev, 'BLEACH')
                bleachList(end+1) = last.k;
            else
                error('Event type should have already been checked.');
            end
        end
    end
end