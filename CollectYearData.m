% First bleaching times
% Build a table, one row per reef, with the following information.
% 1 reef number
% 2 reef lon
% 3 reef lat
% 4 first massive bleach date
% 5 first branching bleach date
% 6 first bleaching of both at once
% 7 5-year mortality end date
% 8 5-year mortality start date
% 9 delta from 6 to 2035
% 10 delta from 8 to 2035
%
% Going for readability over speed - this should be used rarely.
% Must be called after a run for access to data.
yd = NaN(reefsThisRun, 10); % preallocate
allB = bEvents(strcmp({bEvents.event}, 'BLEACH'));
lm = longMortYears; % Only available after variable-year supersymbiont introduction
lm(lm == 0) = NaN;
% Try sorting events by reef.  Speed wasn't supposed to be important but
% time drops from 13.168 seconds to 0.353.
oldK = 0;
for ev = allB
    k = [ev.k];
    if k ~= oldK
        oldK = k;
        n = 1;
    end
    % only need these two fields
    byReef{k}(n).year = ev.year;
    byReef{k}(n).coral = ev.coral;
    n = n + 1;
end
for k = 1:reefsThisRun
    yd(k, 1) = toDo(k);
    yd(k, 2) = Reefs_latlon(k, 1);
    yd(k, 3) = Reefs_latlon(k, 2);
    %b = allB([allB.k] == k);    % any kind of bleaching for reef k
    b = byReef{k};
    if isempty(b)
        mb = [];
        bb = [];
        rb = [];
    else
        mb = b(strcmp({b.coral}, 'MASS')); % massive
        bb = b(strcmp({b.coral}, 'MASS')); % branching
        rb = b(strcmp({b.coral}, 'MASS')); % reef (both)
    end
    % These should arrive in sorted order, but use min to be safe.
    if ~isempty(mb)
        yd(k, 4) = min([mb.year]);
    end
    if ~isempty(bb)
        yd(k, 5) = min([bb.year]);
    end
    if ~isempty(rb)
        yd(k, 6) = min([rb.year]);
    end
    yd(k, 7) = lm(k);  
    yd(k, 8) = yd(k, 7) - 5;
    delta = 2035 - yd(k, 6);
    delta(delta < 0) = NaN;
    yd(k, 9) = delta;
    delta = 2035 - yd(k, 8);
    delta(delta < 0) = NaN;
    yd(k, 10) = delta;
    if mod(k, 250) == 0
        fprintf('%d... ', reefsThisRun - k);
    end
end
fprintf('\n');


