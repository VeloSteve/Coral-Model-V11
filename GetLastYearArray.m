%% Make Maps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% last updated: 1-6-15                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lastYear] = GetLastYearArray(bEvents, maxReefs)

% Get last full-reef mortality events:
events = bEvents([bEvents.last]==1 & strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'MORTALITY'));
lastYear = NaN(maxReefs, 1);

    for e = events
        lastYear(e.k) = e.year;
    end

end