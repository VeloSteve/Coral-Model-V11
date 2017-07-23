%% Repeatedly run the model for all cases at once.
%  Variables from the supported set have defaults in the main code, but
%  it's best to set them all here to avoid surprises.

% First those we aren't changing in this run:
scriptVars.superMode = 0;
scriptVars.superAdvantage = 0.0;

%rcpList = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
%rcpList = {'rcp26', 'rcp45', 'rcp60'};
rcpList = {'rcp26', 'rcp45'}; %, 'rcp60', 'rcp85'};

% Now all the cases
for ooo = 0:1
    for eee = 0  %0:1
        for rrr = rcpList
            scriptVars.E = eee;
            scriptVars.OA = ooo;
            scriptVars.RCP = rrr{1};
            fprintf('Starting model with E = %d, OA = %d, RCP %s\n', eee, ooo, rrr{1});
            A_Coral_Model_170118
        end
    end
end
disp('All cases are complete.');
% If this isn't done, the next one-off model run will mysteriously pick
% up values from the last iteration here!
clearvars scriptVars

