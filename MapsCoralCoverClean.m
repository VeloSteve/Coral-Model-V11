%% Make Maps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% last updated: 1-6-15                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = MapsCoralCoverClean(fullDir, Reefs_latlon, activeReefs, ...
    lastYearAlive, events85_2010, eventsAllYears, yearRange, fullYearRange, ...
    modelChoices, filePrefix )
% Add paths and load mortality statistics
%load(strcat('~/Dropbox/Matlab/SymbiontGenetics/',filename,'/201616_testNF_1925reefs.mat'),'Mort_stats')
format shortg;
% filename = '201616_figs'; %filename = strcat(dateString,'_figs'); mkdir(filename); % location to save files
% map %% NOTE: worldmap doesn't seem to be working on work computer

%% Make map of last mortality event recorded

% We need to map a spot for all reefs, to show those that never bleached.
% Not every reef has a last mortality, but all have BLEACH8510 stats.
allReefs(1:length(activeReefs), 1) = Reefs_latlon(activeReefs, 1);
allReefs(1:length(activeReefs), 2) = Reefs_latlon(activeReefs, 2);

tName = strcat(modelChoices,'. Year Corals No Longer Persist');
fileBase = strcat(fullDir, filePrefix,'_LastYrCoralMap');
% Green points everywhere
oneMap(12, allReefs(:, 1), allReefs(:, 2), [0 0.8 0], yearRange, parula, tName,'', false);

% Color-scaled points where there is a last year
outFile = strcat(fileBase, '.pdf');
if any(lastYearAlive)
    ind = find(lastYearAlive);
    customColors = customScale();
    oneMap(12, Reefs_latlon(ind, 1), Reefs_latlon(ind, 2), lastYearAlive(ind), yearRange, customColors, tName, outFile, true);
end

% This one may be post-processed, so save .fig
if verLessThan('matlab', '8.2')
    saveas(gcf, fileBase, 'fig');
else
    savefig(strcat(fileBase,'.fig'));
end

%% Make map showing # all bleaching events bn 1985-2010
tName = 'Bleaching Events Between 1985-2010';
outFile = strcat(fullDir, filePrefix,'_MortEvents8510Map','.pdf');
oneMap(13, allReefs(:, 1), allReefs(:, 2), events85_2010(activeReefs), [], jet, tName, outFile, false);


%% Figure 14 Make map showing # all bleaching events
rangeText = sprintf('%d-%d',fullYearRange);
tName = strcat('Bleaching Events Between ', rangeText);
outFile = strcat(fullDir, filePrefix,'_AllMortEventsMap','.pdf');
oneMap(14, allReefs(:, 1), allReefs(:, 2), eventsAllYears(activeReefs), [], jet, tName, outFile, false);


%% Figure 15  Same as 14 but with restricted color scale
cRange = [0, 20];
outFile = strcat(fullDir, filePrefix,'_AllMortEventsMap_Scale20','.pdf');
oneMap(15, allReefs(:, 1), allReefs(:, 2), eventsAllYears(activeReefs), cRange, jet, tName, outFile, false);


%% Figure 16  Same as 14 but with log2 of the number of events
%{
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
tName = 'Bleaching Events Between 1861-2100 (log base 2)';
outFile = strcat(fullDir, filePrefix, '_AllMortEventsMap_log2', '.pdf');
oneMap(16, allReefs(:, 1), allReefs(:, 2), log2(eventsAllYears(activeReefs)), [], jet, tName, outFile, false);
%}




end

% Arguments:
% n     figure number
% lons  longitudes (was [events.lon])
% lats  latitudes
% values what to plot at each position
% cRange man and max data values for color scale
% t title
% outFile pdf output file
function [] = oneMap(n, lons, lats, values, cRange, cMap, t, outFile, add)
    figure(n);
    if add
        hold on;
    else
        clf;
        % first pass only:
        m_proj('miller'); % , 'lon', 155.0); - offsets map, but drops some data!
        m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
        m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
    end

    % Points with last-year mortality values:
    [LONG,LAT] = m_ll2xy(lons,lats); hold on % convert reef points to M-Map lat long

    scatter(LONG,LAT,5, values) ; % plot bleaching events onto map
    
    if isempty(cMap)
        colormap default;
    else
        colormap(cMap)
    end
    if ~isempty(cRange)
        caxis(cRange);
    end
    colorbar
    title(t)
    
    if ~isempty(outFile)
        print('-dpdf', '-r200', outFile);
    end
    hold off;
end

function [cMap] = customScale()
    % Map code posted by Stephen Cobeldick at https://www.mathworks.com/matlabcentral/fileexchange/25536-red-blue-colormap
    m = 20;
    n = fix(m/2);
    x = n~=(m/2); 
    b = [(0:1:n-1)/n,ones(1,n+x)]; 
    g = [(0:1:n-1)/n/2,ones(1,x),(n-1:-1:0)/n/2]; % Extra "/2" so light blue doesn't fade into ocean blue
    r = [ones(1,n+x),(n-1:-1:0)/n]; 
    cMap = [r(:),g(:),b(:)];
end