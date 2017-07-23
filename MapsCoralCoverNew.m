%% Make Maps
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% last updated: 1-6-15                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = MapsCoralCoverNew(fullDir, bEvents, yearRange, modelChoices, filePrefix )
% Add paths and load mortality statistics
%load(strcat('~/Dropbox/Matlab/SymbiontGenetics/',filename,'/201616_testNF_1925reefs.mat'),'Mort_stats')
format shortg;
% filename = '201616_figs'; %filename = strcat(dateString,'_figs'); mkdir(filename); % location to save files
% map %% NOTE: worldmap doesn't seem to be working on work computer

%% Make map of last mortality event recorded
% Get last full-reef mortality events:
events = bEvents([bEvents.last]==1 & strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'MORTALITY'));
% We need to map a spot for all reefs, to show those that never bleached.
% Not every reef has a last mortality, but all have BLEACH8510 stats.
forBG = bEvents(strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'BLEACH8510'));

tName = strcat(modelChoices,'. Year Corals No Longer Persist');
fileBase = strcat(fullDir, filePrefix,'_LastYrCoralMap');
% Green points everywhere
oneMap(12, [forBG.lon], [forBG.lat], [0 0.8 0], yearRange, '', tName,'', false);

% Color-scaled points where there is a last year
outFile = strcat(fileBase, '.pdf');
oneMap(12, [events.lon], [events.lat], [events.year], yearRange, customScale(), tName, outFile, true);

% This one may be post-processed, so save .fig
if verLessThan('matlab', '8.2')
    saveas(gcf, fileBase, 'fig');
else
    savefig(strcat(fileBase,'.fig'));
end

%% Make map showing # all bleaching events bn 1985-2010
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT8510'));
tName = 'Bleaching Events Between 1985-2010';
outFile = strcat(fullDir, filePrefix,'_MortEvents8510Map','.pdf');
oneMap(13, [events.lon], [events.lat], [events.count], [], jet, tName, outFile, false);


%% Figure 14 Make map showing # all bleaching events
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
tName = 'Bleaching Events Between 1861-2100';
outFile = strcat(fullDir, filePrefix,'_AllMortEventsMap','.pdf');
oneMap(14, [events.lon], [events.lat], [events.count], [], jet, tName, outFile, false);


%% Figure 15  Same as 14 but with restricted color scale
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
cRange = [0, 20];
tName = 'Bleaching Events Between 1861-2100';
outFile = strcat(fullDir, filePrefix,'_AllMortEventsMap_Scale20','.pdf');
oneMap(15, [events.lon], [events.lat], [events.count], cRange, jet, tName, outFile, false);


%% Figure 16  Same as 14 but with log2 of the number of events
%{
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
tName = 'Bleaching Events Between 1861-2100 (log base 2)';
outFile = strcat(fullDir, filePrefix, '_AllMortEventsMap_log2', '.pdf');
oneMap(16, [events.lon], [events.lat], log2([events.count]), [], jet, tName, outFile, false);
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
        fprintf('Hold back on for fig %d to add %d points.\n', n, length(lons));
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