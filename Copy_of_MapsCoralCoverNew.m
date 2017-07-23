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
figure(12);
m_proj('miller'); % , 'lon', 155.0); - offsets map, but drops some data!
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
% Get last full-reef mortality events:
events = bEvents([bEvents.last]==1 & strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'MORTALITY'));
% We need to map a black spot for all reefs, to show those that never
% bleached.  Not every reef has a last mortality, but all have BLEACH8510
% stats.
forBG = bEvents(strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'BLEACH8510'));
% Points with last-year mortality values:
[LONG,LAT] = m_ll2xy([events.lon],[events.lat]); hold on % convert reef points to M-Map lat long
caxis(yearRange);
% All reefs:
[LONG2,LAT2] = m_ll2xy([forBG.lon],[forBG.lat]);  % convert reef points to M-Map lat long
scat12a = scatter(LONG2, LAT2, 5, [0 0 0]); %[.7 .7 .7]) % plot reefs  onto map
%scatter(LONG,LAT,5, (psw2_var_all(:,3)/2)+.7) ; % plot bleaching events onto map
scat12b = scatter(LONG,LAT,5, [events.year]) ; % plot bleaching events onto map
%colormap(hot); %(flipud(jet))
colorbar
title(strcat(modelChoices,'. Year Corals No Longer Persist'))
%print('-dpdf', '-r200', strcat(fullDir, filePrefix,'_LastYrCoralMap','.pdf'));
% This one may be post-processed, so save .fig
if verLessThan('matlab', '8.2')
    saveas(gcf, strcat(fullDir, filePrefix,'_LastYrCoralMap'), 'fig');
else
    savefig(strcat(fullDir, filePrefix,'_LastYrCoralMap','.fig'));
end
hold off;

%% Make map showing # all bleaching events bn 1985-2010
figure(13)
m_proj('miller');
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT8510'));
[LONG,LAT] = m_ll2xy([events.lon],[events.lat]); hold on % convert reef points to M-Map lat long
%scatter(Mort_stats(:,1),Mort_stats(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
scat13a = scatter(LONG,LAT,5, [events.count]) ; % plot bleaching events onto map
colormap(jet)
colorbar
title('Bleaching Events Between 1985-2010')
print('-dpdf', '-r200', strcat(fullDir, filePrefix,'_MortEvents8510Map','.pdf'));
hold off;

%% Make map showing # all bleaching events
figure(14)
m_proj('miller');
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
[LONG,LAT] = m_ll2xy([events.lon],[events.lat]); hold on % convert reef points to M-Map lat long
%scatter(Mort_stats(:,1),Mort_stats(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
scat14a = scatter(LONG,LAT,5, [events.count]) ; % plot bleaching events onto map
colormap(jet)
colorbar
title('Bleaching Events Between 1861-2100')
print('-dpdf', '-r200', strcat(fullDir, filePrefix,'_AllMortEventsMap','.pdf'));
hold off;

figure(15)
m_proj('miller');
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
[LONG,LAT] = m_ll2xy([events.lon],[events.lat]); hold on % convert reef points to M-Map lat long
%scatter(Mort_stats(:,1),Mort_stats(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
scat14a = scatter(LONG,LAT,5, [events.count]) ; % plot bleaching events onto map
colormap(jet)
caxis([0, 20]);
colorbar
title('Bleaching Events Between 1861-2100')
print('-dpdf', '-r200', strcat(fullDir, filePrefix,'_AllMortEventsMap_Scale20','.pdf'));
hold off;

figure(16)
m_proj('miller');
m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
events = bEvents(strcmp({bEvents.event}, 'BLEACHCOUNT'));
[LONG,LAT] = m_ll2xy([events.lon],[events.lat]); hold on % convert reef points to M-Map lat long
%scatter(Mort_stats(:,1),Mort_stats(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
scat14a = scatter(LONG,LAT,5, log2([events.count])) ; % plot bleaching events onto map
colormap(jet)
colorbar
title('Bleaching Events Between 1861-2100 (log base 2)')
print('-dpdf', '-r200', strcat(fullDir, filePrefix,'_AllMortEventsMap_log2','.pdf'));
hold off;


% % % Make map of prop constants 
% load ('~/Dropbox/Matlab/SymbiontGenetics/mat_files/Updated_psw2.mat','psw2_new')
% figure(15)
% m_proj('miller');
% m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
% m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
% [LONG,LAT] = m_ll2xy(Reefs_latlon(:,1),Reefs_latlon(:,2)); hold on % convert reef points to M-Map lat long
% scatter(Reefs_latlon(:,1),Reefs_latlon(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
% scatter(LONG,LAT,5, psw2_new(:,1)) ; % plot bleaching events onto map
% colorbar
% title('Proportionality Constant')
% print('-dpdf',strcat(dateString,'_','PropConstantMap','.pdf')); 
% 
% % Make map of var(SSThist)
% load ('~/Dropbox/Matlab/SymbiontGenetics/mat_files/Updated_psw2.mat','psw2_new')
% figure(16)
% m_proj('miller');
% m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
% m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
% [LONG,LAT] = m_ll2xy(Reefs_latlon(:,1),Reefs_latlon(:,2)); hold on % convert reef points to M-Map lat long
% scatter(Reefs_latlon(:,1),Reefs_latlon(:,2), 5, [.7 .7 .7]) % plot reefs  onto map
% scatter(LONG,LAT,5, psw2_new(:,2)) ; % plot bleaching events onto map
% colorbar
% title('SST Variance')
% print('-dpdf',strcat(dateString,'_','SSTvarMap','.pdf')); 
% 
end