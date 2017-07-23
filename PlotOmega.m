% Plot Omega
% time from 0 to 2880:
tStep = 2800;
showFactor = 1;
figure();

m_coast('patch',[0.7 0.7 0.7],'edgecolor','none');
hold on;
m_proj('miller');
m_grid('box','fancy','linestyle','none','backcolor',[.9 .99 1], 'xticklabels', [], 'yticklabels', []);
[LONG2,LAT2] = m_ll2xy(Reefs_latlon(:,1), Reefs_latlon(:,2));
om = Omega_all(:, tStep);
if showFactor == 1
    om = (1-(4-Omega_all(:, tStep))*0.15);
end
scatter(LONG2, LAT2, 5, om)
colorbar
if showFactor
    % Smallest factor possible is 0.55 (or zero if omega < 1), but
    % 0.7 seems to be close to a lower bound for rcp 6.0.
    caxis([.7 1]);
else
    caxis([2 4]);
end
