% Plot old C versus Prince-Dormand C.

figure()
plot(time, C(:,1), 'm')
hold on;
plot(time, C(:,2), 'b')

% Time for PD is in months, but for old data it's in equally spaced double
% time values.  Convert, noting that both start at the same time.
start = time(1);
for i = length(tPD):-1:1
    tConv(i) = addtodate(start, round(tPD(i)), 'month');
end
datenum(strcat(num2str(startYear),'-01-01'));

plot(tConv, SPD(:,1), 'k')
plot(tConv, SPD(:,2), 'k')