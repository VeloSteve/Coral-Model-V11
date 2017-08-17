% Plot old C versus Dormand-Prince C.

figure()
plot(time, C(:,1), 'DisplayName', 'Massive', 'Color', [1 0 1])
hold on;
plot(time, C(:,2), 'DisplayName', 'Branching', 'Color', [0 0 1])

% Time for PD is in months, but for old data it's in equally spaced double
% time values.  Convert, noting that both start at the same time.
start = time(1);
for i = length(tPD):-1:1
    tConv(i) = addtodate(start, round(tPD(i)), 'month');
end
datenum(strcat(num2str(startYear),'-01-01'));

plot(tConv, CPD(:,1), 'DisplayName', 'D-P', 'Color', [0 0 0])
plot(tConv, CPD(:,2), 'DisplayName', 'D-P', 'Color', [0 0 0])
hold off;
datetick('x', 'keeplimits')
legend('show');
title('Corals - compare original R-K to Dormand-Prince, reef 366');

% Now the symbionts
figure()
plot(time, S(:,1), 'DisplayName', 'Massive', 'Color', [1 1 0])
hold on;
plot(time, S(:,2), 'DisplayName', 'Branching', 'Color', [0 1 0])

% Time for PD is in months, but for old data it's in equally spaced double
% time values.  Convert, noting that both start at the same time.
start = time(1);
for i = length(tPD):-1:1
    tConv(i) = addtodate(start, round(tPD(i)), 'month');
end
datenum(strcat(num2str(startYear),'-01-01'));

plot(tConv, SPD(:,1), 'DisplayName', 'D-P', 'Color', [0 0 0])
plot(tConv, SPD(:,2), 'DisplayName', 'D-P', 'Color', [0 0 0])
hold off;
datetick('x', 'keeplimits')
legend('show');
title('Symbionts - compare original R-K to Dormand-Prince, reef 366');
