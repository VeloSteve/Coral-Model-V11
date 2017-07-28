function BHSCO()
%relPath = '../bleaching_history/';
relPath = 'D:/CoralTest/V11Test_SC/bleaching/';

figure('color', 'w');
set(gcf, 'Units', 'inches', 'Position', [1, 1.5, 17, 11]);

% Make 4 subplots across, each for one rcp scenario.
% Within each, have lines for E=0/1 and OA=0/1
rcpList = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
rcpName = {'RCP 2.6', 'RCP 4.5', 'RCP 6.0', 'RCP 8.5'};

legText = {};
legCount = 1;
for i = 1:4
    rrr = rcpList{i};
    subplot(2,2,i);
    cases = [];
    hFile = '';
    for eee = 0:1
        for ooo = 1:-1:0           
            if ~(eee == 0 && ooo == 1)  % Skip one curve
                hFile = strcat(relPath, 'BleachingHistory', rrr, 'E=', num2str(eee), 'OA=', num2str(ooo), '.mat');
                load(hFile, 'yForPlot');
                cases = [cases; yForPlot];
                legText{legCount} = strcat('E = ', num2str(eee), ' OA = ', num2str(ooo));
                legCount = legCount + 1;
            end
        end
    end
    % x values are the same for all. Use the latest file;
    load(hFile, 'xForPlot');
    oneSubplot(xForPlot, cases, legText, rcpName{i}, i);
    if i == 1
        yHandle = ylabel({'Percent Current Bleaching,'; 'Frequent Bleaching, or Mortality'},'FontSize',22);
        set(yHandle, 'Position', [1925 -10 0])
    end

end
end

function oneSubplot(X, Yset, legText, tText, baseColor) 

    % Create axes
    %axes1 = axes;
    %hold(axes1,'on');
    switch baseColor
        case 1 
            base = [0 0 .9];
            light = [0.5 0.5 1.0];
        case 2 
            base = [0 .7 0];
            light = [0.4 1.0 0.4];
        case 3 
            base = [.6 .6 0];
            light = [0.9 0.9 0];
        case 4 
            base = [.9 0 0];
            light = [1.0 0.5 0.5];
    end
    col{1} = base;    % XXX - NOTE that this line is hidden for now.
    col{2} = base;

    col{3} = light;
    col{4} = light;
    
    % Create multiple lines using matrix input to plot
    plot1 = plot(X,Yset(:, :));
    for i = size(Yset, 1):-1:1
        set(plot1(i),...
            'DisplayName',legText{i}, ...
            'Color', col{i}, ...
            'LineWidth', 2);
        if i == 2 || i == 4
            set(plot1(i), 'LineStyle', '--');
        end     
    end

    % Create xlabel
    xlabel('Year','FontSize',13.2);

    % Create title
    title(tText, 'FontSize', 22);


    xlim([1950 2100]);
    ylim([0 105]);
    box('on');
    grid('on');
    % Set the remaining axes properties
    set(gca, 'FontSize',22,'XTick',[ 1950 2000 2050 2100 ],...
        'YTick',[0  50  100]);
    % Create legend
    legend1 = legend('show');
    set(legend1,'Location','northwest','FontSize',20);


    
end