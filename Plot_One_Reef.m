%% PLOTS FOR TEMP, SYMBIONT DENSITY, SYMBIONT GENOTYPE, CORAL COVER
% Inputs:
% C    - Coral population, downsampled to monthly in advance.
% S    - Symbiont density, downsampled to monthly in advance.
% bleachEvent - sparse boolean array with true for each bleaching event.
%               1 column per coral type.
% psw2 - selectional variance
% time - time axis for plotting
% temp - SST history
% lat, lon - latitude and longitude of current reef area
% RCP  - the current Representative Concentration Pathway name
% hist - historical mean temperature
% C    - Coral population history
% S    - Symbiont population history
% Data - 1 for ESM2M, 2 for HadISST dataset
% I    - index of a reference date in the time array
% gi   - symbiont genotype
% ri   - symbiont growth rate
% SGPath - where to find the DHM mat file
% outputPath  - Directory for m-files and output subdirectories
% k    - number of the current reef grid cell
% pdfDirectory  - Output directory for the current run
% LOC  - location of current cell, like Loc but with underscores
% NF   - normalization factor, currently 1
% E    - evolution on or off
% dateString - today's date as formatted in the calling program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% 12-15-15                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Plot_One_Reef(C, S, bleachEvent, psw2, time, temp, lat, lon, RCP, ...
            hist, Data, SGPath, outputPath, k, pdfDirectory, LOC, ...
            NF, E, dateString, months)
    % Note that the persistent names are the same as in Plot_SST_Decimate.
    % The namespaces should be separate do there's no side effect.
    persistent figHandle DHM Time; %#ok<PSET>
    persistent plotDHM1; % doesn't need to be modified: plotDHM2;
    persistent plotSST1 plotSST2;
    persistent plotSD1 plotSD2;
    persistent plotCC1 scatCC2 scatCC1; 
    persistent timeY;
    if nargin == 0
        % Clear variables for a clean start
        %disp('Clean start in Plot worker');
        close all;
        clearvars figHandle DHM Time plotDHM1 plotDHM2 plotSST1 plotSST2 plotSD1 plotSD2 plotCC1 scatCC2 scatCC1 timeY; 
        return;
    end
    
    fontSize = 18;
    % control the X axis ticks;
    fiftyYears = 365.25*50;  % 1 unit per day in datetime form
    tickStart = floor(1850*365.25);
    tickEnd = floor(2100*365.25);
    tickVals = tickStart:fiftyYears:tickEnd;
    %tickVals = 1850:50:2100;
    
    smoothSymbiont = 6; % number of points +- for moving average. 0 = no smoothing
    
    shrink = true;
    % Shrink all time series to monthly for faster plotting.  Note that
    % since this function does not return the variables, it is safe to
    % modify them in place.
    if shrink
        factor = round(length(time)/months);
        time = decimate(time, factor, 'fir');
        temp = decimate(temp, factor, 'fir');
    end
    
    % Needed every time for output:
    shortprop = sprintf('%.2f', psw2); 
    
    % Set up the figure on the first call, or (for matlab >=8.4) assign
    % new data to the existing plots.
    if isempty(figHandle) || (~verLessThan('matlab', '8.4') && ~isvalid(figHandle))
        %disp('New figure handle for figs');
        %close all; 
        figHandle = figure(3001);
        set(figHandle, 'color', 'w');
        %set(gcf,'Visible', 'off');
        % Because the figures are saved in "fig" format they need to be
        % visible or they won't open easily (openfig will do it, though).
        set(figHandle,'Visible', 'on');
        % Make sure this is the same "headless" or not:
        set(groot,'defaultFigurePaperPositionMode','manual');

        % Make serial and parallel calls match:
        % PP is left, bottom, width, height
        % May not be doing anything!
        % original:
        %set(gcf, 'PaperPosition', [1, 1.5, 6.5, 8]);

        set(gcf, 'PaperSize', [8,8], 'OuterPosition', [100, 100, 1200, 1200], 'PaperPosition', [0.5, 0.5, 12, 5], 'PaperUnits', 'inches');
        first = true;
    else
        % Below fixes figure overwrite, but makes figure visible.  Second
        % line only takes effect after flashing to the screen.
        figure(figHandle);  % use existing
        
        first = false;

        set(plotSST1, 'YData', temp);
        set(plotSST2, 'YData', [hist hist]);
        % Only SST has a variable title.
        subplot(2,2,1); % 1 row, 3 columns
        Loc = strcat(num2str(lat),',',num2str(lon));
        %title(strcat('SST ESM2M',RCP ,';prop ', num2str(shortprop),';  latlon:', Loc));
        title(strcat(RCP ,'  prop ', num2str(shortprop),'  latlon:', Loc));

        if smoothSymbiont > 0
            smooth1 = movmean(S(:,1)./C(:,1), [smoothSymbiont smoothSymbiont]);
            smooth2 = movmean(S(:,2)./C(:,2), [smoothSymbiont smoothSymbiont]);
            set(plotSD1, 'YData', smooth1);
            set(plotSD2, 'YData', smooth2);
        else
            set(plotSD1, 'YData', S(:,1)./C(:,1));
            set(plotSD2, 'YData', S(:,2)./C(:,2));
        end
        set(plotSD1, 'YData', smooth1);
        set(plotSD2, 'YData', smooth2);


        set(plotCC1, {'YData'}, {C(:,1), C(:,2)}');
        
        % For the bleaching markers, the X values change each time, since
        % there is just one even, so we set YData too.
        [subscripts, ~] = find(bleachEvent(:,2));
        tSub = timeY(subscripts);
        subscripts(subscripts > 0) = 1000000;  % Arbitrary value above y=0
        set(scatCC2, 'XData', tSub);
        set(scatCC2, 'YData', subscripts);
        [subscripts, ~] = find(bleachEvent(:,1));
        tSub = timeY(subscripts);
        subscripts(subscripts > 0) = 1000000;  % Arbitrary value above y=0
        set(scatCC1, 'XData', tSub); 
        set(scatCC1, 'YData', subscripts); 

        set(plotDHM1, 'YData', DHM(k,:));
    end
    
    % Awkward logic since this code could have been inside the structure
    % above, but that's how it is for now.
    if first
        subplot(2,2,1); % Three rows and two columns of plots

        %% Plot Temperature (upper left)
        plotSST1 = plot(time,temp,'k');  %(1069:1693) gives 1960-2011
        hold on;
        Loc = strcat(num2str(lat),',',num2str(lon));
        if Data == 1 % ESM2M
            %     if OA ==1
            %         plot(time,Omega,'b');  %(1069:1693) gives 1960-2011
            %         legend('ESM2MSST','Omega' ,'Location','NorthWest');
            %         ylim([0 32]);
            %         title(strcat('SST and Omega at',RCP ,'- ', LOC)) %
            %         ylabel('SST (C) and Omega')
            %     end
            
            % no legend - moved prop to title
            %legend(strcat('ESM2M', RCP, ';prop ',num2str(shortprop) ),'Location','NorthWest');
            hold on;
            %title(strcat('SST ESM2M',RCP ,';  latlon:', Loc)); %
            %title(strcat('SST ESM2M',RCP ,';prop ', num2str(shortprop),';  latlon:', Loc));
            title(strcat(RCP ,'  prop ', num2str(shortprop),  '  latlon:', Loc));
            
            ylabel('SST (C) ');
        elseif Data == 2 % HadISST
            legend('HadISST','Location','NorthWest');
            title(strcat('SST',';prop ', num2str(shortprop),';  latlon:', Loc)); %
            ylabel('SST (C)');
            
        end
        plotSST2 = plot([time(1) time(end)],[hist hist],'r');

        xlabel('Time (years)');
        xlim([time(1) time(end)]);
        ylim([17 35]);
        set(gca, 'FontSize',fontSize);
        %set(gca, 'XTick', tickVals);
        datetick('x','keeplimits');


        %% Plot All Symbionts Density (lower left)
        subplot(2,2,3);
%    pos = get(gca, 'Position');
%    pos(1) = 0.055;
%    pos(3) = 0.9;
%    set(gca, 'Position', pos)
        % smooth the data so plots can be overlaid
        % for the original plots, just use the S/C argument without
        % smoothing.
        if smoothSymbiont > 0
            smooth1 = movmean(S(:,1)./C(:,1), [smoothSymbiont smoothSymbiont]);
            smooth2 = movmean(S(:,2)./C(:,2), [smoothSymbiont smoothSymbiont]);
            plotSD1 = plot(time, smooth1,'color','y','LineWidth',1);
            hold on;
            plotSD2 = plot(time, smooth2,'color','g','LineWidth',1);
        else
            plotSD1 = plot(time, S(:,1)./C(:,1),'color','y','LineWidth',1);
            hold on;
            plotSD2 = plot(time, S(:,2)./C(:,2),'color','g','LineWidth',1);
        end

        xlabel('Time (years)') ;
        ylabel('Mean Symbiont Density (cells/cm^2)');
        title('Symbiont Densities') ;
        legend('massive', 'branching','Location','NorthWest');
        %legend('massive', 'branching',strcat('+',num2str(x), 'C (m)'), strcat('+',num2str(x),'C (b)'),'Location','NorthWest');
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);3
        ylim([0 5e+06]);
        %YTick(linspace(0, 5e+06, 6));
       
        set(gca, 'FontSize',fontSize);
        datetick('x','keeplimits');
        %set(gca,'XTick',tickVals);

        %print -dpdf -r600 fig2.pdf


        %% Plot Coral Cover (upper right)
        subplot(2,2,2);

        plotCC1 = plot(time,C(:,1),'m',time,C(:,2),'b');
        % figure
        % plot(time, M)
        hold on;
        %scatter(time(locs+12), M(locs+12),'r','filled')
        %datestr(time(6719)) % Dec 2010
        %CoralCover_Branch_2010(k,1) = C(6718/dt*.25+1,2);
        
        % For the scatter data, we're just drawing a circle on the x axis.
        % approximating to min-year will be close enough.
        timeY(1:length(bleachEvent)) = time(6:12:length(time));
        [subscripts, ~] = find(bleachEvent(:,2));
        tSub = timeY(subscripts);
        subscripts(subscripts > 0) = 1000000;  % Arbitrary value above y=0
        scatCC2 = scatter(tSub,subscripts,'b'); datetick % branching mort event
       
        [subscripts, ~] = find(bleachEvent(:,1));
        tSub = timeY(subscripts);
        subscripts(subscripts > 0) = 1000000;
        scatCC1 = scatter(tSub, subscripts,'m');  % massive mort event
       
        xlabel('Time (years)')   ; 

        ylabel('Coral Cover (cm^2)');
        ylim([0 1.0e+08]);
        title('Coral Population Size') ;
        legend('massive', 'branching','Location','best');

        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        set(gca, 'FontSize',fontSize);

        set(gca,'XTick',tickVals);

        datetick('x','keeplimits');

        %legend('massive','branching','Location','NorthWest');
        %print -dpdf -r600 fig4.pdf

       


        %% Plot DHMs from RCP85 noadapt MMMmax method (lower right)
        subplot(2,2,4);
        if isempty(DHM)
            % Select the matching file:
            rcpDigits = RCP(end-1:end);
            % for control cases it won't be a positive integer.
            rcpNum = str2double(rcpDigits);
            if rcpNum > 0
                 % load DHMs from RCP(rcpNum) noadapt MMMmax method
                load(strcat(SGPath, 'DHM_', rcpDigits, 'noadaptMMMmax.mat'),'DHM','Time'); 
            else
                % Default to RCP85 - better than nothing???
                load(strcat(SGPath, 'DHM_85noadaptMMMmax.mat'),'DHM','Time');
            end
        end
        plotDHM1 = plot(Time(1501:end),DHM(k,:), 'k'); hold on;
        plotDHM2 = plot([time(1) time(end)],[2 2],'r'); %#ok<NASGU>
        datetick;
        legend('DHMs','DHM=2','Location','NorthWest');
        title('Degree Heating Months (MMMmax)');
        xlabel('Time (years)');
        ylabel('DHM (no adapt)');
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        ylim([0 10]);
        set(gca, 'FontSize',fontSize);
        set(gca,'XTick',tickVals);
        datetick('x','keeplimits');

    end
    
    %% Save with Date Stamp
    format shortg;
    % The print line takes over half the time of the whole routine.  Break it
    % down and see why.
    if Data == 1
        name = strcat('SDC_',dateString,'_',num2str(k),'_normSST',RCP,LOC,'prop',num2str(shortprop),'_NF',num2str(NF),'_E',num2str(E));
    else
        name = strcat('SDC_',dateString,'_',num2str(k),'_HADISST',LOC,'prop',num2str(shortprop),'_NF',num2str(NF),'_E',num2str(E));
    end
    fullName = strcat(outputPath, pdfDirectory, name);
    if verLessThan('matlab', '8.2')
        saveas(gcf, name, 'fig');
    else
        savefig(strcat(fullName, '.fig'));
    end
    
end