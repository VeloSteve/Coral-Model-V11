%% PLOTS FOR TEMP, SYMBIONT DENSITY, SYMBIONT GENOTYPE, CORAL COVER
% Inputs:
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
% dCov - mortality events
% ri   - symbiont growth rate
% basePath - where to find the DHM mat file
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
function Plot_SST_Decimate(psw2, time, temp, lat, lon, RCP, ...
            hist, C, S, Data, I, gi, dCov, ri, basePath, outputPath, k, pdfDirectory, LOC, ...
            NF, E, dateString, months)
    persistent figHandle DHM Time;
    persistent plotSG1 plotSG2 ;
    persistent plotSP1 plotSP2;
    persistent plotDHM1 plotDHM2;
    persistent plotSST1 plotSST2;
    persistent plotSD1 plotSD2;
    persistent plotCC1 scatCC2 scatCC1; 
    if nargin == 0
        % Clear variables for a clean start
        %disp('Clean start in Plot worker');
        close all;
        clearvars figHandle DHM Time plotSG1 plotSG2  plotSP1 plotSP2 plotDHM1 plotDHM2 plotSST1 plotSST2 plotSD1 plotSD2 plotCC1 scatCC2 scatCC1; 
        return;
    end
    
    shrink = true;
    % Shrink all time series to monthly for faster plotting.  Note that
    % since this function does not return the variables, it is safe to
    % modify them in place.
    if shrink
        factor = round(length(time)/months);
        time = decimate(time, factor, 'fir');
        temp = decimate(temp, factor, 'fir');
        temp1 = decimate(C(:, 1), factor, 'fir');
        temp2 = decimate(C(:, 2), factor, 'fir');
        C = [temp1 temp2];
        temp1 = decimate(S(:, 1), factor, 'fir');
        temp2 = decimate(S(:, 2), factor, 'fir');
        S = [temp1 temp2];
        temp1 = decimate(gi(:, 1), factor, 'fir');
        temp2 = decimate(gi(:, 2), factor, 'fir');
        gi = [temp1 temp2];
        % dCov is an array of logicals, not continuous values, so it needs
        % special treatment.
        % Group each part of dCov into "factor" rows and check for a 1 in
        % any row.  This yields a logical, and we must then get back to
        % double for plotting.
        dCov1 = reshape(dCov(:, 1), factor, []);
        dCov2 = reshape(dCov(:, 2), factor, []);
        dCov = double([any(dCov1(:, :))' any(dCov2(:, :))']);
        dCov(dCov == 0) = NaN ; % make zero's NaN so they don't plot
        % old:
        %{
        temp1 = decimate(dCov(:, 1), factor, 'fir');
        temp2 = decimate(dCov(:, 2), factor, 'fir');
        dCov = [temp1 temp2];
        %}
        temp1 = decimate(ri(:, 1), factor, 'fir');
        temp2 = decimate(ri(:, 2), factor, 'fir');
        ri = [temp1 temp2];
        I = I/4;  % XXXX should be I*factor?
    end
    
    % Needed every time for output:
    shortprop = sprintf('%.2f', psw2); 
    
    % Plots
    %figure(Data)
    % gcf refers to or creates the current figure handle
    % close all ensures that older visible plots don't pop up.  This is
    % also very important for preventing memory leak when using parallel
    % workers to do the pdf conversion.
    % The verLessThan check is meant to make this compatible with the old
    % and new graphics systems which changed at version 8.4 (2014b).  Needs
    % testing.
    if isempty(figHandle) || (~verLessThan('matlab', '8.4') && ~isvalid(figHandle))
        %disp('New figure handle for figs');
        %close all; 
        figHandle = figure(3001);
        %set(gcf,'Visible', 'off');
        set(figHandle,'Visible', 'off');
        % Make sure this is the same "headless" or not:
        set(groot,'defaultFigurePaperPositionMode','manual');

        % Summary of tests attempting to get consistent results in parallel
        % and serial.
        % Settings                      Parallel file size (kb)     Serial file size (kb)
        % minimal, pdf                  1010                        103
        % -bestfit (nicer aspect ratio) 1010                        96  
        % -r200  (only affects serial)  1010                        141
        % -painters (meant for vector)  999                         1002
        % -opengl (meant for raster)    1087                        141
        % -deps (requires -bestfit off, parallel doubles, serial same, harder to work with).
        % -djpeg - dpainters -r200      218                         223 <-- best yet!
        % Tradeoffs:
        % With pdf and -rXXX can improve serial quality, but not cut back parallel.
        % With jpeg can control both but it takes about 1/3 longer than pdf.
        
        %set(gcf, 'FigurePosition', [200 200 600 500]);
        % Make serial and parallel calls match:
        set(gcf, 'PaperPosition', [1, 1.5, 6.5, 8]);
        first = true;
    else
        %fprintf('Old figure handle %d for figs\n', [figHandle.Number]);
        %set(groot, 'CurrentFigure', figHandle);
        % Below fixes figure overwrite, but makes figure visible.  Second
        % line only takes effect after flashing to the screen.
        figure(figHandle);  % use existing
        %set(gcf,'Visible', 'off'); % but "figure" toggles it to visible.
        
         first = false;

        set(plotSST1, 'YData', temp);
        set(plotSST2, 'YData', [hist hist]);
        % Only SST has a variable title.
        subplot(3,2,1);
        Loc = strcat(num2str(lat),',',num2str(lon));
        title(strcat('SST ESM2M',RCP ,';prop ', num2str(shortprop),';  latlon:', Loc));


        set(plotSD1, 'YData', S(:,1)./C(:,1));
        set(plotSD2, 'YData', S(:,2)./C(:,2));


        set(plotCC1, {'YData'}, {C(:,1), C(:,2)}');
        set(scatCC2, 'YData', dCov(:,2));
        set(scatCC1, 'YData', dCov(:,1)); 
        set(plotSP1, 'YData', ri(:,1));
        set(plotSP2, 'YData', ri(:,2));

        set(plotSG1, 'YData', gi(:,1));
        set(plotSG2, 'YData', gi(:,2));
        set(plotDHM1, 'YData', DHM(k,:));
    end
    


    % Awkward logic since this code could have been inside the structure
    % above, but that's how it is for now.
    if first

        subplot(3,2,1); % Three rows and two columns of plots

        %% Plot Temperature & Omega (upper left)
        % figure(5) % 1
        plotSST1 = plot(time,temp,'k');  %(1069:1693) gives 1960-2011
        hold on;
        Loc = strcat(num2str(lat),',',num2str(lon));
        if Data == 1; % ESM2M
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
            title(strcat('SST ESM2M',RCP ,';prop ', num2str(shortprop),';  latlon:', Loc));

            ylabel('SST (C) ');
        elseif Data == 2; % HadISST
            legend('HadISST','Location','NorthWest');
            title(strcat('SST',';prop ', num2str(shortprop),';  latlon:', Loc)); %
            ylabel('SST (C)');
        end
        % "end" below was "I"
        plotSST2 = plot([time(1) time(end)],[hist hist],'r');
        xlabel('Time (years)');
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        ylim([20 32]);
        datetick('x','keeplimits');
        %print -dpdf -r600 fig1.pdf


        %% Plot All Symbionts Density (lower right)
        %figure(2)

         subplot(3,2,6);
        col=colormap;
        plotSD1 = plot(time,S(:,1)./C(:,1),'color','y','LineWidth',1);
        hold on;
        plotSD2 = plot(time,S(:,2)./C(:,2),'color','g','LineWidth',1);
        % hold on; plot(time,S(:,3)./C(:,1),'color',col(45,:),'LineWidth',1)
        % hold on; plot(time,S(:,4)./C(:,2),'color',col(25,:),'LineWidth',1)
        % hold on; plot(time,S(:,5)./C(:,1),'color','r','LineWidth',1)
        % hold on; plot(time,S(:,6)./C(:,2),'color',col(18,:),'LineWidth',1)
        % hold on; plot(time,S(:,7)./C(:,1),'color','k','LineWidth',1)
        % hold on; plot(time,S(:,8)./C(:,2),'color',col(5,:),'LineWidth',1)

        xlabel('Time (years)') ;
        ylabel('Mean Symbiont Density (cells/cm^2)');
        title('Tolerant Symbiont Densities') ;
        legend('massive', 'branching','Location','NorthWest');
        %legend('massive', 'branching',strcat('+',num2str(x), 'C (m)'), strcat('+',num2str(x),'C (b)'),'Location','NorthWest');
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);3
        ylim([0 6e+06]);
        datetick('x','keeplimits');
        %print -dpdf -r600 fig2.pdf

        %% Plot Symbiont Genotype (middle right)
        %figure(3)
        subplot(3,2,4);

        plotSG1 = plot(time,gi(:,1),'color','y','LineWidth',2);
        hold on;
        plotSG2 = plot(time,gi(:,2),'color','g','LineWidth',2);
        % hold on; plot(time,gi(:,3),'color',col(45,:),'LineWidth',2)
        % hold on; plot(time,gi(:,4),'color',col(25,:),'LineWidth',2)
        % hold on; plot(time,gi(:,5),'r','LineWidth',2)
        % hold on; plot(time,gi(:,6),'color',col(18,:),'LineWidth',2)
        % hold on; plot(time,gi(:,7),'k','LineWidth',2)
        % hold on; plot(time,gi(:,8),'color',col(5,:),'LineWidth',2)

        ylim([24 28]);
        %ylim([ mean(gi(:,1))-1  mean(gi(:,2))+1 ])
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        xlabel('Time (years)');
        ylabel('Mean Symbiont Genotype (C)');
        title('Symbiont Genotypes');
        Lh=legend('massive','branching','Location','NorthWest');
        %Lh=legend('massive','branching',strcat('+',num2str(x), 'C (m) '),strcat('+',num2str(x), 'C (b)'),'Location','NorthWest');
        set(Lh,'color','none');
        datetick('x','keeplimits');
        %print -dpdf -r600 fig3.pdf
        xlim([time(1) time(end)]);

        %% Plot Coral Cover (middle left)
        %figure(4)
        subplot(3,2,3) ;

        plotCC1 = plot(time,C(:,1),'m',time,C(:,2),'b');
        % figure
        % plot(time, M)
        hold on;
        %scatter(time(locs+12), M(locs+12),'r','filled')
        %datestr(time(6719)) % Dec 2010
        %CoralCover_Branch_2010(k,1) = C(6718/dt*.25+1,2);
        scatCC2 = scatter(time,dCov(:,2),'b'); datetick % branching mort event
        scatCC1 = scatter(time,dCov(:,1),'m');  % massive mort event

        xlabel('Time (years)')   ; 

        ylabel('Coral Cover (cm^2)');
        ylim([0 1.0e+08]);
        title('Coral Population Size') ;
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        datetick('x','keeplimits');

        %legend('massive','branching','Location','NorthWest');
        %print -dpdf -r600 fig4.pdf

        %% Plot symbiont population growth rate (upper right)
        subplot(3,2,2);
        plotSP1 = plot(time,ri(:,1),'color','y','LineWidth',2);
        hold on; plotSP2 = plot(time,ri(:,2),'color','g','LineWidth',2); datetick;
        title('Symbiont Population Growth Rate');
        xlabel('Time (years)');
        ylabel('Symbiont Population Growth Rate');
        xlim([time(1) time(end)]);


        %% Plot DHMs from RCP85 noadapt MMMmax method (lower left)
        subplot(3,2,5);
        if isempty(DHM)
            load('DHM_85noadaptMMMmax.mat','DHM','Time'); % load DHMs from RCP85 noadapt MMMmax method
        end
        plotDHM1 = plot(Time(1501:end),DHM(k,:)); hold on;
        plotDHM2 = plot([time(1) time(end)],[2 2],'r');
        datetick;
        legend('DHMs','DHM=2','Location','NorthWest');
        title('Degree Heating Months (MMMmax)');
        xlabel('Time (years)');
        ylabel('DHM (no adapt)');
        xlim([time(1) time(end)]);
        %xlim([time(1) time(960)]);
        ylim([0 10]);
        datetick('x','keeplimits');

    end
    
    %% Save with Date Stamp
    format shortg;
    % The print line takes over half the time of the whole routine.  Break it
    % down and see why.
    if Data == 1;
        name = strcat(dateString,'_',num2str(k),'_normSST',RCP,LOC,'JD_B_hist_prop',num2str(shortprop),'_NF',num2str(NF),'_E',num2str(E),'min0v2.pdf');
    else
        name = strcat(dateString,'_',num2str(k),'_HADISST',LOC,'JD_B_hist_prop',num2str(shortprop),'_NF',num2str(NF),'_E',num2str(E),'min0v2.pdf');
    end
    fullName = strcat(outputPath, pdfDirectory, name);
    if ~verLessThan('matlab', '8.4')
        print('-dpdf', '-r200', fullName);  % Older versions don't have "-bestfit".
    else
        print('-dpdf', '-bestfit', '-r200',  fullName);
    end

    % print('-djpeg', '-painters', '-r200',  name);
    % Very slow, but gives high-quality vector output: export_fig(name, '-pdf', handle);
    
    
    % change x limits for a subset of years
    %{  
    May be used in the future
    subplot(3,2,1);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    subplot(3,2,2);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    subplot(3,2,3);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    subplot(3,2,4);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    subplot(3,2,5);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    subplot(3,2,6);
    xlim([time(4850) time(end)]); % Jan 1982-1986
    %}
end