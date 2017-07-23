%% PROPORTIONALITY CONSTANT CALCULATIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% last updated: 8/26/16                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[basePath, outputPath, sstPath, matPath, Computer] = useComputer(3);


% SST DATASET?
Data = 1; % 1=ESM2M_norm
%Data = 2; % 2=HADISST (through 3/16/16)
if Data == 1, dataset = 'ESM2M'; else dataset = 'HadISST'; end

%% DEFINE CLIMATE CHANGE SCENARIO (from normalized GFDL-ESM2M; J Dunne)
%for i={'rcp26','rcp85'}  
%RCP = char(i)
RCP = 'rcp85'; % options; 'rcp26', 'rcp45', 'rcp60', 'rcp85'
format shortg; c = clock; date = strcat(num2str(c(1)),num2str(c(2)),num2str(c(3))); % today's date stamp

%% LOAD JOHN'S NORMALIZED SSTS FROM EARTH SYSTEM CLIMATE MODEL OR HADISST
%GetSST_norm_GFDL_ESM2M % sub m-file for extracting SSTs for a ALL reef grid cells
[SST, Reefs_latlon, TIME, startYear] = GetSST_norm_GFDL_ESM2M(sstPath, matPath, Data, RCP);

SST_1861_2000 = SST(:,1:1680);
SSThist = SST_1861_2000;

%% LOAD OLD SELECTIONAL VARIANCE (psw2) 
%load ('~/Dropbox/Matlab/SymbiontGenetics/mat_files/psw2_trials.mat','psw2_var_allv2')

%% Store optimizer inputs from propInputValues with any constant values to be computed.
pswInputs(:,1) = propInputValues';
pswInputs(:,2) = [0.35; 2.0; 0.5; 2.517];
pswInputs(:,3) = [0.35; 2.0; 0.5; 3.104];
pswInputs(:,4) = [0.35; 2.0; 0.5; 2.390];
pswInputs(:,5) = [0.35; 2.0; 0.5; 3.063];
% Newer sets with case 5 constants and target 2% bleaching
pswInputs(:,6) = [0.35; 1.8; 0.5; 2.9925];  % RCP 2.6, E=0
pswInputs(:,7) = [0.35; 1.8; 0.5; 3.1475];  % RCP 8.5, E=0
pswInputs(:,8) = [0.35; 1.8; 0.5; 3.7800];  % RCP 2.6, E=1
pswInputs(:,9) = [0.35; 1.8; 0.5; 3.8000];  % RCP 8.5, E=1
% 1/25/2017 sets with case 5 constants and target 5% bleaching
pswInputs(:,10) = [0.35; 1.5; 0.26; 2.0833];  % RCP 2.6, E=0
pswInputs(:,11) = [0.35; 1.5; 0.26; 2.1111];  % RCP 8.5, E=0
pswInputs(:,12) = [0.35; 1.5; 0.46; 4.4500];  % RCP 2.6, E=1
pswInputs(:,13) = [0.35; 1.5; 0.46; 4.4667];  % RCP 8.5, E=1
% 1/26/2017 sets with case 5 constants and target 10% bleaching
% For 10% it seems impossible to keep any of the input values constant.
pswInputs(:,14) = [0.3667; 1.3; 0.375; 4.0000];  % RCP 2.6, E=0
pswInputs(:,15) = [0.4222; 1.2222; 0.4556; 5.1111];  % RCP 8.5, E=0
pswInputs(:,16) = [0.3556; 1.3; 0.3417; 12.2223];  % RCP 2.6, E=1
pswInputs(:,17) = [0.6; 1.3; 0.3777; 4.5556];  % RCP 8.5, E=1

% 1/27/2017 alternate values for E=1 at 10% bleaching.  This produces
% less T-shaped variance parameter distributions than those just above, 
% and better parameter matches to Baskett 2009, but not more exact matches
% to the 10.0% bleaching target.
pswInputs(:,18) = [0.36; 1.2000; 0.4167; 5.0000];  % RCP 2.6, E=1
pswInputs(:,19) = [0.36; 1.2222; 0.4778; 6.3333];  % RCP 8.5, E=1
% 2/4/2017 sets with case 5 constants and target 5% bleaching
% This has the fix to make ri variable in Runge-Kutta.
pswInputs(:,20) = [0.36; 1.5; 0.46; 4.1850];  % RCP 2.6, E=0
pswInputs(:,21) = [0.36; 1.5; 0.46; 4.2550];  % RCP 8.5, E=0
pswInputs(:,22) = [0.36; 1.5; 0.46; 4.7100];  % RCP 2.6, E=1
pswInputs(:,23) = [0.36; 1.5; 0.46; 4.7156];  % RCP 8.5, E=1
% 2/28/17 Add rcp4.5 and 6.0
pswInputs(:,24) = [0.36; 1.5; 0.46; 4.3617];  % RCP 4.5, E=0
pswInputs(:,25) = [0.36; 1.5; 0.46; 4.2077];  % RCP 6.0, E=0
pswInputs(:,26) = [0.36; 1.5; 0.46; 4.8444];  % RCP 4.5, E=1
pswInputs(:,27) = [0.36; 1.5; 0.46; 4.7156];  % RCP 6.0, E=1

[~, pswCount] = size(pswInputs);


%% CALULATE PROP CONSTANT FOR EACH GRID CELL
psw2_new = nan(length(Reefs_latlon),1); % initialize matrix
for reef = 1:length(Reefs_latlon)
    SSThistReef = SST_1861_2000(reef,:)';  % extract SSTs for grid cell bn 1861-2000
    % XXX Note that SSThistReef was created but not used in the code I
    % received!
    %psw2(k,1) = max(0.7,min(1.3,(mean(exp(0.063*SSThistReef))/var(exp(0.063*SSThistReef)))^0.5 -1.2));% John's new eqn 7/19/16
    
    %psw2_new(reef,2) = max(0.6,min(1.3,()); % John's new eqn 7/25/16
    % break up to find error:
    %middle = mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef)).^0.25/2;
    %psw2_new(reef,2) = max(0.6,min(1.3,middle));
    %{
    psw2_new(reef,1) = max(propInputValues(1),min(propInputValues(2),(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^propInputValues(3)/propInputValues(4))); % John's new eqn 7/25/16
    % Fixed values - this costs a little time on each iteration, but keeps
    % them available.
    % Next 2 were used between about 1/13 and 1/17/2017
    %psw2_new(reef,2) = max(0.35,min(1.8,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.355/2.0167)); % 1/13/17 rcp 8.5 optimized
    %psw2_new(reef,3) = max(0.35,min(1.8,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.4556/2.556)); % 1/13/17 rcp 2.6 optimized
    % Re-optimized on 1/17 to a 1985-2010 bleaching level of 3 percent.  2%
    % did not seem feasible
    psw2_new(reef,2) = max(0.35,min(2.0,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.5/2.517)); % 1/17/17 rcp 8.5 E=0 optimized
    psw2_new(reef,3) = max(0.35,min(2.0,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.5/3.104)); % 1/17/17 rcp 8.5 E=1 optimized
    psw2_new(reef,4) = max(0.35,min(2.0,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.5/2.390)); % 1/17/17 rcp 2.6 E=0 optimized
    psw2_new(reef,5) = max(0.35,min(2.0,(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^0.5/3.063)); % 1/17/17 rcp 2.6 E=1 optimized
    %}
    for j = 1:pswCount
        psw2_new(reef,j) = max(pswInputs(1, j),min(pswInputs(2, j),(mean(exp(0.063.*SSThistReef))./var(exp(0.063.*SSThistReef))).^pswInputs(3, j)/pswInputs(4, j))); % John's new eqn 7/25/16
    end
    clear SSThistReef
end

%% Save new psw2 values
cd(matPath);
save('Optimize_psw2.mat','psw2_new','Reefs_latlon','pswInputs'); %% CAL 10-3-16 based on hist SSTs!
cd(basePath);
