%% CORAL AND SYMBIONT POPULATION DYNAMICS
% MAIN FILE TO TEST PROPORTIONALITY CONSTANTS AND VARIANCE
% 14.9, 14.3, 14.8 V8 speed, everyx=1 one key reef, pdfs off
% 19.2, 18.1, 18.6 V6 speed, with two ways of computing bleaching.
% 11.7 seconds, V8 after more cleanup and removing large broadcast
% variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evolutionary model for coral cover (from Baskett et al. 2009)     %
% modified by Cheryl Logan (clogan@csumb.edu)                       %
% last updated: 5-3-16                                              %
% Performance and structural changes 9/2016 by Steve Ryan (jaryan)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SET WORKING DIRECTORY AND PATH
timerStart = tic;
Computer = 3; % 1=office; 2=laptop; 3=Steve; 4=Steve laptop; 0 = autodetect
[basePath, outputPath, sstPath, SGPath, matPath, Computer, defaultThreads] ...
    = useComputer(Computer);
cd (basePath);

%% Clear variables which I want to examine between runs, but not carry over.
clearvars bEvents resultSimilarity Omega_all Omega_factor;

%% Most-used case settings
% DEFINE CLIMATE CHANGE SCENARIO (from normalized GFDL-ESM2M; J Dunne)
RCP = 'rcp60'; % options; 'rcp26', 'rcp45', 'rcp60', 'rcp85', 'control', 'control400'
E = 1;  % EVOLUTION ON (1) or OFF (0)?
OA = 0; % Ocean Acidification ON (1) or OFF (0)?
maxReefs = 1925;  %never changes, but used below.
%% Variables for plotting, debugging, or speed testing
skipPostProcessing = false;     % Don't do final stats and plots when timing.
everyx = 20; % 1;                % run code on every x reefs, plus "keyReefs"
                                % if everyx is one of 'eq', 'lo', 'hi' it
                                % selects reefs for abs(latitude) bins [0,7],
                                % (7, 14], or (14,90] respectively.
allPDFs = false;                % if false, just prints for keyReefs.
neverPlot = false;               % For optimization runs, turn off all plotting.
saveEvery = 5000;               % How often to save stillrunning.mat. Not related to everyx.
saveVarianceStats = false;      % Only when preparing to plot selV, psw2, and SST variance.

% Super symbiont options
newMortYears = 0; % If true, save a fresh set of "long mortality" years based on this run.
superMode = 0;  % 0 = add superAdvantage temperature to standard "hist" value.
                % 1 = use mean of temperatures in the superInitYears time range.
                % 2 = use max of temperatures in the superInitYears time range.
                %     Modes 0 to 2 all start at 2035.
                % 3 = use mean of last 10 years before the end of the first 5-year mortality.
                % 4 = use max of last 10 years before the end of the first 5-year mortality.
                % 5 = use fixed delta like option 0, but start as in
                %     options 3 and 4.
                % 6 = use fixed delta like option 0, but start according to
                %     first year of bleaching.

superAdvantage = 0.123;           % Degrees C above native symbionts.
startSymFractions = [1.0 0.0];  % Starting fraction for native and super symbionts.

% If this code is called from a script, we want some of the variables above
% to be overwritten.  Do that here before they are used in the code below.
if exist('scriptVars', 'var')
    disp('Called by a script');
    % Only these variables are supported now, but it's easy to add more.
    if isfield(scriptVars, 'RCP') RCP = scriptVars.RCP; end;
    if isfield(scriptVars, 'E') E = scriptVars.E; end;
    if isfield(scriptVars, 'OA') OA = scriptVars.OA; end;
    if isfield(scriptVars, 'superMode') superMode = scriptVars.superMode; end;
    if isfield(scriptVars, 'superAdvantage'); superAdvantage = scriptVars.superAdvantage; end;
    fprintf('After update in model E = %d, OA = %d, RCP = %s Mode = %d Adv = %d \n', E, OA, RCP, superMode, superAdvantage);
end


if superMode >= 3 && superMode <=5
    fn = strcat('longMortYears_', RCP, '_', num2str(E));
    load(fn, 'longMortYears');
    suppressSuperUntil = longMortYears;
    % XXX - next two lines for testing only!!!
    %subFrom = (suppressSuperUntil > 1882);
    %suppressSuperUntil = suppressSuperUntil - 20*subFrom;
    assert(length(suppressSuperUntil) == maxReefs, 'By-year symbiont starts should have a start year for every reef.');
elseif superMode >= 6
    fn = strcat('firstBleachYears_', RCP, '_', num2str(E));
    load(fn, 'firstBleachYears');
    suppressSuperUntil = firstBleachYears;
    assert(length(suppressSuperUntil) == maxReefs, 'By-year symbiont starts should have a start year for every reef.');
else
    suppressSuperUntil = 2035*ones(maxReefs, 1); % Until this date set S(:, 3:4) to zero.  Use 1861 to start from the beginning.
end
%superInitYears = [2025 2035]; % Years from which adaptation temp is selected.
superSeedFraction = 10^-4;      % Fraction of S_seed to use for seeding super symbionts.
% sizing notes: for massive,
% KSm = 3000000
% seed = 100000
% KSm/seed = 30
% the seed fraction is a fraction of this value.
% for 0.01 -> KSm/introduced = 3000
% for 0.0001 -> KSm/introduced = 300000
oneShot = true;  % After supersymbiont introduction, set its seed to zero.
assert((length(suppressSuperUntil) == 1 && suppressSuperUntil<=1861)  || startSymFractions(2) == 0, 'Start with no symbionts if they are to be suppressed at first.');
assert(sum(startSymFractions) == 1.0, 'Start fractions should sum to 1.');

%% Specify a set of reefs for diagnostic plotting and set up parameters for the plot routine.
%probreefs = [24,49,50,55,60,70,86];  % used in the past
% Moorea, Curacao, St John VI, Ko Phuket, Australia (-16.92, 145.78)
% 144     402      420         793        1541  % Mentioned in Baskett
% 2009, but there it's Heron Island, which is at -23.4, -151.9
% 106  Jarvis Island, the first to show mortality.
% Others to test based on variance...       lat     lon     Area
% Reef 610, highest variance = 26.269       26.52,  52.5    Persian Gulf
% Reef 1463, lowest variance = 0.32899      -1.94,  135.5   NW New Guinea
% Reef 1638, coldest mean = 20.108          -34.5,  151.5   SE of Sydney
% (is this right? Lord Howe Island is said to have the southernmost
% fringing reef at -31.5 S.
% Reef 1354, hottest mean = 29.617           0.50,  128.5   North Maluku, Indonesia
%
% Galapagos
% San Cristobal -0.8574, -89.4364
% 246 [-89.5000000000000,-0.843542042560123] is the closest stored value is it the cell center?
% Wolf 1.38159, -91.816438
% 234 [-91.5000000000000,1.56020553279879]
% Darwin 1.678434, -92.003086
% 225 [-92.5000000000000,1.56020553279879]  note that cell 234 is equally
% close since 9.00 is halfway between.
%
% Rarotonga
% 103 [
%
%keyReefs = [106 144 402 420 610 793 1354 1463 1541 1638];
% Reefs that give good matches between 12/13 and original bleaching: 951, 1001
% Reefs with the most bleaching: 952, 968
% Reefs that are terrible matches: 106 and 610 (already listed) are tied with many others
%keyReefs = [1:10:maxReefs 30 130 106 144 402 420 610]; %, 130];
%keyReefs = [130 71 511 1691 391 1301 1061 106 144 402 420 610];
%keyReefs = [maxReefs 130 71 1691 1301 106 402 420 610 793 1354 1463 1541 ];
%keyReefs = 1:maxReefs;
% Large set to look at
%keyReefs = [71 106 144 261 402 420 511 521 610 793 1301 1354 1463 1541 1638 1691 1701];
%keyReefs = [keyReefs [99 420 1738]];
%keyReefs = [402 420];
%keyReefs = [225 230 231 232 233 234 238 239 240 241 244 245 246 247 248];
%keyReefs = [610 1463];
keyReefs = [1 103];

% Reefs with the earliest mortality in the rcp85, E=1 case are listed below.  All
% experi3.53ence 5 years of mortality by 2012.
%  29  30  47  48  68  69  98  99 482 1903 1904
% These reefs reach the same condition in 2099
% 199         990        1738

% Put them in order without duplicates.  Shouldn't matter, but why not?
keyReefs = unique(keyReefs);
% too good to be true matches: keyReefs = [561 574 762 1271 1298 1299 1325 1326 1228 402];
%keyReefs = []; % can be set empty to minimize plotting and other work for some tests


%% Use of parallel processing on the local machine.
% no argument: uses the existing parallel pool if any.
% 0:  runs sequentially.
% > 0 tries to start the requested number of workers, failing if the value
%     is greater than the value specified in Matlab's Parallel Preferences.
cd (basePath)

if exist('useTestThreads', 'var')
    % Scripts can provide this variable to compare speed versus threads.
    [multiThread, queueMax] = parallelSetup(useTestThreads);
else
    % defaultThreads is set in useComputer for each user's computers.
    % 1 is not normally used.  0 means to run sequentially.
    [multiThread, queueMax] = parallelSetup(defaultThreads);
end

if multiThread
    queue =  parallel.FevalFuture.empty;
    % Call plot function with no arguments.  If the workers have been used
    % before this clears their variables.  If not it gets the plot routine
    % loaded onto them.
    spmd
        %Plot_SST_Decimate();
        Plot_SSTDensityCover_Decimate();
        graphCompare();
    end
end
fprintf('Starting (after parallel setup) at %s\n', datestr(now));

%% Less frequently changed model parameters
dt = 0.125;  % The fraction of a month for R-K time steps

% SST DATASET?
Data = 1; % 1=ESM2M_norm;  2=HADISST (through 3/16/16)
if Data == 1
    dataset = 'ESM2M';
else
    dataset = 'HadISST';
end

% NORMALIZATION FACTOR
NF = 1 ; % modify reach 10% bleaching freq bn 1985-2010

% Selection of variance column from psw2_new
% As of 1/13/17, column 1 is for optimization
% 2 - rcp 8.5 E0 cases, opt for bleaching = 3.0
% 3 - rcp 8.5 E1 cases, opt for bleaching = 3.0
% 4 - rcp 2.6 E0 cases, opt for bleaching = 3.0
% 5 - rcp 2.6 E1 cases, opt for bleaching = 3.0
% for  5% target use 10 to 13
% for 10% target use 14 to 17.
% after 2/4/17 use 20 to 23 for 5% target.
% Choose the correct match according to the E and RCP values.
if exist('optimizerMode', 'var')
    propTest = 1; % 16 was used before optimization and code changes
elseif E == 0 && (strcmp(RCP, 'rcp26') || strcmp(RCP, 'control400'))
        propTest = 20; % 6;
elseif E == 0 && strcmp(RCP, 'rcp85')
        propTest = 21; % 7;  % 2;
elseif E == 1 && (strcmp(RCP, 'rcp26') || strcmp(RCP, 'control400'))
        propTest = 22; % 8;
elseif E == 1 && strcmp(RCP, 'rcp85')
        propTest = 23; % 9;
elseif E == 0 && strcmp(RCP, 'rcp45')
        propTest = 24; % 9;
elseif E == 0 && strcmp(RCP, 'rcp60')
        propTest = 25; % 9;
elseif E == 1 && strcmp(RCP, 'rcp45')
        propTest = 26; % 9;
elseif E == 1 && strcmp(RCP, 'rcp60')
        propTest = 27; % 9;
else
    error('No propTest value is defined for E = %d and RCP = %s\n', E, RCP);
end

initYear = '2001';  % SelV and hist will be initialized from 1861 to this year.
%% Parameters to the new bleaching model (massive, branching)
%  Original parameters are shown after the comments.  Current defaults
%  (used whenever optimization is not running) are assigned to each
%  variable.
% Values in the comments to the right were used for 4PM 1/9/2017 outputs.
% Keep for comparison.
% BLEACHING
bleachParams.sBleach = [0.3 0.3];         % Dropping to this value is bleaching   [0.22 .32]
bleachParams.cBleach = [0.1 0.1];          % Dropping to this value is bleaching   [0.1 0.1]
% New on 1/10/17
bleachParams.sRecoverySeedMult = [4 4];     % Required for recovery.
bleachParams.cRecoverySeedMult = [4 4];     % Required for recovery. Should be greater than cSeedThresholdMult
% MORTALITY
bleachParams.cSeedThresholdMult = [2 2];   % Seed multiplier for mortality (a)     [3 20]
bleachParams.yearsToMortality = 5;          % Years of continuous bleaching before mortality.   5


% Above values may be overridden by the optimizer.
if exist('optimizerMode', 'var') && exist('optimizerBleachParams', 'var')
    % Only change the parts of the struct which are given.
    bleachParams = updateIfGiven(bleachParams, optimizerBleachParams);
end



%% Some useful paths and test strings to be used later.
format shortg; c = clock; dateString = strcat(num2str(c(1)),num2str(c(2),'%02u'),num2str(c(3),'%02u')); % today's date stamp
pdfDirectory = strcat(dataset, RCP,'_E',num2str(E), '_OA',num2str(OA),'_NF',num2str(NF),...
    '_sM',num2str(superMode),'_sA',num2str(superAdvantage),'_',dateString,'_figs/');
%pdfDirectory = strcat(dataset, RCP,'_E',num2str(E), '_OA',num2str(OA),'_NF',num2str(NF),'_',dateString,'_figs/');
mkdir(outputPath, pdfDirectory);

mkdir(strcat(outputPath, 'bleaching'));
% Map directory is used for maps, but also for miscellaneous text and
% plotted output, since it is the least cluttered of the directories.
mapDirectory = strrep(pdfDirectory, '_figs', '_maps');
mkdir(outputPath, mapDirectory);
fullMapDir = strcat(outputPath, mapDirectory);

% Initialize a file for logging most of what goes to the console.
echoFile = fopen(strcat(outputPath, mapDirectory, 'console.txt'), 'w+');
logTwo(echoFile); % Required first call to set output path.

%% LOAD JOHN'S NORMALIZED SSTS FROM EARTH SYSTEM CLIMATE MODEL OR HADISST
% Extract SSTs for a ALL reef grid cells
[SST, Reefs_latlon, TIME, startYear] = GetSST_norm_GFDL_ESM2M(sstPath, matPath, Data, RCP);
lenTIME = length(TIME);

%% LOAD Omega (aragonite saturation) values if needed

if OA == 1
    [Omega_all] = GetOmega(SGPath, RCP);
    % just for debugging: Omega_all = 5.0*ones(maxReefs, 2880);
    % Convert omegas to growth-factor multipliers so there's
    % less logic inside the time interations.
    [Omega_factor] = omegaToFactor(Omega_all);
else
    % Wasteful to make a big empty array, but it makes entering the
    % parallel loop simpler.  Note that only the last value is set.
    %Omega_all(maxReefs, lenTIME) = 0.0;
    Omega_factor(maxReefs, lenTIME) = 0.0;
end

%% SUB-SAMPLE REEF GRID CELLS [%DIFFERs FROM PREV VERSION]
numTestReefs = length(Reefs_latlon);
% Since just iterating with "everyx" won't hit all keyReefs, build a list
% of reefs for the current run.
if isnumeric(everyx)
    toDo = 1:everyx:numTestReefs;   % as specified by everyx
else
    % everyx can specify a reef area
    if strcmp(everyx, 'eq')
        toDo = find(abs(Reefs_latlon(:, 2))<=7)';
    elseif strcmp(everyx, 'lo')
        toDo = find(abs(Reefs_latlon(:, 2))<=14 & abs(Reefs_latlon(:, 2)) > 7)';
    elseif strcmp(everyx, 'hi')
        toDo = find(abs(Reefs_latlon(:, 2))>14)';
    else
        disp('WARNING: everyx was not a number or one of a few allowed strings.  Using 1.');
        toDo = 1:1:numTestReefs;
    end
end
toDo = unique([toDo keyReefs]); % add keyReefs defined above
% Before adding keyReefs, the code used length(Reefs_latlon)/everyx
reefsThisRun = length(toDo);
lastReef = toDo(end);
% Note that this sizes the array for all reefs, even when doing a subset -
% better for parfor (maybe?).
logTwo('Modeling %d reefs.\n', reefsThisRun);


%% INITIALIZE MATRICES AND DEFINE 'NORMALIZATION FACTOR' (NF)

% TODO see if these could be initialized to reefsThisRun to save memory:
Mort_85_10_all = nan(numTestReefs,1);
CoralCover_Branch_2010 = nan(numTestReefs,1);


%% LOAD SELECTIONAL VARIANCE (psw2)
load (strcat(matPath, 'Optimize_psw2.mat'),'psw2_new', 'pswInputs')
% pswInputs are not used in computations, but they are recorded to document
% each run.
pswInputs = pswInputs(:, propTest);
% XXX this changes results - don't do it! psw2_new = psw2_new(:, propTest);  % no need to carry columns we never use.
%% Constants moved outside of the loop so they are only defined once:
C_seed = [0.1*10^7 0.01*10^7]; % mass branch
S_seed = [0.1*10^6 0.1*10^6 0.1*10^6 0.1*10^6];  % This was a single value until 2/17/2017
% A string frequently used in titles and file names:
modelChoices = strcat(dataset, RCP,'.E',num2str(E),'.OA',num2str(OA));

% Load .mat file for Coral and Symbiont genetics constants
% As of 1/3/2017 the file contains a structure rather than individual
% variables.
load(strcat(matPath, 'Coral_Sym_constants_2.mat')); % default is evolution OFF
% XXX override Sn for testing.
%coralSymConstants.Sn = 1;
assert(length(startSymFractions) == coralSymConstants.Sn, 'Symbiont start fractions should match number of symbionts.');


% Mutational Variance w (E=1) and w/o Evolution (E=0)
if E==0; vM = 0 ;     % Mutational variance (degC^2/yr) (convert to months/12)
else vM = coralSymConstants.ve*.001/12; end %%(1.142*10^-5)/12 ;    % Mutational variance (degC^2/yr) (convert to months/12)
vMT   = vM;                    % Mutational variance (degC^2/yr) (convert to months/12)
MutV  = [vM vM];               % Mutational variance matrix for symbiont calcs
MutVx = repmat(MutV,1,coralSymConstants.Sn);     % Mutational variance matrix for coral calcs
% January 2016, more variables need replication when Sn > 1
coralSymConstants.EnvVx = repmat(coralSymConstants.EnvV,1,coralSymConstants.Sn);     % Environmental variance matrix
%XXX
coralSymConstants.KSx = repmat(coralSymConstants.KS,1,coralSymConstants.Sn);     % Environmental variance matrix



%% Set up indexing and time arrays before entering the main loop.

% Moved from Interp_data since it doesn't change between reefs.
% Requires a dummy read of SSTHist for sizing.
SST_LOC = SST(1, :);                     % Reef grid cell location
SSThist = SST_LOC';

% NOTE: the next 7 lines of code just create a duplicate of TIME,
% called tim.  Why?
months = length(SSThist);
assert(mod(months, 12) == 0, 'Calculations assume a time span of whole years.');
years = months/12;
stepsPerYear = 12/dt;
% Array approach to creating the list of dates is over 100 times faster
yyy = repmat(startYear:startYear+years-1, 12, 1); % all the years, repeated in 12 rows
mmm = repmat(1:12, years, 1)';  % months 1-12 repeated for each year
ddd = 15*ones(12, years);        % day 15
tim = datenum(yyy, mmm, ddd);   % the big saving is here
tim = tim(:)';                  % turn from a 2D array to 1D as needed below.

time = interp(tim,1/dt,1,0.0001)'; % Interpolate time 4X times (8X when dt = 0.125)
% Set index for mean historical Temp between 1861-2000 (ESM2M_historical; yrs used in Baskett et al 2009)
% Note: original code had a different end point for Data=2 than Data=1.
initIndex = findDateIndex(strcat('30-Dec-', initYear), strcat('31-Dec-', initYear), time);
if length(suppressSuperUntil) > 1
    % This could be done more efficiently...
    % suppressSuperUntil units are years.
    for i = length(suppressSuperUntil):-1:1
        ssY = suppressSuperUntil(i);
        if ssY == 0
            suppressSuperIndex(i) = 0;
        else
            suppressSuperIndex(i) = findDateIndex(strcat('14-Jan-', num2str(ssY)), strcat('16-Jan-',num2str(ssY)), time);
        end
    end

else
    error('Code no longer gets here, right?');
    % This is not reached currently.
    %{
    suppressSuperIndex = findDateIndex(strcat('14-Jan-', num2str(suppressSuperUntil)), strcat('16-Jan-',num2str(suppressSuperUntil)), time);
    superInitRange(1) = findDateIndex(strcat('14-Jan-', num2str(superInitYears(1))), strcat('16-Jan-',num2str(superInitYears(1))), time);
    superInitRange(2) = findDateIndex(strcat('14-Dec-', num2str(superInitYears(2))), strcat('16-Dec-',num2str(superInitYears(2))), time);
    %}
end
% max so it's always a valid index, BUT note that suppressSuperIndex
% can be zero when a super symbiont is never needed!
suppressSuperIndexM10 = max(1, suppressSuperIndex - 10*stepsPerYear);

% Same for the SST and TIME arrays, but they are coarser, with dates at the
% 15th of each month.
initSSTIndex = findDateIndex(strcat('14-Dec-', initYear), strcat('16-Dec-',initYear), TIME);
iteratorHandle = selectIteratorFunction(length(time), Computer);
%XXX DEBUG:
% timeSteps = 96;
timeSteps = length(time) - 1;      % number of time steps to calculate

% This should cause the "parfor" to run as serial code when plotting
% more than one reef per plot.  The single-reef plot options work fine in
% parallel.  Note that this does NOT remove the overhead of copying arrays
% for the parallel code.
% Also override queueMax so the arrays are handled correctly.

%% Split up reefs into batches for parallel computation.  With one core
%  specified, it simply uses one batch.
[parSwitch, queueMax, chunkSize, toDoPart] = parallelInit(queueMax, toDo);

% Several arrays are built in the parallel loop and then used for
% later analysis.  Parfor doesn't like indexing into part of an array.  The
% trick is to make a local array for each iteration inside the parfor, and
% then assemble them into the desired shape afterwards.  Note: all the "nan(1,1)"
% entries are there because the parfor needs to have the "empty" output
% arrays defined before the loop.  Instead of creating the contents at full
% size and passing big arrays of nan to the workers, it make more sense to
% pass these dummy arrays and allocate the working memory in each worker.
for i = queueMax:-1:1
    resultSim_chunk{i} = nan(1,1);
    % ACM_chunk{i} = nan(1, 1);
    % ACB_chunk{i} = nan(1, 1);
    bEvents_chunk{i} = [];
    LatLon_chunk{i} = Reefs_latlon(min(toDoPart{i}):max(toDoPart{i}),1:2);
    %suppressYears_chunk{i} = suppressSuperUntil(min(toDoPart{i}):max(toDoPart{i}));
    SST_chunk{i} =SST(min(toDoPart{i}):max(toDoPart{i}), :);
    Omega_chunk{i} = Omega_factor(min(toDoPart{i}):max(toDoPart{i}), :);
    kOffset(i) = min(toDoPart{i});
    C_cum_chunk{i} = zeros(length(time), coralSymConstants.Sn*coralSymConstants.Cn); % Sum coral cover for all reefs.
    % 3D array of time by reef by coral type.  Note that we don't care
    % about the identity of the reefs in this case, so we just need enough
    % columns for all reefs actually calculated, ignoring those which are
    % skipped.
    C_year_chunk{i} = zeros(years, chunkSize, coralSymConstants.Cn); % Coral cover for all reefs, but just 2 columns.
    Massive_dom_chunk{i} = zeros(length(time), 1);

end
% TODO input arrays such as SST and Reefs_latlon are sent at full size to
% each worker.  Consider sending just the correct subset to each.

%% RUN EVOLUTIONARY MODEL
%parfor (parSet = 1:queueMax, parSwitch)
for parSet = 1:queueMax
    %  pause(1); % Without this pause, the fprintf doesn't display immediately.
    %  fprintf('In parfor set %d\n', parSet);
    reefCount = 0;
    % Variables to collect and return, since parfor won't allow direct
    % insertion into an output array.
    % TODO see if length(toDoPart(parSet)) should be used instead of
    % chunkSize.
    par_bEvents = [];
    par_SST = SST_chunk{parSet};
    par_Omega = Omega_chunk{parSet};
    par_LatLon = LatLon_chunk{parSet};
    %par_SuppressYears = suppressYears_chunk{parSet};
    par_kOffset = kOffset(parSet);
    par_C_cum = C_cum_chunk{parSet};
    par_C_year = C_year_chunk{parSet};
    par_Massive_dom = Massive_dom_chunk{parSet};
    par_HistSuperSum = 0.0;
    par_HistOrigSum = 0.0;
    par_HistOrigEvolvedSum = 0.0;
    for k = toDoPart{parSet}
        reefCount = reefCount + 1;
        SST_LOC = par_SST(1+ k - par_kOffset, :);                     % Reef grid cell location
        SSThist = SST_LOC';                      % Transpose SST matrix
        Omega_hist = par_Omega(1+ k - par_kOffset, :);
        psw2 = NF*psw2_new(k, propTest) ;         % UPDATED** max(0.3,min(1.3,(mean(exp(0.063.*SSThist))./var(exp(0.063.*SSThist))).^0.5/7)); % John's new eqn 8/10/16** try this
        reefLatlon = par_LatLon(1 + k - par_kOffset, :);
        lat = num2str(round(reefLatlon(2)));
        lon = num2str(round(reefLatlon(1)));
        LOC = strcat('_', num2str(lat),'_',num2str(lon),'_');

        % Interpolate data and create time stamp
        %% Set timestep and interpolate temperature, omega, and time stamp
        temp = interp(SSThist,1/dt); % Resample temp 4X times higher rate using lowpass interpolation
        omega = interp(Omega_hist, 1/dt);
        % new vector is 4x length of orginal


        %% Make Histogram of psw2 and map of var(SSThist)
        % hist_prop_fig     % run sub m-file to make map & histogram

        % Use a limited range of SSTHist for selectional variance, so
        % that we don't include modern climate change.
        SelV = [1.25 1]*psw2*var(SSThist(1:initSSTIndex));
        SelVx = repmat(SelV,1,coralSymConstants.Sn);     %#ok<PFBNS> % Selectional variance matrix for coral calcs
        %SelV = [1.25 1]*psw2*var(SSThist_anom(:))

        % Initialize symbiont genotype, sym/coral population sizes, carrying capacity

        %ssss = findDateIndex(strcat('14-Jan-', num2str(par_SuppressYears(1+ k - par_kOffset)-10)), strcat('16-Jan-',num2str(par_SuppressYears(1+ k - par_kOffset)-10)), time);
        %eeee = findDateIndex(strcat('14-Dec-', num2str(par_SuppressYears(1+ k - par_kOffset))), strcat('16-Dec-',num2str(par_SuppressYears(1+ k - par_kOffset))), time);

        [vgi, gi, S, C, hist, ri] = Init_genotype_popsize(Data, time, initIndex, temp, coralSymConstants, ...
            E, vM, SelV, superMode, superAdvantage, startSymFractions, ...
            [suppressSuperIndexM10(k) suppressSuperIndex(k)], k);



        if length(suppressSuperIndex) > 1
            ssi = suppressSuperIndex(k);
        else
            ssi = suppressSuperIndex;
        end
        %fprintf('super will start at index %d\n', ssi);
        %% MAIN LOOP: Integrate Equations 1-5 through 2100 using Runge-Kutta method
        [S, C, ri, gi, vgi, origEvolved] = iteratorHandle(timeSteps, S, C, dt, ri, temp, OA, omega, vgi, gi, ...
                    MutVx, SelVx, C_seed, S_seed, ssi, ...
                    superSeedFraction, oneShot, coralSymConstants); %#ok<PFBNS>
        % These, with origEvolved, compare the average native and
        % supersymbiont genotypes with the evolved state of the native
        % symbionts just before the supersymbionts are introduced.
        origHist = gi(1,1);
        superHist = gi(1,3);

        par_HistSuperSum = par_HistSuperSum + superHist;
        par_HistOrigSum = par_HistOrigSum + origHist;
        par_HistOrigEvolvedSum = par_HistOrigEvolvedSum + origEvolved;

        if any(keyReefs == k)  % temporary genotype diagnostic
            %{
            suff = '';
            if superMode && superMode ~= 5
                suff = sprintf('_%s_E%d_SymStrategy%d_Reef%d', RCP, E, superMode, k);
            elseif superMode == 0 || superMode == 5
                suff = sprintf('_%s_E%d_SymStrategy%dAdv%0.2fC_Reef%d', RCP, E, superMode, superAdvantage, k);
            end
            genotypeFigure(fullMapDir, suff, k, time, gi, ssi);

            % Growth rate vs. T as well
            % TODO: dies when ssi = 0...
            growthRateFigure(fullMapDir, suff, datestr(time(ssi), 'yyyy'), ...
                k, temp, gi, vgi, ssi, ...
                coralSymConstants, SelVx);

           %}
        end

        par_C_cum = par_C_cum + C;
        par_Massive_dom = par_Massive_dom + C(:, 1) > C(:, 2);
        % Time and memory will be consumed, but we need stats on coral
        % cover.
        par_C_year(:, reefCount, 1) =  decimate(C(:, 1), stepsPerYear, 'fir');
        par_C_year(:, reefCount, 2) =  decimate(C(:, 2), stepsPerYear, 'fir');
        
        %% New clean stats section
        
        % Call to check for fatal errors, but don't use the outputs yet.
        [ C_monthly, S_monthly, C_yearly, bleachEvent ] = Clean_Bleach_Stats(C, S, C_seed, S_seed, dt, TIME, bleachParams, coralSymConstants);
        
        plotFlag = ~neverPlot && (allPDFs || any(keyReefs == k));     % may only be plotting keyReefs
        if plotFlag
            % Now that we have new stats, reproduce the per-reef plots.
            Plot_One_Reef(C_monthly, S_monthly, bleachEvent, psw2, time, temp, lat, lon, RCP, ...
                  hist, Data, initIndex, gi, ri, SGPath, outputPath, k, ...
                  pdfDirectory, LOC, NF, E, dateString, lenTIME);
        end

        
        %% Back to pre-cleanup code
        % Testing a separate routine to get bleaching from symbiont density.
        %{
        [bleachOneReef, dCov] = Get_bleach_freq(C, C_seed, S, S_seed, k, time, reefLatlon, ...
                    TIME, dt, initIndex, startYear, bleachParams); %, bleachFrac);
        %}
        if ~isempty(bleachOneReef)
            par_bEvents = [par_bEvents bleachOneReef];
        end
        % GET MORTALITY FREQUENCY (bleaching estimate across all reefs)

        bleachLimit = 500;  % Just an arbitrary limit to catch extreme cases.
        if length(bleachOneReef) > bleachLimit
            % Something's wrong with this case. Send a warning and get out.
            bleachParams
            length(bleachOneReef)
            k
            error('ExcessiveBleaching.  Aborting one parallel set.  Reef %d has over %d bleaching events.', k, bleachLimit);
        end

        % MAKE PLOTS (SST, Symbiont Genotype, Symbiont Density, Coral Cover)

        % Plot the current reef - call is different when using threads.
        plotFlag = ~neverPlot && (allPDFs || any(keyReefs == k));     % may only be plotting keyReefs
        if plotFlag
            %Plot_SST_Decimate(psw2, time, temp, lat, lon, RCP, ...
            %  hist, C, S, Data, initIndex, gi, dCov, ri, basePath, outputPath, k, ...
            %  pdfDirectory, LOC, NF, E, dateString, lenTIME);
            %{
            Plot_SSTDensityCover_Decimate(psw2, time, temp, lat, lon, RCP, ...
              hist, C, S, Data, initIndex, gi, dCov, ri, SGPath, outputPath, k, ...
              pdfDirectory, LOC, NF, E, dateString, lenTIME);
            %}
        end

        printFreq = max(10, ceil(length(toDoPart{parSet})/4)); % The last digit is the number of pieces to report.
        if parSwitch && mod(reefCount, printFreq) == 0
            fprintf('Set %d is %3.0f percent complete.\n', parSet, (100*reefCount/length(toDoPart{parSet})));
        end

    end % End of reef areas for one parallel chunk

    bEvents_chunk{parSet} = par_bEvents;
    C_cum_chunk{parSet} = par_C_cum;
    C_year_chunk{parSet} = par_C_year(:, 1:reefCount, :);
    Massive_dom_chunk{parSet} = par_Massive_dom;
    histSuper_chunk(parSet) = par_HistSuperSum;
    histOrig_chunk(parSet) = par_HistOrigSum;
    histOrigEvolved_chunk(parSet) = par_HistOrigEvolvedSum;


end % End of parfor loop
% Build these variables from the chunks.
clearvars SST_chunk;
bEvents = horzcat(bEvents_chunk{:});
clearvars bEvents_chunk; % release some memory.

C_yearly = horzcat(C_year_chunk{:});
% Total coral cover across all reefs, for ploting shift of dominance.
C_cumulative = zeros(length(time), coralSymConstants.Sn*coralSymConstants.Cn);
Massive_dom_cumulative = zeros(length(time), 1);
superSum = 0.0;
histSum = 0.0;
histEvSum = 0.0;
for i = 1:queueMax
    C_cumulative = C_cumulative + C_cum_chunk{i};
    Massive_dom_cumulative = Massive_dom_cumulative + Massive_dom_chunk{i};
    superSum = superSum + histSuper_chunk(i);
    histSum = histSum + histOrig_chunk(i);
    histEvSum = histEvSum + histOrigEvolved_chunk(i);
end
clearvars C_cum_chunk;
superSum = superSum/reefsThisRun;
histSum = histSum/reefsThisRun;
histEvSum = histEvSum/reefsThisRun;
logTwo('Super symbiont genotype = %5.2f C.  Base genotype %5.2f C (advantage %5.2f), Evolved base %5.2f (advantage %5.2f).\n', ...
    superSum, histSum, (superSum-histSum), histEvSum, (superSum-histEvSum));

if saveVarianceStats
    assert(length(toDo) == maxReefs, 'Only save variance data when running all reefs!');
    % Save selectional variance and last year of cover for binned plotting
    % by case.  Note that these numbers are computed inside the parallel
    % loop, but it's easier to recompute them here than to build and
    % extract arrays from the parallel code.
    selVariance(maxReefs) = 0;
    tVariance(maxReefs) = 0;
    for k = 1:maxReefs
        SSThist = SST(k, :);
        tVariance(k) = var(SSThist(1:initSSTIndex));
        selVariance(k) = psw2_new(k)*tVariance(k);
    end

    lastYearOfCover = GetLastYearArray(bEvents, maxReefs);
    save(strcat(basePath, 'LastYear', '_selV_', RCP, 'E=', num2str(E), 'OA=', num2str(OA), '.mat'), 'psw2_new', 'selVariance', 'tVariance', 'lastYearOfCover', 'RCP', 'OA', 'E');
end

if ~skipPostProcessing

    justSummaries = bEvents(strcmp({bEvents.event}, 'BLEACH8510') & strcmp({bEvents.coral}, 'REEF'));
    Bleaching_85_10 = 100 * sum(cat(1, justSummaries.count)) / (reefsThisRun*(1+2010-1985));
    % NEW 1/10/2017 count bleaching _events_
    bleachEvents = bEvents(strcmp({bEvents.event}, 'BLEACH') & [bEvents.year] >=1985 & [bEvents.year] <= 2010);
    count852010 = length(bleachEvents);
    Bleaching_85_10_By_Event = 100*count852010/reefsThisRun/(2010-1985+1);
    fprintf('Bleaching by duration = %6.4f and by event = %6.4f\n', ...
        Bleaching_85_10, Bleaching_85_10_By_Event);

    if strcmp(RCP, 'control400')
        PRGYears = 1950:25:2250;
    else
        PRGYears = [1950 2000 2050 2100 2016];
    end
    ct = 1;
    % Last modeled year, so we don't count a coral which lives to then as
    % dead.
    lastYear = str2double(datestr(TIME(end), 'yyyy'));



    format shortg;
    filePrefix = strcat(modelChoices,'_NF',num2str(NF),'_',dateString);
    % Don't save all this data if we're just optimizing.  TODO: better name
    % for the boolean.
    if ~neverPlot
        % fname = strcat(filePrefix,'.mat');
        fname = strcat(outputPath, pdfDirectory, filePrefix, '.mat');
        save(fname, 'toDo', ...
            'E','OA','pdfDirectory','dataset', ...
            'Reefs_latlon','everyx','CoralCover_Branch_2010','NF','RCP','numTestReefs');

        % Get the year range for last year maps first so both functions use
        % the same range.
        lastFull = bEvents([bEvents.last]==1 & strcmp({bEvents.coral},'REEF') & strcmp({bEvents.event}, 'MORTALITY'));
        if isempty(lastFull)
            rangeBoth = [1960 2100];
        else
            rangeNew = [min([lastFull.year]) max([lastFull.year])];
            % Plotting chokes if the values are equal.
            if rangeNew(1) == rangeNew(2)
                rangeNew(2) = rangeNew(2) + 1;
            end
            rangeBoth = rangeNew;
        end

        MapsCoralCoverNew(fullMapDir, bEvents, rangeBoth, modelChoices, filePrefix);


        % New dominance graph
        % Get stats based on C_yearly - try getting quantiles per row.
        % C_yearly has year/reef/coral type
        % C_quant will have year/quantiles/coraltype
        C_quant = quantile(C_yearly, [0.1 0.25 0.5 0.75 0.9], 2);
        C_quant(:, :, 1) =  100 * C_quant(:, :, 1) / coralSymConstants.KCm;
        C_quant(:, :, 2) =  100 * C_quant(:, :, 2) / coralSymConstants.KCb;
        % Create vertexes around area to shade, running left to right and
        % the right to left.
        % 5% and 95%
        span = [startYear:startYear+years-1]';

        if superMode && superMode ~= 5
            suffix = sprintf('_%s_E%d_SymStrategy%d', RCP, E, superMode);
        elseif superMode == 0 || superMode == 5
            suffix = sprintf('_%s_E%d_SymStrategy%dAdv%0.2fC', RCP, E, superMode, superAdvantage);
        end
        coralCoverFigure(fullMapDir, suffix, [span;flipud(span)], ...           % x values for areas
            [C_quant(:,1,1);flipud(C_quant(:,5,1))], ...    % wide quantiles, massive
            [C_quant(:,1,2);flipud(C_quant(:,5,2))], ...    % wide quantiles, branching
            [C_quant(:,2,1);flipud(C_quant(:,4,1))], ...    % narrow quantiles, massive
            [C_quant(:,2,2);flipud(C_quant(:,4,2))], ...    % narrow quantiles, branching
            span, squeeze(C_quant(:, 3, 1:2))); %  values for lines

        %{
        works for just average cover...
        figureCover = figure('Name','Global Coral Cover');
        % Massive +- the percentage of reefs not in the majority.
        % Scale by carrying capacity, just using the first two columns.
        % We want a percentage of capacity, so also divide by reefs.
        C_cumulative(:, 1) =  100 * C_cumulative(:, 1) / coralSymConstants.KCm / reefsThisRun;
        C_cumulative(:, 2) =  100 * C_cumulative(:, 2) / coralSymConstants.KCb / reefsThisRun;
        dy = Massive_dom_cumulative/reefsThisRun; % fraction massive dominant
        dy = min(dy, 1-dy);  % fraction minority dominant
        dy = 5*dy; % arbitrary magnifier - a few percent is nearly invisible.
        fill([time;flipud(time)], [C_cumulative(:,1).*(1-dy);flipud(C_cumulative(:,1).*(1+dy))],[1.0 0.6 0.6], 'linestyle', 'none');
        hold on;
        fill([time;flipud(time)], [C_cumulative(:,2).*(1-dy);flipud(C_cumulative(:,2).*(1+dy))],[0.6 0.6 1.0], 'linestyle', 'none');
        plot(time, C_cumulative(:, 1), '-r');
        plot(time, C_cumulative(:, 2), '-b');
        yyaxis right
        plot(time, 100*Massive_dom_cumulative/reefsThisRun);
        %}
    end
    % Note that percentMortality is not used in normal runs, but it is
    % examined by the optimizer when it is used.
    [percentBleached, percentMortality] = New_Stats_Bleach(bEvents, toDo, Reefs_latlon, ...
        outputPath, bleachParams, RCP, E, OA);

    % Get the years when reefs first experienced lasting mortality.
    % This isn't wanted every run, and certainly not when super symbionts
    % are introduced in a varable way.
    if superMode == 0 && newMortYears
        lastingEvent = bEvents(strcmp({bEvents.event}, 'LONGMORT'));
        longMortYears = [lastingEvent(:).year];
        fn = strcat('longMortYears_', RCP, '_', num2str(E));
        save(fn, 'RCP', 'E', 'OA', 'longMortYears');
        MapsSymbiontYears(fullMapDir, modelChoices, filePrefix, longMortYears, Reefs_latlon);
        quants = [0.005 .05 .25 .5 .75 .95 0.995];
        qLong = quantile(longMortYears, quants);
        logTwo('Quantile fractions for the new longMortYears array:\n');
        logTwo('Fraction  Year\n');
        fmt = repmat('%8.3f', 1, length(quants));
        logTwo(strcat(fmt,'\n'), quants);
        fmt = repmat('%8d', 1, length(quants));
        logTwo(strcat(fmt,'\n'), qLong);

        % A second set of flags, using the first year of bleaching
        % I'm assuming full-reef bleaching is wanted (both coral types),
        % but massive bleachings is a good proxy for that since branching
        % corals (always?) bleach at the same time or sooner.
        reefBleach = bEvents(strcmp({bEvents.coral}, 'MASS') & strcmp({bEvents.event}, 'BLEACH'));
        reefMort = bEvents(strcmp({bEvents.coral}, 'MASS') & strcmp({bEvents.event}, 'MORTALITY'));
        % Need to find the smallest year value for each reef.  There's
        % probably some syntax to do it in one line, but I can't think of
        % it at the moment.  This code will be called rarely, so
        % performance isn't a big issue.
        firstBleachYears(reefsThisRun) = NaN;
        for k = toDo
            yb = min([ reefBleach([reefBleach.k]==k).year  ]);
            ym = min([ reefMort([reefMort.k]==k).year  ]);
            y = 9999;
            if ~isempty(yb) || ~isempty(ym)
                if yb > 0
                    y = yb;
                end
                if ym > 0
                    y = min(y, ym);
                end
                firstBleachYears(k) = y;
            else
                firstBleachYears(k) = 0;
            end
            %fprintf('%5d x %5d x %5d x %5d \n', k, yb, ym, y);
        end
        fn = strcat('firstBleachYears_', RCP, '_', num2str(E));
        save(fn, 'RCP', 'E', 'OA', 'firstBleachYears');
    end

    logTwo('Bleaching by duration = %7.4f and by event = %7.4f\n', ...
        Bleaching_85_10, Bleaching_85_10_By_Event);
end


elapsed = toc(timerStart);
logTwo('Finished in %7.1f seconds.\n', elapsed);
logTwo('Finished at %s\n', datestr(now));
fclose('all'); % Just the file used by logTwo, normally.

%% After each run, update an excel file with descriptive information.
% 1) There seems to be no easy way to know the number of rows in the file, so
% it must be read each time.  This takes almost 1.5 seconds, even on a
% small file, so it would be best to rename the file occasionally and start
% over with a small current file.
% 2) Much of the code below is there to handle what happens when the file
% is open in Excel.  Writing is block, so user is prompted to skip the
% write or close Excel and retry.  This probably applies to any application
% using the file, not just Excel.
if ~exist('optimizerMode', 'var') || optimizerMode == false
    filename = strcat(basePath, 'RunSummaries.xlsx');
    sheet = 'Run Info';
    oldLines = 0;
    saveToExcel = 0;
    if exist(filename, 'file')
        [~, txt, ~] = xlsread(filename, sheet, 'A:A');
        oldLines = length(txt);
        unsure = 1;
        while unsure
            fid=fopen(filename,'a');
            if fid < 0
                % Construct a questdlg with three options
                choice = questdlg('RunSummaries.xlsx may be open in Excel.  Close it and try again to save results from this run.', ...
                    'File Open Conflict', ...
                    'Skip saving','Try again','Try again');
                % Handle response
                switch choice
                    case 'Skip saving'
                        disp([choice ' not saving.'])
                        unsure = 0;
                        saveToExcel = 0;
                    case 'Try again'
                        disp([choice ' re-checking file.'])
                end
            else
                fclose(fid);
                unsure = 0;
            end
        end
    end
    if saveToExcel && ~skipPostProcessing
        mat = {datestr(now), RCP, E, everyx, queueMax, elapsed, ...
                Bleaching_85_10_By_Event, ...
                bleachParams.sBleach(1), bleachParams.cBleach(1), ...
                bleachParams.sRecoverySeedMult(1), bleachParams.cRecoverySeedMult(1), ...
                bleachParams.cSeedThresholdMult(1), ...
                pswInputs(1), pswInputs(2), pswInputs(3), pswInputs(4)
            };
        if ~oldLines
            matHeader = ...
            {'Date', 'RCP', 'E', 'everyx', 'workers', 'run time', ...
                '85-2010 bleaching', ...
                'sBleach 1', 'cBleach 1', ...
                'sRecSeedMult 1', 'cRecSeedMult 1', ...
                'cSeedThreshMult 1', ...
                'pMin', 'pMax', 'exponent', 'divisor'
                };
            mat = [matHeader; mat];
        end
        range = strcat('A', num2str(oldLines+1));
        xlswrite(filename, mat, sheet, range);
    end

end










%%
% if norm == 1; print('-dpdf',strcat('~/Dropbox/Matlab/SymbiontGenetics/Figures/',LOC, '_arag_normSST_',RCP, '_prop' ,num2str(prop), '.pdf')); end
% if norm == 0; print('-dpdf',strcat('~/Dropbox/Matlab/SymbiontGenetics/Figures/',LOC, '_aragSST_',dataset, RCP, '_prop' ,num2str(prop), '.pdf')); end
