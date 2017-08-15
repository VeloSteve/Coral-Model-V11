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
Computer = 0; % 1=office; 2=laptop; 3=Steve; 4=Steve laptop; 0 = autodetect
[basePath, outputPath, sstPath, SGPath, matPath, Computer, defaultThreads] ...
    = useComputer(Computer);

%% Clear variables which I want to examine between runs, but not carry over.
clearvars bleachEvents bleachState mortState resultSimilarity Omega_all Omega_factor;

%% Most-used case settings
% DEFINE CLIMATE CHANGE SCENARIO (from normalized GFDL-ESM2M; J Dunne)
RCP = 'rcp45'; % options; 'rcp26', 'rcp45', 'rcp60', 'rcp85', 'control', 'control400'
E = 1;  % EVOLUTION ON (1) or OFF (0)?
OA = 1; % Ocean Acidification ON (1) or OFF (0)?
maxReefs = 1925;  %never changes, but used below.
%% Variables for plotting, debugging, or speed testing
skipPostProcessing = false;     % Don't do final stats and plots when timing.
everyx = 100; % 1;   % run code on every x reefs, plus "keyReefs" if everyx is
                    % one of 'eq', 'lo', 'hi' it selects reefs for abs(latitude)
                    % bins [0,7], (7, 14], or (14,90] respectively.  Also,
                    % if everyx >= 10000, only do reefs specifed in
                    % keyReefs
allPDFs = false;                % if false, just prints for keyReefs.
doPlots = true;                 % For optimization runs, turn off all plotting regardless
                                % of the individual flags below.
doCoralCoverMaps = false;        % World maps of cover, survival, etc. 
doCoralCoverFigure = true;     % Plot cover vs. time
doGenotypeFigure = false;
doGrowthRateFigure = false;
saveEvery = 5000;               % How often to save stillrunning.mat. Not related to everyx.
saveVarianceStats = false;      % Only when preparing to plot selV, psw2, and SST variance.

% Super symbiont options
newMortYears = false; % If true, save a fresh set of "long mortality" years based on this run.
                     % A map of those years is generated.  This also saves the years of first full-reef bleaching.
superMode = 0;  % 0 = add superAdvantage temperature to standard "hist" value.
                % 1 = use mean of temperatures in the superInitYears time range.
                % 2 = use max of temperatures in the superInitYears time range.
                %     Modes 0 to 2 all start at superStart, 2035 by default.
                % 3 = use mean of last 10 years before the end of the first 5-year mortality.
                % 4 = use max of last 10 years before the end of the first 5-year mortality.
                % 5 = use fixed delta like option 0, but start as in
                %     options 3 and 4.
                % 6 = use fixed delta like option 0, but start according to
                %     first year of bleaching.

superAdvantage = 0.0;           % Degrees C above native symbionts.
startSymFractions = [1.0 0.0];  % Starting fraction for native and super symbionts.
superStart = 2035;              % Start in this year for modes 0-2 and superAdvantage > 0.

% If this code is called from a script, we want some of the variables above
% to be overwritten.  Do that here before they are used in the code below.
if exist('scriptVars', 'var')
    disp('Called by a script');
    % Only these variables are supported now, but it's easy to add more.
    if isfield(scriptVars, 'RCP'); RCP = scriptVars.RCP; end;
    if isfield(scriptVars, 'E'); E = scriptVars.E; end;
    if isfield(scriptVars, 'OA'); OA = scriptVars.OA; end;
    if isfield(scriptVars, 'superMode'); superMode = scriptVars.superMode; end;
    if isfield(scriptVars, 'superAdvantage'); superAdvantage = scriptVars.superAdvantage; end;
    fprintf('After update in model E = %d, OA = %d, RCP = %s Mode = %d Adv = %d \n', E, OA, RCP, superMode, superAdvantage);
end

if superMode >= 3 && superMode <=5
    fn = strcat('longMortYears_', RCP, '_', num2str(E));
    load(fn, 'longMortYears');
    superStartYear = longMortYears;
    % XXX - next two lines for testing only!!!
    %subFrom = (superStartYear > 1882);
    %superStartYear = superStartYear - 20*subFrom;
elseif superMode >= 6
    fn = strcat('firstBleachYears_', RCP, '_', num2str(E));
    load(fn, 'firstBleachYears');
    superStartYear = firstBleachYears;
elseif superAdvantage == 0.0
    % If there's not advantage it's a case where there's no addition.
    % Delay to the end or we'll end up "seeding" extra symbionts.
    if strcmp(RCP, 'control400')
        last = 1860+400+1;
    else
        last = 2100 + 1;
    end
    superStartYear = last*ones(maxReefs, 1); % Until this date set S(:, 3:4) to zero.  Use 1861 to start from the beginning.
else  
    % Start the symbionts in this fixed year for all reefs.
    superStartYear = superStart*ones(maxReefs, 1); % Until this date set S(:, 3:4) to zero.  Use 1861 to start from the beginning.
end
assert(length(superStartYear) == maxReefs, 'By-year symbiont starts should have a start year for every reef.');

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
assert(startSymFractions(2) == 0, 'Start with no symbionts if they are to be suppressed at first.');
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
% Lee Stocking Island, for comparison to Fitt et al. 2000.
keyReefs = [366 144];
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
%keyReefs = [144];
%keyReefs = [402 610 823 1463];

% Reefs with the earliest mortality in the rcp85, E=1 case are listed below.  All
% experi3.53ence 5 years of mortality by 2012.
%  29  30  47  48  68  69  98  99 482 1903 1904
% These reefs reach the same condition in 2099
% 199         990        1738

% Put them in order without duplicates.  Shouldn't matter, but why not?
keyReefs = unique(keyReefs);
% too good to be true matches: keyReefs = [561 574 762 1271 1298 1299 1325 1326 1228 402];
%keyReefs = []; % can be set empty to minimize plotting and other work for some tests

%% A list of reefs for which to save data at maximum resolution for detailed
%  analysis or plotting.
dataReefs = [366];

%% Use of parallel processing on the local machine.
% no argument: uses the existing parallel pool if any.
% 0:  runs sequentially.
% > 0 tries to start the requested number of workers, failing if the value
%     is greater than the value specified in Matlab's Parallel Preferences.
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
        Plot_One_Reef();
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

initYear = '2001';  % SelV and hist will be initialized from 1861 to this year.
%% Parameters for the bleaching model (massive first, then branching)
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


%% Some useful paths and text strings to be used later.
format shortg; c = clock; dateString = strcat(num2str(c(1)),num2str(c(2),'%02u'),num2str(c(3),'%02u')); % today's date stamp
% A string frequently used in titles and file names:
modelChoices = strcat(dataset, RCP,'.E',num2str(E),'.OA',num2str(OA));
pdfDirectory = strcat(modelChoices,'_NF',num2str(NF),...
    '_sM',num2str(superMode),'_sA',num2str(superAdvantage),'_',dateString,'_figs/');
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
assert(maxReefs == length(Reefs_latlon), 'maxReefs must match the input data');

%% LOAD Omega (aragonite saturation) values if needed

if OA == 1
    
    
    
    % To hardwire OA for testing control cased: Normally just pass RCP!
    % RCPfake = 'rcp60';
    %[Omega_all] = GetOmega(SGPath, RCPfake);
    
    [Omega_all] = GetOmega(SGPath, RCP);
    % Enlarge the array to match the extended control400 array
    %{
    copyLine = Omega_all(:, 2880);
    for iii = lenTIME:-1:2881
        Omega_all(:, iii) = copyLine;
    end
    %}
    
    % Convert omegas to growth-factor multipliers so there's
    % less logic inside the time interations.
    [Omega_factor] = omegaToFactor(Omega_all);
else
    % Wasteful to make a big empty array, but it makes entering the
    % parallel loop simpler.  Note that only the last value is set.
    Omega_factor(maxReefs, lenTIME) = 0.0;
end

%% SUB-SAMPLE REEF GRID CELLS
% Since just iterating with "everyx" won't hit all keyReefs, build a list
% of reefs for the current run.
if isnumeric(everyx)
    if everyx >= 10000
        % Only use keyReefs
        toDo = [];
    else
        toDo = 1:everyx:maxReefs;   % as specified by everyx
    end
else
    % everyx can specify a reef area (eq, lo, hi)
    toDo = latitudeBin(everyx, Reefs_latlon);
end
toDo = unique([toDo keyReefs]); % add keyReefs defined above
if isempty(toDo)
    error('No reefs specified.  Exiting.');
end
reefsThisRun = length(toDo);
logTwo('Modeling %d reefs.\n', reefsThisRun);


%% LOAD SELECTIONAL VARIANCE (psw2)
load (strcat(matPath, 'Optimize_psw2.mat'),'psw2_new', 'pswInputs')
% pswInputs are not used in computations, but they are recorded to document
% each run.
% Selection of variance column from psw2_new.
if exist('optimizerMode', 'var')
    propTest = 1;
else
    propTest = getPropTest(E, RCP);
end
pswInputs = pswInputs(:, propTest);

% XXX this changes results - don't do it! psw2_new = psw2_new(:, propTest);  % no need to carry columns we never use.
%% Constants moved outside of the loop so they are only defined once:
C_seed = [0.1*10^7 0.01*10^7]; % mass branch
S_seed = [0.1*10^6 0.1*10^6 0.1*10^6 0.1*10^6];  % This was a single value until 2/17/2017

% Load .mat file for Coral and Symbiont genetics constants
% As of 1/3/2017 the file contains a structure rather than individual
% variables.
load(strcat(matPath, 'Coral_Sym_constants_2.mat')); % default is evolution OFF
% XXX override Sn for testing.
%coralSymConstants.Sn = 1;
assert(length(startSymFractions) == coralSymConstants.Sn, 'Symbiont start fractions should match number of symbionts.');

% Mutational Variance w (E=1) and w/o Evolution (E=0)
if E==0
    vM = 0;     % Mutational variance (degC^2/yr) (convert to months/12)
else
    vM = coralSymConstants.ve*.001/12;
end %%(1.142*10^-5)/12 ;    % Mutational variance (degC^2/yr) (convert to months/12)
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
% All years in this run:
fullYearRange = [startYear startYear+years-1];
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

% Convert years for symbiont activation to indexes in the time array for
% quicker use later.   superStartYear units are years.
for i = length(superStartYear):-1:1
    ssY = superStartYear(i);
    superStartIndex(i) = findDateIndex(strcat('14-Jan-', num2str(ssY)), strcat('16-Jan-',num2str(ssY)), time);
end

% max so it's always a valid index, BUT note that superStartIndex
% can be zero when a super symbiont is never needed!
superStartIndexM10 = max(1, superStartIndex - 10*stepsPerYear);

% Same for the SST and TIME arrays, but they are coarser, with dates at the
% 15th of each month.
initSSTIndex = findDateIndex(strcat('14-Dec-', initYear), strcat('16-Dec-',initYear), TIME);
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
% pass these dummy arrays and allocate the required memory in each worker.
for i = queueMax:-1:1
    kOffset(i) = min(toDoPart{i}); % Number of first reef in these chunks.
    % Inputs
    LatLon_chunk{i} = Reefs_latlon(min(toDoPart{i}):max(toDoPart{i}),1:2);
    SST_chunk{i} = SST(min(toDoPart{i}):max(toDoPart{i}), :);
    Omega_chunk{i} = Omega_factor(min(toDoPart{i}):max(toDoPart{i}), :);
    suppressSI_chunk{i} = superStartIndex(min(toDoPart{i}):max(toDoPart{i}));
    suppressSIM10_chunk{i} = superStartIndexM10(min(toDoPart{i}):max(toDoPart{i}));

    % Outputs
    resultSim_chunk{i} = nan(1,1);
    bleachEvents_chunk{i} = false(1,1);
    bleachState_chunk{i} = false(1,1);
    mortState_chunk{i} = false(1,1);

    C_cum_chunk{i} = zeros(length(time), coralSymConstants.Sn*coralSymConstants.Cn); % Sum coral cover for all reefs.
    % 3D array sized (time by reef by coral type).  Note that we don't care
    % about the identity of the reefs in this case, so we just need enough
    % columns for all reefs actually calculated, ignoring those which are
    % skipped.
    C_year_chunk{i} = zeros(years, chunkSize, coralSymConstants.Cn); % Coral cover for all reefs, but just 2 columns.
    Massive_dom_chunk{i} = zeros(length(time), 1);

end
% TODO input arrays such as SST and Reefs_latlon are sent at full size to
% each worker.  Consider sending just the correct subset to each.

%% RUN EVOLUTIONARY MODEL
iteratorHandle = selectIteratorFunction(length(time), Computer);
% the last argument in the parfor specifies the maximum number of workers.
parfor (parSet = 1:queueMax, parSwitch)
%for parSet = 1:queueMax
    %  pause(1); % Without this pause, the fprintf doesn't display immediately.
    %  fprintf('In parfor set %d\n', parSet);
    reefCount = 0;
    % How often to print progress.
    printFreq = max(10, ceil(length(toDoPart{parSet})/4)); % The last digit is the number of pieces to report.

    % Variables to collect and return, since parfor won't allow direct
    % insertion into an output array.
    % TODO see if length(toDoPart(parSet)) should be used instead of
    % chunkSize.
    par_bleachEvents = false(length(toDoPart(parSet)), years, coralSymConstants.Cn); %#ok<PFBNS>
    par_bleachState = par_bleachEvents;
    par_mortState = par_bleachEvents;
    par_SST = SST_chunk{parSet};
    par_Omega = Omega_chunk{parSet};
    par_LatLon = LatLon_chunk{parSet};
    par_SupressSI = suppressSI_chunk{parSet};
    par_SupressSIM10 = suppressSIM10_chunk{parSet};
    par_kOffset = kOffset(parSet);
    par_C_cum = C_cum_chunk{parSet};
    par_C_year = C_year_chunk{parSet};
    par_Massive_dom = Massive_dom_chunk{parSet};
    par_HistSuperSum = 0.0;
    par_HistOrigSum = 0.0;
    par_HistOrigEvolvedSum = 0.0;
    for k = toDoPart{parSet}
        reefCount = reefCount + 1;
        kChunk = 1 + k - par_kOffset;
        SST_LOC = par_SST(kChunk, :);                     % Reef grid cell location
        SSThist = SST_LOC';                      % Transpose SST matrix
        Omega_hist = par_Omega(kChunk, :);
        psw2 = NF*psw2_new(k, propTest) ;         % UPDATED** max(0.3,min(1.3,(mean(exp(0.063.*SSThist))./var(exp(0.063.*SSThist))).^0.5/7)); % John's new eqn 8/10/16** try this
        reefLatlon = par_LatLon(kChunk, :);
        suppressSI = par_SupressSI(kChunk);
        suppressSIM10 = par_SupressSIM10(kChunk);
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
        SelVx = repmat(SelV,1,coralSymConstants.Sn);     % Selectional variance matrix for coral calcs
        %SelV = [1.25 1]*psw2*var(SSThist_anom(:))

        % Initialize symbiont genotype, sym/coral population sizes, carrying capacity

        %ssss = findDateIndex(strcat('14-Jan-', num2str(par_SuppressYears(kChunk)-10)), strcat('16-Jan-',num2str(par_SuppressYears(kChunk)-10)), time);
        %eeee = findDateIndex(strcat('14-Dec-', num2str(par_SuppressYears(kChunk))), strcat('16-Dec-',num2str(par_SuppressYears(kChunk))), time);

        [vgi, gi, S, C, hist, ri] = Init_genotype_popsize(time, initIndex, temp, coralSymConstants, ...
            E, vM, SelV, superMode, superAdvantage, startSymFractions, ...
            [suppressSI suppressSIM10]);

       

        %fprintf('super will start at index %d\n', suppressSI);
        %% MAIN LOOP: Integrate Equations 1-5 through 2100 using Runge-Kutta method
        % timeIteration is called here, with the version determined by
        % iteratorHandle.
        [S, C, ri, gi, vgi, origEvolved] = iteratorHandle(timeSteps, S, C, dt, ...
                    ri, temp, OA, omega, vgi, gi, MutVx, SelVx, C_seed, S_seed, suppressSI, ...
                    superSeedFraction, oneShot, coralSymConstants); %#ok<PFBNS>
        % These, with origEvolved, compare the average native and
        % supersymbiont genotypes with the evolved state of the native
        % symbionts just before the supersymbionts are introduced.
        origHist = gi(1,1);
        superHist = gi(1,3);

        par_HistSuperSum = par_HistSuperSum + superHist;
        par_HistOrigSum = par_HistOrigSum + origHist;
        par_HistOrigEvolvedSum = par_HistOrigEvolvedSum + origEvolved;

        if any(dataReefs == k) % Save detailed history
            matName = strcat('DetailedSC_Reef', num2str(k), '_', modelChoices, '.mat');
            saveAsMat(matName, C, S, time, temp);
        end
        if doPlots && (doGrowthRateFigure || doGenotypeFigure) && any(keyReefs == k)  % temporary genotype diagnostic
            suff = '';
            if superMode && superMode ~= 5
                suff = sprintf('_%s_E%d_SymStrategy%d_Reef%d', RCP, E, superMode, k);
            elseif superMode == 0 || superMode == 5
                suff = sprintf('_%s_E%d_SymStrategy%dAdv%0.2fC_Reef%d', RCP, E, superMode, superAdvantage, k);
            end
            if doGenotypeFigure
                genotypeFigure(fullMapDir, suff, k, time, gi, suppressSI); %#ok<UNRCH>
            end
            if doGrowthRateFigure
                % Growth rate vs. T as well
                % TODO: dies when suppressSI = 0
                if strcmp(RCP(1:3), 'rcp') %#ok<UNRCH>
                    growthRateFigure(fullMapDir, suff, datestr(time(suppressSI), 'yyyy'), ...
                        k, temp, fullYearRange, gi, vgi, suppressSI, ...
                        coralSymConstants, SelVx, RCP);         
                end
            end
        end

        par_C_cum = par_C_cum + C;
        par_Massive_dom = par_Massive_dom + C(:, 1) > C(:, 2);
        % Time and memory will be consumed, but we need stats on coral
        % cover.
        par_C_year(:, reefCount, 1) =  decimate(C(:, 1), stepsPerYear, 'fir');
        par_C_year(:, reefCount, 2) =  decimate(C(:, 2), stepsPerYear, 'fir');
        
        %% New clean stats section
        
        [ C_monthly, S_monthly, ~, bleachEventOneReef, bleachStateOne, mortStateOne ] = ...
            Clean_Bleach_Stats(C, S, C_seed, S_seed, dt, TIME, bleachParams, coralSymConstants);
     
        if doPlots && (any(keyReefs == k) || allPDFs)
            % Now that we have new stats, reproduce the per-reef plots.
            Plot_One_Reef(C_monthly, S_monthly, bleachEventOneReef, psw2, time, temp, lat, lon, RCP, ...
                  hist, Data, SGPath, outputPath, k, ...
                  pdfDirectory, LOC, NF, E, dateString, lenTIME);
        end

        if ~isempty(bleachEventOneReef)
            % bleachEventOneReef is returned as a sparse array which is
            % great for plotting and saving space.  Unfortunately, sparse
            % arrays are limited to 2D, so here it gets expanded.  It could
            % be better to return a set or list of sparse arrays, but this
            % is much clearer.
            par_bleachEvents(reefCount, :, :) = full(bleachEventOneReef);
        end
        % Like bleachEventOneReef, but the bleaching and mortality states
        % are never empty.  Also, these aren't stored as sparse since they
        % can have long sequences of "true".
        par_bleachState(reefCount, :, :) = bleachStateOne;
        par_mortState(reefCount, :, :) = mortStateOne;

        if parSwitch && mod(reefCount, printFreq) == 0
            fprintf('Set %d is %3.0f percent complete.\n', parSet, (100*reefCount/length(toDoPart{parSet})));
        end

    end % End of reef areas for one parallel chunk

    bleachEvents_chunk{parSet} = par_bleachEvents;
    bleachState_chunk{parSet} = par_bleachState;
    mortState_chunk{parSet} = par_mortState;
    C_cum_chunk{parSet} = par_C_cum;
    C_year_chunk{parSet} = par_C_year(:, 1:reefCount, :);
    Massive_dom_chunk{parSet} = par_Massive_dom;
    histSuper_chunk(parSet) = par_HistSuperSum;
    histOrig_chunk(parSet) = par_HistOrigSum;
    histOrigEvolved_chunk(parSet) = par_HistOrigEvolvedSum;


end % End of parfor loop
% Clear variables use only as inputs inside the loop.
clearvars SST_chunk Omega_chunk LatLon_chunk;

% Build these variables from the chunks.
clearvars SST_chunk;

% bleachEvents_chunk contains a chunk per worker.  Each worker's chunk
% contains a 3D array where the first index is the sequential number of the
% reef in that to-do chunk.  Here we build a full 3D array for all possible
% reefs, sized (reefs, years, coral types).
% old way: bleachEvents = horzcat(bleachEvents_chunk{:});
bleachEvents = false(maxReefs, years, coralSymConstants.Cn);
bleachState = false(maxReefs, years, coralSymConstants.Cn);
mortState = false(maxReefs, years, coralSymConstants.Cn);
for i = 1:queueMax
    tdp = toDoPart{i}';
    chunkE = bleachEvents_chunk{i};
    chunkB = bleachState_chunk{i};
    chunkM = mortState_chunk{i};
    assert(size(chunkE, 1) == length(tdp), 'Number of bleach event results must match list of reef numbers.');
    assert(size(chunkB, 1) == length(tdp), 'Number of bleach state results must match list of reef numbers.');
    assert(size(chunkM, 1) == length(tdp), 'Number of mortality state results must match list of reef numbers.');
    for chunkIndex = 1:size(chunkE, 1)
        k = toDoPart{i}(chunkIndex);
        bleachEvents(k, :, :) = chunkE(chunkIndex, :, :);
        bleachState(k, :, :) = chunkB(chunkIndex, :, :);
        mortState(k, :, :) = chunkM(chunkIndex, :, :);
    end
end
clearvars bleachEvents_chunk bleachState_chunk mortState_chunk; % release some memory.

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
clearvars C_cum_chunk C_year_chunk Massive_dom_chunk histSuper_chunk histOrig_chunk histOrigEvolved_chunk;
superSum = superSum/reefsThisRun;
histSum = histSum/reefsThisRun;
histEvSum = histEvSum/reefsThisRun;
logTwo('Super symbiont genotype = %5.2f C.  Base genotype %5.2f C (advantage %5.2f), Evolved base %5.2f (advantage %5.2f).\n', ...
    superSum, histSum, (superSum-histSum), histEvSum, (superSum-histEvSum));



if ~skipPostProcessing

    % Count bleaching events between 1985 and 2010 inclusive.
    i1985 = 1985 - startYear + 1;
    i2010 = 2010 - startYear + 1;
    % Count by reef
    events85_2010(maxReefs) = 0;
    eventsAllYears(maxReefs) = 0;
    for k = 1:maxReefs
        events85_2010(k) = nnz(bleachEvents(k, i1985:i2010, :));
        eventsAllYears(k) = nnz(bleachEvents(k, :, :));
    end
    % Count for all reefs over this time period.
    count852010 = sum(events85_2010);
    
    Bleaching_85_10_By_Event = 100*count852010/reefsThisRun/(2010-1985+1);
    fprintf('Bleaching by event = %6.4f\n', ...
        Bleaching_85_10_By_Event);

    % Last modeled year, so we don't count a coral which lives to then as
    % dead.
    lastYear = str2double(datestr(TIME(end), 'yyyy'));

    % Build an array with the last year each reef is alive.
    % First add a column (indexed r) to mortState which is true when all
    % coral types are dead.
    % Also find the last bleaching event here.
    fullReef = coralSymConstants.Cn + 1;
    lastYearAlive = nan(maxReefs, 1);
    lastBleachEvent = nan(maxReefs, fullReef);
    for k = 1:maxReefs
        for i = 1:years
            mortState(k, i, fullReef) = all(mortState(k, i, 1:fullReef-1));
            bleachState(k, i, fullReef) = all(bleachState(k, i, 1:fullReef-1));
            % Now find the last year alive - leave NaN if it ends alive.
        end
        if mortState(k, years, fullReef)
            ind = find(~mortState(k, :, fullReef), 1, 'last');
            assert(~isempty(ind), 'Reef %d should never start out dead.', k);
            lastYearAlive(k) = ind(1) + startYear - 1;
        end
        for rr = 1:fullReef
            if bleachState(k, years, rr)
                ind = find(~bleachState(k, :, rr), 1, 'last');
                assert(~isempty(ind), 'Reef %d coral type %d should never start out bleached.', k, rr);
                lastBleachEvent(k, rr) = ind(1) + startYear - 1;
            end
        end
    end
    frequentBleaching = defineFrequentBleaching(bleachEvents);
    
    if saveVarianceStats
        assert(length(toDo) == maxReefs, 'Only save variance data when running all reefs!'); %#ok<UNRCH>
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
        save(strcat(basePath, 'LastYear', '_selV_', RCP, 'E=', num2str(E), ...
            'OA=', num2str(OA), '.mat'), ...
            'psw2_new', 'selVariance', 'tVariance', 'lastYearAlive', 'RCP', 'OA', 'E');
    end

    format shortg;
    filePrefix = strcat(modelChoices,'_NF',num2str(NF),'_',dateString);
    % Don't save all this data if we're just optimizing.
    if doPlots
        % fname = strcat(filePrefix,'.mat');
        fname = strcat(outputPath, pdfDirectory, filePrefix, '.mat');
        save(fname, 'toDo', ...
            'E','OA','pdfDirectory','dataset', ...
            'Reefs_latlon','everyx','NF','RCP','reefsThisRun');

        if doCoralCoverMaps
            MapsCoralCoverClean(fullMapDir, Reefs_latlon, toDo, lastYearAlive, ...
                events85_2010, eventsAllYears, frequentBleaching, ...
                mortState, bleachState, ...
                fullYearRange, ...
                modelChoices, filePrefix); %#ok<UNRCH>
        end

        if doCoralCoverFigure
            coralCoverFigure(C_yearly, coralSymConstants, startYear, years, RCP, E, OA, superMode, ...
                    superAdvantage, fullMapDir) %#ok<UNRCH>
        end
    end
    % Note that percentMortality is not used in normal runs, but it is
    % examined by the optimizer when it is used.
    [percentBleached, percentMortality] = ...
        Stats_Tables(bleachState, mortState, lastYearAlive, ...
        lastBleachEvent, frequentBleaching, toDo, Reefs_latlon, outputPath, startYear, RCP, E, OA, ...
        bleachParams);

    % Get the years when reefs first experienced lasting mortality and 
    % bleaching.  This isn't wanted in every run, and certainly not when 
    % super symbionts are introduced in a variable way.
    if superMode == 0 && newMortYears
        if everyx ~= 1
            disp('WARNING: saving mortality and bleaching should only be done when all reefs are computed.');
        end
        saveMortYears(mortState, startYear, RCP, E, OA, fullMapDir, ...
            modelChoices, filePrefix, Reefs_latlon, bleachState, maxReefs);
    end

    logTwo('Bleaching by event = %6.4f\n', Bleaching_85_10_By_Event);
end % End postprocessing block.

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
if ~exist('optimizerMode', 'var') || optimizerMode == false && ...
    (saveToExcel && ~skipPostProcessing)

    saveExcelHistory(basePath, now, RCP, E, everyx, queueMax, elapsed, ...
        Bleaching_85_10_By_Event, bleachParams, pswInputs);
end
