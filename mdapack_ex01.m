
% a simple example on how to run MDAPACK

clear variables;
clear import;
close all;

% create a PRNG to use with the model - this makes the experiment
% "deterministic", if passed as option.
s1=RandStream('mt19937ar','seed',1);


%% 1) Choose representation, model and model parameters
parameters.representation = { ...
    % we use an ensemble, since we want to use EnKF
     @(model) representations.Ensemble(model, ...
     'ensSize', 20 ...
     ) ... 
};

parameters.model = {...
    @() models.Lorenz1984(...
    'measurementSchedule' , 10:2:100, ...
    'plot2dDim', 1, ...
    'plotMode', '2d', ...
    'tEnd', 120, ...
    'evidenceNoiseRng', s1 ...
    )...
% Explanation:
% measurementSchedule: 10 days spin-up, then every second day a
% measurement, until day 100. tEnd = 120 means that we run the last 20 days
% without measurements -> "forecast period".
% plot2DDim  = 1 means the first dimension of the attractor is plotted.
% Possible values are 1,2,3. plotMode could be set to 3d, to have a nice,
% but pretty much useless, 3d plot of the assimilation.
};

%% 2) Choose method and method parameters
parameters.method = {...
    % we use a standard EnKF (no options set -> no "special features")
@()filters.Pajonk.AREnKF()
};

%% 3) Set global parameters
parameters.visualize = {true};      % output live graphics?
parameters.progress = {true};       % output progress information to console?

parameters.nRuns = {1};             % number of runs (to create statistics, requires the use of controllers.multiRun)
parameters.statistics = {true};     % compute statistics over the NRUNS or store all runs?
parameters.stat_spectrum = {false}; % compute SVD of ensemble for results? (normally not necessary)
parameters.stat_mean = {true};      % store mean? (this is typically the "best guess")
parameters.stat_relErrors = {true}; % store RMSE/relErr?
parameters.stat_var = {true};       % store variance?
parameters.stat_skewness = {false};
parameters.stat_kurtosis = {false};
parameters.stat_ksdensity = {false};
parameters.stat_summary = {false};  % store quantiles 0.025, 0.25, 0.5, 0.75, 0.975?
parameters.stat_truth = {true};     % store truth?


%% 4) Run
results = tools.batchcall(parameters, @controllers.singleRun);

% all output can be found in the "results" structure, but is not saved to a
% file