
% an example on how to run MDAPACK as batch call

clear variables;
clear import;
close all;

% create a PRNG to use with the model - this makes the experiment
% "deterministic", if passed as option.
s1=RandStream('mt19937ar','seed',1);

%% 1) Choose representation, model and model parameters
parameters.representation = { ...
    % three PCEs, each with a different order
    @(model) representations.PCE( ...
      model, ...
     'pceOrder', 1, ...
     'sampleSize', 10000) ... % also used to compute quantiles
    @(model) representations.PCE( ...
      model, ...
     'pceOrder', 2, ...
     'sampleSize', 10000) ... % also used to compute quantiles
    @(model) representations.PCE( ...
      model, ...
     'pceOrder', 3, ...
     'sampleSize', 10000) ... % also used to compute quantiles
};

parameters.model = {... % bind all parameters of the model generator call
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
@()filters.Pajonk.LPCU()
};


%% 3) Set global parameters
parameters.visualize = {true};      % output live graphics?
parameters.progress = {true};       % output progress information to console?

parameters.nRuns = {1000};          % number of runs (to create statistics, requires the use of controllers.multiRun)
parameters.statistics = {true};     % compute statistics over the NRUNS or store all runs?
parameters.stat_spectrum = {false}; % compute SVD of ensemble for results? (normally not necessary)
parameters.stat_mean = {true};      % store mean?
parameters.stat_skewness = {false};
parameters.stat_kurtosis = {false};
parameters.stat_ksdensity = {false};
parameters.stat_relErrors = {true}; % store RMSE/relErr?
parameters.stat_var = {true};       % store variance? (crude measure for the "uncertainty")
parameters.stat_summary = {false};  % store quantiles 0.025, 0.25, 0.5, 0.75, 0.975?
parameters.stat_truth = {true};     % store truth?


%% 4) Run
tools.batchcall(parameters, @controllers.multiRun, ['results/current/', mfilename, '_%02d.mat']);

% The output will be written to a file with the name
% "results/current/mdapack_ex03_XX.mat", where XX is replaced by a
% consecutive number. In this case we will have three experiments, since we
% have three representations, one model and one method. The statistics for
% each experiment will be computed from 1000 repetitions (see nRuns above).
