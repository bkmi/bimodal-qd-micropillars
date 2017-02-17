function [ soln, figFirst, figLast ] = sweeper( ind_parSweep, parBound, param, varargin )
%This function is designed to sweep a parameter by doing a turn on time
%series at the first parameter value then sweep the time series to the next
%parameter value.
%
%   Input:
%       ind_parSweep, ...
%       parBounds, ...
%       param, ...
%       varargin
%
%   Output:
%       soln
%
%   Options:
%       'durTimeSeries' = 10
%           Duration of time series for swept time series. Turn on is
%           longer by 50%
%       'numSweepSteps' = 50
%           Number of steps to go from parBound(1) to parBound(2).
%       'oscStepsCheck' = 10
%           Number of timesteps backwards to check for oscillations.
%       'oscTol' = 1e-6
%           Maximum difference between the variable values at each step
%           back when checking for oscillations. Larger values are less
%           sensitive to oscillation.
%
%   master_options:
%       'save' = 0, 1
%           By default, this is set to 0. When 'save' = 0, the function
%           does not try to save anything. When 'save' = 1, the function 
%           tries to save ________.
%       'datadir_specific' = '../data_bimodal-qd-micropillars/'
%           By default, this is set as above.
%       'dimensional' = 0, 1
%           By default, this is set to 0. When 'dimensional' = 0, the
%           function applies a non-dimensionalized system. When
%           'dimensional' = 1, the function applies a dimensionalized
%           system.

%% Defaults + inputParser + Organize behavior

p = inputParser;

% General option defaults
p.addParameter('durTimeSeries', 10)
p.addParameter('oscStepsCheck', 10)
p.addParameter('oscTol', 1e-6)
p.addParameter('numSweepSteps', 50)


% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% first parse to set 'par' variable
parse(p,varargin{:})

options = p.Results;


%% Populate soln struct
% Extract the structure for the parameter which is being swept and also
% create an array with each of the swept values.
parSweep = param.(param.var_names{ind_parSweep});
parArray = linspace(parBound(1),parBound(2),options.numSweepSteps);

% preallocate for soln
soln = repmat(...
    struct(...
    parSweep.var_name,0,... % Which parameter is being swept
    'param',struct, ...
    'timeSeries',struct, ...
    'calcTime',0), ...
    [options.numSweepSteps,1]);

for i = 1:options.numSweepSteps
    % fill the soln structure up with the relevant parameters at that
    % particular sweet step.
    soln(i).(parSweep.var_name) = parArray(i);
    soln(i).param = setup_params_nonDim_CnstCplRatio( ...
        parSweep.var_name, parArray(i), ...
        'populate_wrkspc',0, ...
        'save',0, ...
        'clear',0);
end


%% Calculate the turn on time series
tic; % Start timer

[soln(1).timeSeries,figFirst] = solver( ...
    [1e-9;0;1e-9;0;0;0], ...
    [0,1.5*options.durTimeSeries], ...
    soln(1).param, ...
    'plot',1, ...
    'quiet', 1);
% 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)

soln(1).calcTime = toc; % Stop Timer


%% Check for oscillations then calculate the swept time series
for i = 2:options.numSweepSteps
    % Check first or oscillations in the last time series.
    for j = 1:options.oscStepsCheck
        % Alert the user if the past options.oscStepsCheck number of
        % timesteps are changing values up to a tolerance set by
        % options.oscTol
        if (soln(i-1).timeSeries.y(:,end) - soln(i-1).timeSeries.y(:,end-j)) ...
                >= options.oscTol
            error(['Your soln(', num2str(i-1),') is not steady state.'])
        end
    end
    
    % Now sweep the parameter and calculate another time series.
    tic; % Start timer
    if i ~= options.numSweepSteps
        % Everything but the last solver
        soln(i).timeSeries = solver( ...
            soln(i-1).timeSeries, ...
            [soln(i-1).timeSeries.x(end), ... % last value from last time series.
            soln(i-1).timeSeries.x(end)+options.durTimeSeries], ...
            soln(1).param,...
            'quiet', 1);
    elseif i == options.numSweepSteps
        % For the last solver only.
        [soln(i).timeSeries,figLast] = solver( ...
            soln(i-1).timeSeries, ...
            [soln(i-1).timeSeries.x(end), ... % last value from last time series.
            soln(i-1).timeSeries.x(end)+options.durTimeSeries], ...
            soln(1).param,...
            'quiet', 1, ...
            'plot', 1);
    end
    soln(i).calcTime = toc; % Stop Timer
end

