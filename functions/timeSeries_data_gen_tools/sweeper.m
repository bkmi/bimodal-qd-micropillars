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
%       'plot' = 0, 1
%           plot = 0 means nothing will be plotted
%           plot = 1 means the first and last time series will be plotted.
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
%       'save_name' = 'dde23_soln_name'
%           The solver will save the dde23_soln as 'dde23_soln_name' in 
%           a datadir_specific given by master_options. It will overwrite.
%       'quiet' = 0, 1
%           0 -> usual output.
%           1 -> no disp output.
%       'hist' = [1e-9;0;1e-9;0;0;0] OR dde23_soln
%           Allows the user to specify the history vector. If the given
%           hist is a dde23_soln then sweeper will continue from the
%           previous solution.
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
p.addParameter('plot',0)
p.addParameter('durTimeSeries', 10)
p.addParameter('oscStepsCheck', 10)
p.addParameter('oscTol', 1e-6)
p.addParameter('numSweepSteps', 50)
p.addParameter('save_name', 'dde23_soln')
p.addParameter('quiet', 1)
p.addParameter('hist',[1e-9;0;1e-9;0;0;0])


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
    'fbAmp', param.feed_ampli.value, ...
    parSweep.var_name,0,... % Which parameter is being swept
    'param',struct, ...
    'timeSeries',struct, ...
    'calcTime',0), ...
    [options.numSweepSteps,1]);

for i = 1:options.numSweepSteps
    % fill the soln structure up with the relevant parameters at that
    % particular sweet step.
    soln(i).(parSweep.var_name) = parArray(i);
    soln(i).param = updateParams(param, parSweep.var_name, parArray(i));
    
end


%% Calculate the turn on time series
tic; % Start timer


% Using the default value for hist and an input value for hist works the
% same way thanks to the design of 'solver'
if isa(options.hist, 'double')
    % When hist is a vector
    [soln(1).timeSeries,figFirst] = solver( ...
        options.hist, ...
        [0,1.5*options.durTimeSeries], ...
        soln(1).param, ...
        'plot',options.plot, ...
        'quiet', options.quiet);
    % ,'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot));
elseif isa(options.hist, 'struct')
    % When it's a struct, AKA continuing the previous sweep
    [soln(1).timeSeries,figFirst] = solver( ...
        options.hist, ...
        [options.hist.x(end), ...
         options.hist.x(end) + 1.5*options.durTimeSeries], ...
        soln(1).param, ...
        'plot',options.plot, ...
        'quiet', options.quiet);
    % ,'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot));
end

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
            soln(i).param,...
            'quiet', 1);
    elseif i == options.numSweepSteps
        % For the last solver only.
        [soln(i).timeSeries,figLast] = solver( ...
            soln(i-1).timeSeries, ...
            [soln(i-1).timeSeries.x(end), ... % last value from last time series.
            soln(i-1).timeSeries.x(end)+options.durTimeSeries], ...
            soln(i).param,...
            'quiet', 1, ...
            'plot', options.plot);
    end
    soln(i).calcTime = toc; % Stop Timer
end


%% Save
% Save, if necessary
datadir_specific = options.datadir_specific;

% Where will it save?
if options.quiet == 0;
    if options.save == 1
        fprintf(strcat('\n\n Saving in subfolder:\n', datadir_specific,'\n'))
    end
end
    
dde23_soln = soln(end).timeSeries;
sweepSoln = soln;

if options.save == 1 && ...
        ~exist(strcat(datadir_specific,options.save_name,'.mat'),'file')
    save(strcat(datadir_specific,options.save_name), ...
        'dde23_soln', 'sweepSoln')
elseif options.save == 1 && ...
        exist(strcat(datadir_specific,options.save_name,'.mat'),'file')
    warning('That file %s already exists. Overwriting.', ...
        strcat(datadir_specific,options.save_name) )
    save(strcat(datadir_specific,options.save_name), ...
        'dde23_soln', 'sweepSoln')
end

end