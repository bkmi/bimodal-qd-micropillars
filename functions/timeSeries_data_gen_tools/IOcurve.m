function [ sweepUp, sweepDown, turnOn ] = IOcurve( feedAmpMat, fbAmp, ...
    varargin)
%This function produces an input output curve with a sweep up and sweep
%down.
%   Detailed explanation goes here
%% Defaults + inputParser + Organize behavior

p = inputParser;

% General option defaults
p.addParameter('Jmin',43e-6)
p.addParameter('Jmax',750e-6)
p.addParameter('numPoints',150)
p.addParameter('feed_phase',0)


% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% first parse to set 'par' variable
parse(p,varargin{:})

options = p.Results;


%% Setup values
Jarray = linspace(...
    options.Jmin, ...
    options.Jmax, ...
    options.numPoints);

basicStruct = struct(...
    'fbAmp', 0, ...
    'J',0,...
    'param',struct, ...
    'timeSeries',struct, ...
    'calcTime',0);

solnTurnOn = repmat(...
    basicStruct,[options.numPoints,1]);

%% Turn On Solver
% produces the 'turnOn' structure, used for calculating Sweep Down as well.

% Current is constant in a row
% fbAmp   is constant in a column
turnOn = repmat( ...
    basicStruct,...
    options.numPoints, ...   % rows, current is constant 
    numel(fbAmp));           % columns, fbAmp is constant

for j = 1:numel(fbAmp)
    for i = 1:options.numPoints
        tic; % Start timer
        
        param = setup_params_nonDim_CnstCplRatio(...
            'save',0, ...
            'J', Jarray(i), ...
            'feed_ampli',fbAmp(j), ...
            'feed_ampliMatrix', feedAmpMat, ...
            'feed_phase',options.feed_phase, ...
            'clear',0,...
            'populate_wrkspc', 0);
        dde23_soln = solver([1e-9;0;1e-9;0;0;0], [0,9], param, 'quiet',1);
        % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)

        solnTurnOn(i).fbAmp = fbAmp(j);
        solnTurnOn(i).J = param.J.value;
        solnTurnOn(i).param = param;
        solnTurnOn(i).timeSeries = dde23_soln;
        solnTurnOn(i).calcTime = toc;
    end
    
    % Add to turnOn
    turnOn(:,j) = solnTurnOn;
end


%% Sweep Up Solver

% Current is constant in a row
% fbAmp   is constant in a column
sweepUp = repmat( ...
    basicStruct,...
    options.numPoints, ...   % rows, current is constant 
    numel(fbAmp));           % columns, fbAmp is constant

for j = 1:numel(fbAmp)

    param = setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(1), ...
        'feed_ampli',fbAmp(j), ...
        'feed_ampliMatrix', feedAmpMat, ...
        'feed_phase',options.feed_phase, ...
        'populate_wrkspc', 0, ...
        'clear',0);

    solnSweepUp = sweeper( ...
        param.J.index, ...
        [options.Jmin, options.Jmax], ...
        param, ...
        'plot',0, ...
        'numSweepSteps', options.numPoints);

    % Add to sweep
    sweepUp(:,j) = solnSweepUp;
end


%% Sweep Down Solver

% Current is constant in a row
% fbAmp   is constant in a column
sweepDown = repmat( ...
    basicStruct,...
    options.numPoints, ...   % rows, current is constant 
    numel(fbAmp));           % columns, fbAmp is constant

for j = 1:numel(fbAmp)
    initStruct = turnOn(:,j);
    
    solnSweepDown = sweeper(...
        initStruct(end).param.J.index, ...
        [options.Jmax, options.Jmin], ...
        initStruct(end).param, ...
        'plot',0, ...
        'numSweepSteps', options.numPoints, ...
        'hist', initStruct(end).timeSeries);
    
    % Add to downSweep
    sweepDown(:,j) = solnSweepDown;
end

end

