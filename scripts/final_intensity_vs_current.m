%% Report final intensity from a turn-on time series at different currents

clear;

% Threshold is ~43e-6 Amps
Jmin = 43e-6;
Jmax = 200e-6;
numPoints = 50;

solnTurnOn = repmat(...
    struct(...
    'timeSeries',struct, ...
    'J',0,...
    'params',struct, ...
    'calcTime',0),[numPoints,1]);
Jarray = linspace(...
    Jmin, ...
    Jmax, ...
    numPoints); 


%% Turn On Solver
timeTotal = 0;

for i = 1:numPoints
    tic; % Start timer
    
    setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(i), ...
        'feed_ampli',0, ...
        'feed_ampliMatrix', [0, 0; 0, 0], ...
        'feed_phase',0, ...
        'clear',0);
    dde23_soln = solver([1e-9;0;1e-9;0;0;0], [0,15], param, master_options);
    % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)
    
    solnTurnOn(i).timeSeries = dde23_soln;
    solnTurnOn(i).J = param.J.value;
    solnTurnOn(i).param = param;
    solnTurnOn(i).calcTime = toc;
    timeTotal = timeTotal + toc;
end


%% Plot final intensity versus current amplitude

strongFinalIntensity = zeros(numPoints,1);
weakFinalIntensity = zeros(numPoints,1);
for i = 1:numPoints
    strongFinalIntensity(i) = norm( ...
        [solnTurnOn(i).timeSeries.y(1,end),...
        solnTurnOn(i).timeSeries.y(2,end)] );
    weakFinalIntensity(i) = norm( ...
        [solnTurnOn(i).timeSeries.y(3,end),...
        solnTurnOn(i).timeSeries.y(4,end)] );
end

fig1 = figure;
% strong red
semilogy(Jarray,strongFinalIntensity,'r')
hold on
% weak blue
semilogy(Jarray,weakFinalIntensity,'b')
title('Turn on final intensity vs current amplitude')
hold off


%% Sweep Solver
timeTotalSweep = 0;

solnSweep = repmat(...
    struct(...
    'timeSeries',struct, ...
    'J',0,...
    'params',struct, ...
    'calcTime',0),numPoints);

for i = 1
    tic; % Start timer

    setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(i), ...
        'feed_ampli',0, ...
        'feed_ampliMatrix', [0, 0; 0, 0], ...
        'feed_phase',0, ...
        'clear',0);
    dde23_soln = solver([1e-9;0;1e-9;0;0;0], [0,10], param, master_options);
    % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)

    solnSweep(i).timeSeries = dde23_soln;
    solnSweep(i).J = param.J.value;
    solnSweep(i).param = param;
    solnSweep(i).calcTime = toc;
    timeTotal = timeTotal + toc;
end

for i = 2:numPoints
    tic; % Start timer
    
    setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(i), ...
        'feed_ampli',0, ...
        'feed_ampliMatrix', [0, 0; 0, 0], ...
        'feed_phase',0, ...
        'clear',0);
    dde23_soln = solver(solnSweep(i-1).timeSeries, ...
        [10*(i-1),10*(i-1)+10], param, master_options);
    % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)
    
    solnSweep(i).timeSeries = dde23_soln;
    solnSweep(i).J = param.J.value;
    solnSweep(i).param = param;
    solnSweep(i).calcTime = toc;
    timeTotal = timeTotal + toc;
end

disp('It took')
disp(timeTotal)

%% Plot final intensity versus current amplitude

strongFinalIntensitySweep = zeros(numPoints,1);
weakFinalIntensitySweep = zeros(numPoints,1);
for i = 1:numPoints
    strongFinalIntensitySweep(i) = norm( ...
        [solnSweep(i).timeSeries.y(1,end),...
        solnSweep(i).timeSeries.y(2,end)] );
    weakFinalIntensitySweep(i) = norm( ...
        [solnSweep(i).timeSeries.y(3,end),...
        solnSweep(i).timeSeries.y(4,end)] );
end


fig2 = figure;
% strong red
semilogy(Jarray,strongFinalIntensitySweep,'r')
hold on
% weak blue
semilogy(Jarray,weakFinalIntensitySweep,'b')
title('Final intensity vs Current amplitude, "Sweep"')
hold off

