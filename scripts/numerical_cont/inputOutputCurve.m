%% Report final intensity from a turn-on time series at different currents

clear;

% Threshold is ~43e-6 Amps
Jmin = 43e-6;
Jmax = 750e-6;
numPoints = 150;

turnOn = struct;

solnTurnOn = repmat(...
    struct(...
    'timeSeries',struct, ...
    'J',0,...
    'param',struct, ...
    'calcTime',0),[numPoints,1]);
Jarray = linspace(...
    Jmin, ...
    Jmax, ...
    numPoints); 

% feedback amp settings
feedAmpMat = [1, 0; 0, 0];
fbAmp = [0, 0.1, 0.2, 0.3, 0.4, 0.5];

% save location
datadir = '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/inputOutput/';

% Save extras: fbAmp, numPoints
save([datadir, ...
    'extras.mat'], ...
    'fbAmp', 'numPoints', 'Jarray','datadir')

%% Turn On Solver
timeTotal = 0;

for j = fbAmp

    for i = 1:numPoints
        tic; % Start timer

        param = setup_params_nonDim_CnstCplRatio(...
            'save',0, ...
            'J', Jarray(i), ...
            'feed_ampli',j, ...
            'feed_ampliMatrix', feedAmpMat, ...
            'feed_phase',0, ...
            'clear',0,...
            'populate_wrkspc', 0);
        dde23_soln = solver([1e-9;0;1e-9;0;0;0], [0,9], param, 'quiet',1);
        % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)

        solnTurnOn(i).timeSeries = dde23_soln;
        solnTurnOn(i).J = param.J.value;
        solnTurnOn(i).param = param;
        solnTurnOn(i).calcTime = toc;
        timeTotal = timeTotal + toc;
    end
    
    disp('It took')
    disp(timeTotal)
    
    % Add to turnOn
    turnOn.(['fb',num2str(j*10)]) = solnTurnOn;

    % Plot final intensity versus current amplitude

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
    semilogy(Jarray,strongFinalIntensity,'Color','r')
    hold on
    % weak blue
    semilogy(Jarray,weakFinalIntensity,'b')
    title(['Turn On Time Series: Steady state intensity vs Current amplitude, k_{ss}=',num2str(j)])
    set(fig1,'PaperType','a4')
    set(fig1,'PaperOrientation','landscape');
    set(fig1,'PaperUnits','normalized');
    set(fig1,'PaperPosition', [0 0 1 1]);
    hold off
end

% Save your TurnOn
save([datadir, ...
    'turnOn.mat'], ...
    'turnOn')

%% Sweep Solver

sweep = struct;

for j = fbAmp

    param = setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(1), ...
        'feed_ampli',j, ...
        'feed_ampliMatrix', feedAmpMat, ...
        'feed_phase',0, ...
        'populate_wrkspc', 0, ...
        'clear',0);

    solnSweep = sweeper(param.J.index, [Jmin, Jmax], param, 'plot',0, ...
        'numSweepSteps', numPoints);

    timeTotal = 0;
    for i = 1:numel(solnSweep)
        timeTotal = timeTotal + solnSweep(i).calcTime;
    end

    disp('It took')
    disp(timeTotal)

    % Add to sweep
    sweep.(['fb',num2str(j*10)]) = solnSweep;
    
    % Plot final intensity versus current amplitude

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
    title(['Sweep Up Time Series: Final intensity vs Current amplitude, k_{ss}=',num2str(j)])
    set(fig2,'PaperType','a4')
    set(fig2,'PaperOrientation','landscape');
    set(fig2,'PaperUnits','normalized');
    set(fig2,'PaperPosition', [0 0 1 1]);
    hold off
end

% Save your sweep
save([datadir, ...
    'sweep.mat'], ...
    'sweep', '-v7.3')


%% Sweep down

downSweep = struct;

%YOU NEED TO MAKE SURE YOU'RE SWEEPING DOWN FROM THE RIGHT VALUE

for j = fbAmp
    initStruct = turnOn.(['fb',num2str(j*10)]);
    
    solnSweepDown = sweeper(...
        initStruct(end).param.J.index, ...
        [Jmax, Jmin], ...
        initStruct(end).param, ...
        'plot',1, ...
        'numSweepSteps', numPoints, ...
        'hist', initStruct(end).timeSeries);
    
    % Add to downSweep
    downSweep.(['fb',num2str(j*10)]) = solnSweepDown;
    
    % Plot final intensity versus current amplitude
    strongFinalIntensityDown= zeros(numPoints,1);
    weakFinalIntensityDown = zeros(numPoints,1);
    for i = 1:numPoints
        strongFinalIntensityDown(i) = norm( ...
            [solnSweepDown(i).timeSeries.y(1,end),...
             solnSweepDown(i).timeSeries.y(2,end)] );
        weakFinalIntensityDown(i) = norm( ...
            [solnSweepDown(i).timeSeries.y(3,end),...
             solnSweepDown(i).timeSeries.y(4,end)] );
    end
    
    fig2 = figure;
    % strong red
    semilogy(Jarray,strongFinalIntensityDown,'r')
    hold on
    % weak blue
    semilogy(Jarray,weakFinalIntensityDown,'b')
    title(['Sweep Down Time Series: Final intensity vs Current amplitude, k_{ss}=',num2str(j)])
    set(fig2,'PaperType','a4')
    set(fig2,'PaperOrientation','landscape');
    set(fig2,'PaperUnits','normalized');
    set(fig2,'PaperPosition', [0 0 1 1]);
    hold off
end

% Save your sweep
save([datadir, ...
    'downSweep.mat'], ...
    'downSweep', '-v7.3')

