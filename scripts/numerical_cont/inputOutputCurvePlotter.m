%% General plot settings!!
% plot settings
tdeco={'fontsize',13.2,'fontweight','bold'};

% colors
colorNum = 9;
bluecol = brewermap(colorNum,'Blues');
redcol = brewermap(colorNum,'Reds');
greycol = brewermap(colorNum,'Greys');
greencol = brewermap(colorNum,'Greens');

% linewidth
linewidth = 2;

% save location
datadir = '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/inputOutput/';

% Update Jarray
JarrayScaledTOMicro = Jarray * 1e6;

%% More complicated plotter for TurnOn

turnOnFig = figure('color','w');
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Turn On Time Series: Steady state intensity vs Current amplitude',tdeco{:})
xlabel('Current amplitude (\muA)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

strongPlots = cell(numel(fbAmp),1);
weakPlots = cell(numel(fbAmp),1);

for j = 1:numel(fbAmp)
    solnTurnOn = turnOn.(['fb',num2str(fbAmp(j)*10)]);
    
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

    % strong red
    strongPlots{j} = plot(JarrayScaledTOMicro, strongFinalIntensity,'-', ...
    'Color',redcol(colorNum-8+j,:),'LineWidth',linewidth);

    % weak blue
    weakPlots{j} = plot(JarrayScaledTOMicro, weakFinalIntensity,'-', ...
        'Color',bluecol(colorNum-8+j,:),'LineWidth',linewidth);
    
end


% Page settings
set(turnOnFig,'PaperType','a4')
set(turnOnFig,'PaperOrientation','landscape');
set(turnOnFig,'PaperUnits','normalized');
set(turnOnFig,'PaperPosition', [0 0 1 1]);


% Legend
lgnd1 = legend( [strongPlots{:}], ...
    'I_{s}, k_{ss}=0.0', ...
    'I_{s}, k_{ss}=0.1', ...
    'I_{s}, k_{ss}=0.2', ...
    'I_{s}, k_{ss}=0.3', ...
    'I_{s}, k_{ss}=0.4', ...
    'I_{s}, k_{ss}=0.5');
lgnd1Copy = copyobj(lgnd1,turnOnFig);
set(lgnd1Copy, 'Position', [0.75 0.045 0.2 0.5]); 
lgnd1 = legend( [strongPlots{5}, weakPlots{5}], ...
    'I_{s}', ...
    'I_{w}', ...
    'Location',[0.7 0.11 0.1 0.2]); % 'SouthWest'



hold off


turnOnName = [datadir, 'turnOnMulti','.pdf'];
print(turnOnFig, ...
    turnOnName, ...
    '-dpdf')


%% More complicated plotter for Sweep

sweepUpFig = figure('color','w');
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Sweep Up Time Series: Final intensity vs Current amplitude',tdeco{:})
xlabel('Current amplitude (\muA)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

strongPlots = cell(numel(fbAmp),1);
weakPlots = cell(numel(fbAmp),1);

for j = 1:numel(fbAmp)
    solnSweep = sweep.(['fb',num2str(fbAmp(j)*10)]);
    
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

    % strong red
    strongPlots{j} = plot(JarrayScaledTOMicro, strongFinalIntensitySweep,'--', ...
    'Color',redcol(colorNum-8+j,:),'LineWidth',linewidth);

    % weak blue
    weakPlots{j} = plot(JarrayScaledTOMicro, weakFinalIntensitySweep,'--', ...
        'Color',bluecol(colorNum-8+j,:),'LineWidth',linewidth);
    
    
end

% Page Settings
set(sweepUpFig,'PaperType','a4')
set(sweepUpFig,'PaperOrientation','landscape');
set(sweepUpFig,'PaperUnits','normalized');
set(sweepUpFig,'PaperPosition', [0 0 1 1]);

% Legend
lgnd1 = legend( [strongPlots{:}], ...
    'I_{s}, k_{ss}=0.0', ...
    'I_{s}, k_{ss}=0.1', ...
    'I_{s}, k_{ss}=0.2', ...
    'I_{s}, k_{ss}=0.3', ...
    'I_{s}, k_{ss}=0.4', ...
    'I_{s}, k_{ss}=0.5');
lgnd1Copy = copyobj(lgnd1,sweepUpFig);
set(lgnd1Copy, 'Position', [0.75 0.045 0.2 0.5]); 
lgnd1 = legend( [strongPlots{5}, weakPlots{5}], ...
    'I_{s}', ...
    'I_{w}', ...
    'Location',[0.7 0.11 0.1 0.2]); % 'SouthWest'


hold off


sweepName = [datadir, 'sweepUpMulti','.pdf'];
print(sweepUpFig, ...
    sweepName, ...
    '-dpdf')

%% Sweep Down

sweepDownFig = figure('color','w');
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Sweep Down Time Series: Final intensity vs Current amplitude',tdeco{:})
xlabel('Current amplitude (\muA)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

strongPlots = cell(numel(fbAmp),1);
weakPlots = cell(numel(fbAmp),1);

for j = 1:numel(fbAmp)
    solnDownSweep = downSweep.(['fb',num2str(fbAmp(j)*10)]);
    
    strongFinalIntensityDownSweep = zeros(numPoints,1);
    weakFinalIntensityDownSweep = zeros(numPoints,1);
    for i = 1:numPoints
        strongFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(1,end),...
             solnDownSweep(i).timeSeries.y(2,end)] );
        weakFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(3,end),...
             solnDownSweep(i).timeSeries.y(4,end)] );
    end

    % strong red
    strongPlots{j} = plot(JarrayScaledTOMicro(end:-1:1), strongFinalIntensityDownSweep,'--', ...
    'Color',greycol(colorNum-8+j,:),'LineWidth',linewidth);

    % weak blue
    weakPlots{j} = plot(JarrayScaledTOMicro(end:-1:1), weakFinalIntensityDownSweep,'--', ...
        'Color',greencol(colorNum-8+j,:),'LineWidth',linewidth);
    
    
end

% Page Settings
set(sweepDownFig,'PaperType','a4')
set(sweepDownFig,'PaperOrientation','landscape');
set(sweepDownFig,'PaperUnits','normalized');
set(sweepDownFig,'PaperPosition', [0 0 1 1]);

% Legend
lgnd1 = legend( [strongPlots{:}], ...
    'I_{s}, k_{ss}=0.0', ...
    'I_{s}, k_{ss}=0.1', ...
    'I_{s}, k_{ss}=0.2', ...
    'I_{s}, k_{ss}=0.3', ...
    'I_{s}, k_{ss}=0.4', ...
    'I_{s}, k_{ss}=0.5');
lgnd1Copy = copyobj(lgnd1,sweepDownFig);
set(lgnd1Copy, 'Position', [0.75 0.045 0.2 0.5]); 
lgnd1 = legend( [strongPlots{5}, weakPlots{5}], ...
    'I_{s}', ...
    'I_{w}', ...
    'Location',[0.7 0.11 0.1 0.2]); % 'SouthWest'


hold off

sweepDownName = [datadir, 'sweepDownMulti','.pdf'];
print(sweepDownFig, ...
    sweepDownName, ...
    '-dpdf')


%% Turn On and Sweep Down

turnOnSweepDown = figure('color','w');
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Hysteresis Plot: Steady state intensity vs Current amplitude',tdeco{:})
xlabel('Current amplitude (\muA)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

strongPlotsTurnOn = cell(numel(fbAmp),1);
weakPlotsTurnOn = cell(numel(fbAmp),1);

% turn on
for j = 1:numel(fbAmp)
    solnTurnOn = turnOn.(['fb',num2str(fbAmp(j)*10)]);
    
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

    % strong red
    strongPlotsTurnOn{j} = plot(JarrayScaledTOMicro, strongFinalIntensity,'-', ...
    'Color',redcol(colorNum-8+j,:),'LineWidth',linewidth);

    % weak blue
    weakPlotsTurnOn{j} = plot(JarrayScaledTOMicro, weakFinalIntensity,'-', ...
        'Color',bluecol(colorNum-8+j,:),'LineWidth',linewidth);
    
    
end

strongPlotsSweepDown = cell(numel(fbAmp),1);
weakPlotsSweepDown = cell(numel(fbAmp),1);

% sweep down
for j = 1:numel(fbAmp)
    solnDownSweep = downSweep.(['fb',num2str(fbAmp(j)*10)]);
    
    strongFinalIntensityDownSweep = zeros(numPoints,1);
    weakFinalIntensityDownSweep = zeros(numPoints,1);
    for i = 1:numPoints
        strongFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(1,end),...
             solnDownSweep(i).timeSeries.y(2,end)] );
        weakFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(3,end),...
             solnDownSweep(i).timeSeries.y(4,end)] );
    end

    % strong red
    strongPlotsSweepDown{j} = plot(JarrayScaledTOMicro(end:-1:1), ...
        strongFinalIntensityDownSweep,'--', ...
        'Color',greycol(colorNum-8+j,:),'LineWidth',linewidth);

    % weak blue
    weakPlotsSweepDown{j} = plot(JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep,'--', ...
        'Color',greencol(colorNum-8+j,:),'LineWidth',linewidth);
    
    
end

% Page settings
set(turnOnSweepDown,'PaperType','a4')
set(turnOnSweepDown,'PaperOrientation','landscape');
set(turnOnSweepDown,'PaperUnits','normalized');
set(turnOnSweepDown,'PaperPosition', [0 0 1 1]);


% Legend
lgnd1 = legend( [strongPlotsTurnOn{:}], ...
    'I_{s}, k_{ss}=0.0', ...
    'I_{s}, k_{ss}=0.1', ...
    'I_{s}, k_{ss}=0.2', ...
    'I_{s}, k_{ss}=0.3', ...
    'I_{s}, k_{ss}=0.4', ...
    'I_{s}, k_{ss}=0.5');
lgnd1Copy = copyobj(lgnd1,turnOnSweepDown);
set(lgnd1Copy, 'Position', [0.75 0.045 0.2 0.5]); 
lgnd1 = legend( ...
    [strongPlotsTurnOn{5}, weakPlotsTurnOn{5}, ...
    strongPlotsSweepDown{5}, weakPlotsSweepDown{5}], ...
    'I_{s} TO', ...
    'I_{w} TO', ...
    'I_{s} SD', ...
    'I_{w} SD', ...
    'Location',[0.69 0.142 0.1 0.2]); % 'SouthWest'








hold off

hysteresisName = [datadir, 'hysteresisMulti','.pdf'];
print(turnOnSweepDown, ...
    hysteresisName, ...
    '-dpdf')

%% Combine??

unix(['pdftk ', ...
    ['''',turnOnName,''''],' ', ...
    ['''',sweepName,''''], ' ', ...
    ['''',sweepDownName,''''], ' ', ...
    ['''',hysteresisName,''''], ' ', ...
    'cat output ', ...
    ['''',datadir,''''], 'combi.pdf']);


