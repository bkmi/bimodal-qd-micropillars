%% General plot settings!!
% plot settings
tdeco={'fontsize',13.2,'fontweight','bold'};

% colors
colorNum = 9;
bluecol = brewermap(colorNum,'Blues');
redcol = brewermap(colorNum,'Reds');

% linewidth
linewidth = 2;

% save location
datadir = '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/inputOutput/';

%% More complicated plotter for TurnOn

turnOnFig = figure;
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Turn On Time Series: Steady state intensity vs Current amplitude',tdeco{:})
xlabel('Current Amplitude (A)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

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
    plot(Jarray, strongFinalIntensity,'--', ...
    'Color',redcol(colorNum-8+j,:),'LineWidth',linewidth)

    % weak blue
    plot(Jarray, weakFinalIntensity,'--', ...
        'Color',bluecol(colorNum-8+j,:),'LineWidth',linewidth)
    
    
end

% lgnd = legend(...
%     '|E_{s}|^2, k_{ss}=0', ...
%     '|E_{w}|^2, k_{ss}=0', ...
%     '|E_{s}|^2, k_{ss}=0.1', ...
%     '|E_{w}|^2, k_{ss}=0.1', ...
%     '|E_{s}|^2, k_{ss}=0.2', ...
%     '|E_{w}|^2, k_{ss}=0.2', ...
%     '|E_{s}|^2, k_{ss}=0.3', ...
%     '|E_{w}|^2, k_{ss}=0.3', ...
%     '|E_{s}|^2, k_{ss}=0.4', ...
%     '|E_{w}|^2, k_{ss}=0.4', ...
%     '|E_{s}|^2, k_{ss}=0.5', ...
%     '|E_{w}|^2, k_{ss}=0.5');
% set(lgnd,'FontSize',10,'FontWeight','normal')

set(turnOnFig,'PaperType','a4')
set(turnOnFig,'PaperOrientation','landscape');
set(turnOnFig,'PaperUnits','normalized');
set(turnOnFig,'PaperPosition', [0 0 1 1]);
hold off


turnOnName = [datadir, 'turnOnMulti','.pdf'];
print(turnOnFig, ...
    turnOnName, ...
    '-dpdf')


%% More complicated plotter for Sweep

sweepUpFig = figure;
set(gca,'yscale','log')
set(gca,tdeco{:});

hold on
title('Sweep Up Time Series: Final intensity vs Current amplitude',tdeco{:})
xlabel('Current Amplitude (A)',tdeco{:})
ylabel('Steady state intensity (a.u.)',tdeco{:})

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
    plot(Jarray, strongFinalIntensitySweep,'--', ...
    'Color',redcol(colorNum-8+j,:),'LineWidth',linewidth)

    % weak blue
    plot(Jarray, weakFinalIntensitySweep,'--', ...
        'Color',bluecol(colorNum-8+j,:),'LineWidth',linewidth)
    
    
end

set(sweepUpFig,'PaperType','a4')
set(sweepUpFig,'PaperOrientation','landscape');
set(sweepUpFig,'PaperUnits','normalized');
set(sweepUpFig,'PaperPosition', [0 0 1 1]);
hold off


sweepName = [datadir, 'sweepUpMulti','.pdf'];
print(sweepUpFig, ...
    sweepName, ...
    '-dpdf')

%% Combine??

unix(['pdftk ', ...
    ['''',turnOnName,''''],' ', ...
    ['''',sweepName,''''], ' ', ...
    'cat output ', ...
    ['''',datadir,''''], 'combi.pdf']);
