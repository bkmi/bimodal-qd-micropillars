%% General plot settings!!
% You need to generate some data, use inputOutput script in 'old'

% plot settings
tdeco={'fontsize',20,'fontweight','bold'}; % 
ldeco={'fontsize',18}; % ,'fontweight','bold'
lgndFontSize = 15;

% colors
colorNum = 6;
bluecol = brewermap(colorNum,'Blues');
redcol = brewermap(colorNum,'Reds');
greencol = brewermap(colorNum,'Greens');

greycol = brewermap(colorNum,'Greys');
orangecol = brewermap(colorNum,'Oranges');

% linewidth
linewidth = 3;

% save location
datadir = '../data_bimodal-qd-micropillars/inputOutput/';

% Update Jarray
JarrayScaledTOMicro = Jarray * 1e6;

% Relevant fbAmps
relFBamps = fbAmp([1,3,5]);

%% Two plots for up and down, no FB

noFB = figure('color','w');
ySclFactor = 0.79;
axisVector = ...
    [0 ...  % xmin
    800 ... % xmax
    1e-2 ...% ymin
    1e0];

% Sweep up side
ax1 = subplot(2,1,1);

set(ax1,'yscale','log')
set(ax1,ldeco{:});
hold(ax1,'on')

title('Steady state intensity vs Current amplitude: Sweep Up',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepUp = cell(numel(relFBamps),1);
weakPlotsSweepUp = cell(numel(relFBamps),1);

% sweep up
for j = 1 % :numel(relFBamps) % NO FB
    solnSweep = sweep.(['fb',num2str(relFBamps(j)*10)]);
    
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

    % strong grey
    strongPlotsSweepUp{j} = plot( ...
        JarrayScaledTOMicro, ...
        strongFinalIntensitySweep * ySclFactor, ...
        '-', ...
        'Color',greycol(4,:),'LineWidth',linewidth);

    % weak orange
    weakPlotsSweepUp{j} = plot(...
        JarrayScaledTOMicro, ...
        weakFinalIntensitySweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% Legend
lgnd1 = legend( ax1, ...
    [strongPlotsSweepUp{3}, weakPlotsSweepUp{3}], ...
    'I_{s}', ...
    'I_{w}', ...
    'Location','NorthWest'); % [0.7 0.11 0.1 0.2]
set(lgnd1,'FontSize',lgndFontSize)




% Sweep down side
ax2 = subplot(2,1,2);

set(ax2,'yscale','log')
set(ax2,ldeco{:});
hold(ax2,'on')

title('Steady state intensity vs Current amplitude: Sweep Down',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepDown = cell(numel(relFBamps),1);
weakPlotsSweepDown = cell(numel(relFBamps),1);

% sweep down
for j = 1 % :numel(relFBamps) % NO FB
    solnDownSweep = downSweep.(['fb',num2str(relFBamps(j)*10)]);
    
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

    % strong grey
    strongPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        strongFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',greycol(4,:),'LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% Legend
% lgnd2 = legend( ax2, ...
%     [strongPlotsSweepDown{:}], ...
%     'k=0.0', ...
%     'k=0.2', ...
%     'k=0.4', ...
%     'Location','NorthWest');
% set(lgnd2,'FontSize',lgndFontSize)






% Page settings
set(noFB,'PaperType','a4')
set(noFB,'PaperOrientation','landscape');
set(noFB,'PaperUnits','normalized');
set(noFB,'PaperPosition', [0 0 1 1]);

noFBName = [datadir, 'noFB','.pdf'];
print(noFB, ...
    noFBName, ...
    '-dpdf')


%% Two plots for up and down, 3 FB levels

twoPlots = figure('color','w');
ySclFactor = 0.79;
axisVector = ...
    [0 ...  % xmin
    800 ... % xmax
    1e-2 ...% ymin
    1e0];

% Sweep up side
ax1 = subplot(2,1,1);

set(ax1,'yscale','log')
set(ax1,ldeco{:});
hold(ax1,'on')

title('Steady state intensity vs Current amplitude: Sweep Up',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepUp = cell(numel(relFBamps),1);
weakPlotsSweepUp = cell(numel(relFBamps),1);

% sweep up
for j = 1:numel(relFBamps)
    solnSweep = sweep.(['fb',num2str(relFBamps(j)*10)]);
    
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

    % strong grey
    strongPlotsSweepUp{j} = plot( ...
        JarrayScaledTOMicro, ...
        strongFinalIntensitySweep * ySclFactor, ...
        '-', ...
        'Color',greycol(j+2,:),'LineWidth',linewidth);

    % weak orange
    weakPlotsSweepUp{j} = plot(...
        JarrayScaledTOMicro, ...
        weakFinalIntensitySweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(j+2,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% Legend
lgnd1 = legend( ax1, ...
    [strongPlotsSweepUp{3}, weakPlotsSweepUp{3}], ...
    'I_{s}', ...
    'I_{w}', ...
    'Location','NorthWest'); % [0.7 0.11 0.1 0.2]
set(lgnd1,'FontSize',lgndFontSize)




% Sweep down side
ax2 = subplot(2,1,2);

set(ax2,'yscale','log')
set(ax2,ldeco{:});
hold(ax2,'on')

title('Steady state intensity vs Current amplitude: Sweep Down',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepDown = cell(numel(relFBamps),1);
weakPlotsSweepDown = cell(numel(relFBamps),1);

% sweep down
for j = 1:numel(relFBamps)
    solnDownSweep = downSweep.(['fb',num2str(relFBamps(j)*10)]);
    
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

    % strong grey
    strongPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        strongFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',greycol(j+2,:),'LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(j+2,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% Legend
lgnd2 = legend( ax2, ...
    [strongPlotsSweepDown{:}], ...
    'k=0.0', ...
    'k=0.2', ...
    'k=0.4', ...
    'Location','NorthWest');
set(lgnd2,'FontSize',lgndFontSize)






% Page settings
set(twoPlots,'PaperType','a4')
set(twoPlots,'PaperOrientation','landscape');
set(twoPlots,'PaperUnits','normalized');
set(twoPlots,'PaperPosition', [0 0 1 1]);

twoPlotName = [datadir, 'twoPlotMulti','.pdf'];
print(twoPlots, ...
    twoPlotName, ...
    '-dpdf')




%% Combine??

unix(['pdftk ', ...
    ['''',noFBName,''''],' ', ...
    ['''',twoPlotName,''''], ' ', ...
    'cat output ', ...
    ['''',datadir,''''], 'combi.pdf']);