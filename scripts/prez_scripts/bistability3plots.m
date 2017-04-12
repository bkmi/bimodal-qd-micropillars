%% General plot settings!!
% You need to generate some data, use inputOutput script in 'old'

% plot settings
tdeco={'fontsize',25,'fontweight','bold'}; % 
ldeco={'fontsize',23}; % ,'fontweight','bold'
lgndFontSize = 18;

% colors
colorNum = 6;
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


%% Two plots for up and down, NO FB

noFB = figure('color','w');
ySclFactor = 0.79;
axisVector = ...
    [0 ...  % xmin
    750 ... % xmax
    1e-2 ...% ymin
    1e0];

% Sweep up side
ax1 = subplot(2,1,1);

set(ax1,'yscale','log')
set(ax1,ldeco{:});
hold(ax1,'on')

title('No feedback - Sweep Up',ldeco{:})
% xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

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
        'Color','k','LineWidth',linewidth);

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

title('Sweep Down',ldeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

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
    
    % switch btwn end-18 and end-19
    % strong grey
    strongPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        strongFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    % Add color rectangle
%     xRectnoFB = ...
%         [(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % bot left
%         (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % top left
%         750, ... % bot right
%         750]; % top right
%     yRectnoFB = ...
%         [1e-3, 1e0, 1e0, 1e-3];
%     rec1 = patch(ax1, ...
%         xRectnoFB, yRectnoFB, ...
%         'blue', ...
%         'FaceColor','blue', ...
%         'FaceAlpha',.5);
%     rec2 = patch(ax2, ...
%         xRectnoFB, yRectnoFB, ...
%         'blue', ...
%         'FaceColor','blue', ...
%         'FaceAlpha',.5);
    
%     vertNOFB = ...
%         [ (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, 1e-3; ...
%         (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, 1e0; ...
%         750, 1e0; ...
%         750, 1e-3 ];
%     faceNOFB = [1 2 3 4];
%     patch(ax1, ...
%         'Faces',faceNOFB, ...
%         'Vertices',vertNOFB, ...
%         'FaceColor','blue', ...
%         'FaceAlpha',.1)
%     patch(ax2, ...
%         'Faces',faceNOFB, ...
%         'Vertices',vertNOFB, ...
%         'FaceColor','blue', ...
%         'FaceAlpha',.1)
    
    
end

% Set axis
axis(axisVector)

%Legend
lgnd2 = legend( ax2, ...
    [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
    'k_{s}=0.0', ...
    'k_{w}=0.0', ...
    'Location','NorthWest');
set(lgnd2,'FontSize',lgndFontSize)



set(gcf,'NextPlot','add');
axes;
h = ylabel('Steady state intensity (a.u.)',ldeco{:});
set(gca,'Visible','off');
set(h,'Visible','on'); 
pos = get(h,'pos'); % Read position [x y z]
set(h,'pos',pos+[pos(1)*.4 0 0]) % Move label

% t = title('No feedback',tdeco{:}); 
% set(t,'Visible','on'); 
% pos = get(t,'pos'); % Read position [x y z]
% set(t,'pos',pos+[0 0 0]) % Move title 


% Page settings
set(noFB,'PaperType','a4')
set(noFB,'PaperOrientation','landscape');
set(noFB,'PaperUnits','normalized');
set(noFB,'PaperPosition', [0 0 1 1]);

midFBName = [datadir, 'noFB','.pdf'];
print(noFB, ...
    midFBName, ...
    '-dpdf')


%% Two plots for up and down, MID FB

midFB = figure('color','w');
ySclFactor = 0.79;
axisVector = ...
    [0 ...  % xmin
    750 ... % xmax
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
for j = 2 % :numel(relFBamps) % NO FB
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
        'Color','k','LineWidth',linewidth);

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
for j = 2 % :numel(relFBamps) % NO FB
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
        'Color','k','LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

%Legend
lgnd2 = legend( ax2, ...
    [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
    'k_{s}=0.2', ...
    'k_{w}=0.0', ...
    'Location','NorthWest');
set(lgnd2,'FontSize',lgndFontSize)






% Page settings
set(midFB,'PaperType','a4')
set(midFB,'PaperOrientation','landscape');
set(midFB,'PaperUnits','normalized');
set(midFB,'PaperPosition', [0 0 1 1]);

midFBName = [datadir, 'midFB','.pdf'];
print(midFB, ...
    midFBName, ...
    '-dpdf')




%% Two plots for up and down, MAX FB

maxFB = figure('color','w');
ySclFactor = 0.79;
axisVector = ...
    [0 ...  % xmin
    750 ... % xmax
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
for j = 3 % :numel(relFBamps) % NO FB
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
        'Color','k','LineWidth',linewidth);

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
for j = 3 % :numel(relFBamps) % NO FB
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
        'Color','k','LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        JarrayScaledTOMicro(end:-1:1), ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

%Legend
lgnd2 = legend( ax2, ...
    [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
    'k_{s}=0.4', ...
    'k_{w}=0.0', ...
    'Location','NorthWest');
set(lgnd2,'FontSize',lgndFontSize)






% Page settings
set(maxFB,'PaperType','a4')
set(maxFB,'PaperOrientation','landscape');
set(maxFB,'PaperUnits','normalized');
set(maxFB,'PaperPosition', [0 0 1 1]);

maxFBName = [datadir, 'maxFB','.pdf'];
print(maxFB, ...
    maxFBName, ...
    '-dpdf')







%% Combine??

unix(['pdftk ', ...
    ['''',noFBName,''''],' ', ...
    ['''',midFBName,''''], ' ', ...
    ['''',maxFBName,''''], ' ', ...
    'cat output ', ...
    ['''',datadir,''''], 'combi.pdf']);
