% This is a script which first calculates a sweep up, then a sweep down
% given strong to strong feedback. Then it plots the same
%% str2strFB Calc

clear;

% Threshold is ~43e-6 Amps
Jmin = 43e-6;
Jmax = 750e-6;
numCurrentPoints = 150;

Jarray = linspace(...
    Jmin, ...
    Jmax, ...
    numCurrentPoints); 

% feedback amp settings
feedAmpMat = [1, 0; 0, 0];
fbAmp = [0, 0.1, 0.2, 0.3, 0.4, 0.5];

% save location
foldername  = 'str2str_IO';
datadir = ['/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/',foldername,'/'];
mkdir(datadir)

% Save extras: fbAmp, numPoints
save([datadir, ...
    'extras.mat'], ...
    'fbAmp', 'numCurrentPoints', 'Jarray','datadir')





% turn on
turnOn = struct;

solnTurnOn = repmat(...
    struct(...
    'timeSeries',struct, ...
    'J',0,...
    'param',struct, ...
    'calcTime',0),[numCurrentPoints,1]);

for j = fbAmp
    for i = 1:numel(Jarray)
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
    end
    
    % Add to turnOn
    turnOn.(['fb',num2str(j*10)]) = solnTurnOn;
end

% Save your TurnOn
save([datadir, ...
    'turnOn.mat'], ...
    'turnOn')

% This is a better version of organizing the data, but I'm lazy now.
%
% % turnOn is organized s.t. a row is for a particular value of current and
% % each column is for a particular feedback amplitude. Like this:
% % (fbAmp, j)
% % [ (0.0,43e-6), (0.01,43e-6), ... ; ...
% %   (0.0,45e-6), (0.01,45e-6), ... 
% %   ...                            ]
% turnOn = repmat(...
%     struct(...
%     'timeSeries',struct, ...
%     'J',0,...
%     'fbAmp',0, ...
%     'param',struct, ...
%     'calcTime',0),[numel(Jarray),numel(fbAmp)]);
%
% % turn on
% for ind_fbAmp = fbAmp
%     for ind_cur = 1:numel(Jarray)
%         tic; % Start timer
% 
%         param = setup_params_nonDim_CnstCplRatio(...
%             'save',0, ...
%             'J', Jarray(ind_cur), ...
%             'feed_ampli',ind_fbAmp, ...
%             'feed_ampliMatrix', feedAmpMat, ...
%             'feed_phase',0, ...
%             'clear',0,...
%             'populate_wrkspc', 0);
%         dde23_soln = solver([1e-9;0;1e-9;0;0;0], [0,9], param, 'quiet',1);
%         % 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)
% 
%         turnOn(ind_cur, ind_fbAmp).timeSeries = dde23_soln;
%         turnOn(ind_cur, ind_fbAmp).J = param.J.value;
%         turnOn(ind_cur, ind_fbAmp).fbAmp = param.feed_ampli.value;
%         turnOn(ind_cur, ind_fbAmp).param = param;
%         turnOn(ind_cur, ind_fbAmp).calcTime = toc;
%     end
% end





% sweep up
sweep = struct;

for ind_fbAmp = fbAmp
    param = setup_params_nonDim_CnstCplRatio(...
        'save',0, ...
        'J', Jarray(1), ...
        'feed_ampli',ind_fbAmp, ...
        'feed_ampliMatrix', feedAmpMat, ...
        'feed_phase',0, ...
        'populate_wrkspc', 0, ...
        'clear',0);

    solnSweep = sweeper(param.J.index, [Jmin, Jmax], param, 'plot',0, ...
        'numSweepSteps', numCurrentPoints);

    timeTotal = 0;
    for ind_cur = 1:numel(solnSweep)
        timeTotal = timeTotal + solnSweep(ind_cur).calcTime;
    end

    % Add to sweep
    sweep.(['fb',num2str(ind_fbAmp*10)]) = solnSweep;
end

% Save your sweep
save([datadir, ...
    'sweep.mat'], ...
    'sweep', '-v7.3')










% sweep down
downSweep = struct;

%YOU NEED TO MAKE SURE YOU'RE SWEEPING DOWN FROM THE RIGHT VALUE

for j = fbAmp
    initStruct = turnOn.(['fb',num2str(j*10)]);
    
    solnSweepDown = sweeper(...
        initStruct(end).param.J.index, ...
        [Jmax, Jmin], ...
        initStruct(end).param, ...
        'plot',1, ...
        'numSweepSteps', numCurrentPoints, ...
        'hist', initStruct(end).timeSeries);
    
    % Add to downSweep
    downSweep.(['fb',num2str(j*10)]) = solnSweepDown;
end

% Save your sweep
save([datadir, ...
    'downSweep.mat'], ...
    'downSweep', '-v7.3')











%% str2strFB Plot

% plot settings
tdeco={'fontsize',30,'fontweight','bold'}; % 
ldeco={'fontsize',28}; % ,'fontweight','bold'
lgndFontSize = 24;

% colors
colorNum = 6;
greycol = brewermap(colorNum,'Greys');
orangecol = brewermap(colorNum,'Oranges');

% linewidth
linewidth = 3;

% Update Jarray
JarrayScaledTOMicro = Jarray * 1e6;
flipedJ = fliplr(JarrayScaledTOMicro);

% Relevant fbAmps
relFBamps = fbAmp([1,3,5]);





































% Two plots for up and down, NO FB

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
set(ax1,'FontSize',lgndFontSize);
hold(ax1,'on')

title('No feedback - Sweep Up',ldeco{:})
% xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepUp = cell(numel(relFBamps),1);
weakPlotsSweepUp = cell(numel(relFBamps),1);

% sweep up
for j = 1 % :numel(relFBamps) % NO FB
    solnSweep = sweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensitySweep = zeros(numCurrentPoints,1);
    weakFinalIntensitySweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
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
set(ax2,'FontSize',lgndFontSize);
hold(ax2,'on')

title('Sweep Down',ldeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepDown = cell(numel(relFBamps),1);
weakPlotsSweepDown = cell(numel(relFBamps),1);

% sweep down
for j = 1 % :numel(relFBamps) % NO FB
    solnDownSweep = downSweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    weakFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
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
        flipedJ, ...
        strongFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        flipedJ, ...
        weakFinalIntensityDownSweep * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    

%     % patch option
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

% %Legend
% lgnd2 = legend( ax2, ...
%     [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
%     'k_{s}=0.0', ...
%     'k_{w}=0.0', ...
%     'Location','NorthWest');
% set(lgnd2,'FontSize',lgndFontSize)


% Add color rectangle
xRectnoFB = ...
    [(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % bot left
    (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % top left
    750, ... % bot right
    750]; % top right
yRectnoFB = ...
    [1e-3, 1e0, 1e0, 1e-3];
rec1 = patch(ax1, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
rec2 = patch(ax2, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
text(ax1,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,8e-1, ...
    'Bistable', ...
    'FontSize',lgndFontSize)
text(ax2,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,8e-1, ...
    'Bistable', ...
    'FontSize',lgndFontSize)


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

noFBName = [datadir, 'noFB'];
print(noFB, ...
    noFBName, ...
    '-dpdf')
print(noFB, ...
    noFBName, ...
    '-dpng')






















% Two plots for up and down, MID FB

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

title('Low Feedback - Sweep Up',tdeco{:})
% xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepUp = cell(numel(relFBamps),1);
weakPlotsSweepUp = cell(numel(relFBamps),1);

% sweep up
for j = 2 % :numel(relFBamps) % NO FB
    solnSweep = sweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensitySweep = zeros(numCurrentPoints,1);
    weakFinalIntensitySweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
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

title('Sweep Down',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepDown = cell(numel(relFBamps),1);
weakPlotsSweepDown = cell(numel(relFBamps),1);

% sweep down
for j = 2 % :numel(relFBamps) % NO FB
    solnDownSweep = downSweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    weakFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
        strongFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(1,end),...
             solnDownSweep(i).timeSeries.y(2,end)] );
        weakFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(3,end),...
             solnDownSweep(i).timeSeries.y(4,end)] );
    end

    % switch between 99 and 100
    cutit = 99;
    % strong grey
    strongPlotsSweepDown{j} = plot( ...
        flipedJ(1:cutit), ...
        strongFinalIntensityDownSweep(1:cutit) * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);
    plot( ...
        flipedJ(cutit+1:end), ...
        strongFinalIntensityDownSweep(cutit+1:end) * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);
    
    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        flipedJ(1:cutit), ...
        weakFinalIntensityDownSweep(1:cutit) * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    plot( ...
        flipedJ(cutit+1:end), ...
        weakFinalIntensityDownSweep(cutit+1:end) * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% %Legend
% lgnd2 = legend( ax2, ...
%     [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
%     'k_{s}=0.2', ...
%     'k_{w}=0.0', ...
%     'Location','NorthWest');
% set(lgnd2,'FontSize',lgndFontSize)




% Add first color rectangle
xRectnoFB = ...
    [(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % bot left
    (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % top left
    750, ... % bot right
    750]; % top right
yRectnoFB = ...
    [1e-3, 1e0, 1e0, 1e-3];
rec1 = patch(ax1, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
rec2 = patch(ax2, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
% text(ax1,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,2e-2, ...
%     'Without feedback', ...
%     'FontSize',lgndFontSize)
% text(ax2,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,2e-2, ...
%     'Without feedback', ...
%     'FontSize',lgndFontSize)

% Add second color rectangle
xRectnoFB = ...
    [(flipedJ(99)+flipedJ(100))/2, ... % bot left
    (flipedJ(99)+flipedJ(100))/2, ... % top left
    750, ... % bot right
    750]; % top right
yRectnoFB = ...
    [1e-3, 1e0, 1e0, 1e-3];
rec1 = patch(ax1, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
rec2 = patch(ax2, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
text(ax1,(flipedJ(99)+flipedJ(100))/2,8e-1, ...
    'Bistable', ...
    'FontSize',lgndFontSize)
text(ax2,(flipedJ(99)+flipedJ(100))/2,8e-1, ...
    'Bistable', ...
    'FontSize',lgndFontSize)




set(gcf,'NextPlot','add');
axes;
h = ylabel('Steady state intensity (a.u.)',ldeco{:});
set(gca,'Visible','off');
set(h,'Visible','on'); 
pos = get(h,'pos'); % Read position [x y z]
set(h,'pos',pos+[pos(1)*.4 0 0]) % Move label




% Page settings
set(midFB,'PaperType','a4')
set(midFB,'PaperOrientation','landscape');
set(midFB,'PaperUnits','normalized');
set(midFB,'PaperPosition', [0 0 1 1]);

midFBName = [datadir, 'midFB'];
print(midFB, ...
    midFBName, ...
    '-dpdf')
print(midFB, ...
    midFBName, ...
    '-dpng')































% Two plots for up and down, MAX FB

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

title('High Feedback - Sweep Up',tdeco{:})
% xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepUp = cell(numel(relFBamps),1);
weakPlotsSweepUp = cell(numel(relFBamps),1);

% sweep up
for j = 3 % :numel(relFBamps) % NO FB
    solnSweep = sweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensitySweep = zeros(numCurrentPoints,1);
    weakFinalIntensitySweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
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

title('Sweep Down',tdeco{:})
xlabel('Current amplitude (\muA)',ldeco{:})
% ylabel('Steady state intensity (a.u.)',ldeco{:})

strongPlotsSweepDown = cell(numel(relFBamps),1);
weakPlotsSweepDown = cell(numel(relFBamps),1);

% sweep down
for j = 3 % :numel(relFBamps) % NO FB
    solnDownSweep = downSweep.(['fb',num2str(relFBamps(j)*10)]);
    
    strongFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    weakFinalIntensityDownSweep = zeros(numCurrentPoints,1);
    for i = 1:numCurrentPoints
        strongFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(1,end),...
             solnDownSweep(i).timeSeries.y(2,end)] );
        weakFinalIntensityDownSweep(i) = norm( ...
            [solnDownSweep(i).timeSeries.y(3,end),...
             solnDownSweep(i).timeSeries.y(4,end)] );
    end

    % between 14 and 15
    cutit = 14;
    % strong grey
    strongPlotsSweepDown{j} = plot( ...
        flipedJ(1:cutit), ...
        strongFinalIntensityDownSweep(1:cutit) * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);
    plot( ...
        flipedJ(cutit+1:end), ...
        strongFinalIntensityDownSweep(cutit+1:end) * ySclFactor, ...
        '-', ...
        'Color','k','LineWidth',linewidth);

    % weak orange
    weakPlotsSweepDown{j} = plot( ...
        flipedJ(1:cutit), ...
        weakFinalIntensityDownSweep(1:cutit) * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    plot( ...
        flipedJ(cutit+1:end), ...
        weakFinalIntensityDownSweep(cutit+1:end) * ySclFactor, ...
        '-', ...
        'Color',orangecol(4,:),'LineWidth',linewidth);
    
    
end

% Set axis
axis(axisVector)

% %Legend
% lgnd2 = legend( ax2, ...
%     [strongPlotsSweepDown{3}, weakPlotsSweepDown{3}], ...
%     'k_{s}=0.4', ...
%     'k_{w}=0.0', ...
%     'Location','NorthWest');
% set(lgnd2,'FontSize',lgndFontSize)




% Add first color rectangle
xRectnoFB = ...
    [(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % bot left
    (JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2, ... % top left
    750, ... % bot right
    750]; % top right
yRectnoFB = ...
    [1e-3, 1e0, 1e0, 1e-3];
rec1 = patch(ax1, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
rec2 = patch(ax2, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
% text(ax1,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,2e-2, ...
%     'Without feedback', ...
%     'FontSize',lgndFontSize)
% text(ax2,(JarrayScaledTOMicro(18)+JarrayScaledTOMicro(19))/2,2e-2, ...
%     'Without feedback', ...
%     'FontSize',lgndFontSize)

% Add second color rectangle
xRectnoFB = ...
    [(flipedJ(14)+flipedJ(15))/2, ... % bot left
    (flipedJ(14)+flipedJ(15))/2, ... % top left
    750, ... % bot right
    750]; % top right
yRectnoFB = ...
    [1e-3, 1e0, 1e0, 1e-3];
rec1 = patch(ax1, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');
rec2 = patch(ax2, ...
    xRectnoFB, yRectnoFB, ...
    'blue', ...
    'FaceColor','blue', ...
    'FaceAlpha',.25, ...
    'EdgeColor', 'none');

axText = axes('Position',[0 0 1 1],'Visible','off');
text(axText, .8, 0.45, ...
    'Bistable', ...
    'FontSize',lgndFontSize)
text(axText, .8, 0.925, ...
    'Bistable', ...
    'FontSize',lgndFontSize)


% text(ax1,(flipedJ(14)+flipedJ(15))/2+4,4e-1, ...
%     'Bistable', ...
%     'FontSize',lgndFontSize)
% text(ax2,(flipedJ(14)+flipedJ(15))/2+4,4e-1, ...
%     'Bistable', ...
%     'FontSize',lgndFontSize)




set(gcf,'NextPlot','add');
axes;
h = ylabel('Steady state intensity (a.u.)',ldeco{:});
set(gca,'Visible','off');
set(h,'Visible','on'); 
pos = get(h,'pos'); % Read position [x y z]
set(h,'pos',pos+[pos(1)*.4 0 0]) % Move label




% Page settings
set(maxFB,'PaperType','a4')
set(maxFB,'PaperOrientation','landscape');
set(maxFB,'PaperUnits','normalized');
set(maxFB,'PaperPosition', [0 0 1 1]);

maxFBName = [datadir, 'maxFB'];
print(maxFB, ...
    maxFBName, ...
    '-dpdf')
print(maxFB, ...
    maxFBName, ...
    '-dpng')
























% Combine??

unix(['pdftk ', ...
    ['''',noFBName,''''],' ', ...
    ['''',midFBName,''''], ' ', ...
    ['''',maxFBName,''''], ' ', ...
    'cat output ', ...
    ['''',datadir,''''], 'combi.pdf']);

















