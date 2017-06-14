%% This script will create a turn-on time series
% calc

clear;

% calculate
feedPhaseMat = [1, 1; 1, 1];
feedAmpMat = [0, 0; 0, 1];

param = setup_params_nonDim_CnstCplRatio(...
    'save',0, ...
    'populate_wrkspc', 0, ...
    'alpha_par',0, ...
    'feed_ampli', 0.4, ...
    'feed_ampliMatrix', feedAmpMat, ...
    'feed_phase',0, ...
    'feed_phaseMatrix', feedPhaseMat, ...
    'clear',0,...
    'J',160e-6,...
    'tau_fb', 0.8);

sol = solver([1e-9;0;1e-9;0;0;0], ...
    [0,4], ...
    param, ...
    'plot',0, ...
    'dde23_options',ddeset('RelTol',10^-8), ...
    'quiet', 1);

ef_units = '(\epsilon_{ss} \epsilon_{tilda})^{-1/2}';
time_units = '(\tau_{sp})';
n_units = '(S^{in} \tau_{sp})^{-1}';


%% plot
fig = figure;
clf;

% plot settings
tdeco={'fontsize',30,'fontweight','bold'}; % 
ldeco={'fontsize',26}; % ,'fontweight','bold'
lgndFontSize = 24;

% colors
colorNum = 6;
greycol = brewermap(colorNum,'Greys');
orangecol = brewermap(colorNum,'Oranges');

% linewidth
linewidth = 3;


ax1 = subplot(2,2,[1,2]);
set(ax1,ldeco{:});
ylim([0,0.601])

hold on
plot(sol.x,arrayfun(@(x)norm(x),sol.y(1,:)+1i*sol.y(2,:)), ...
    '-','Color',greycol(5,:),'LineWidth',linewidth)
plot(sol.x,arrayfun(@(x)norm(x),sol.y(3,:)+1i*sol.y(4,:)), ...
    '-', 'Color',orangecol(5,:),'LineWidth',linewidth)
hold off

lgnd1 = legend('E_{s}','E_{w}',ldeco{:});
set(lgnd1,'FontSize',lgndFontSize)

% title({'Electric Field Amplitude vs Time '; ...
%     strcat('with J=',num2str(param.J.value,'%1.1e'),'A');...
%     strcat('|E_s| = Green, |E_w| = red')});
% title('Electric Field Amplitude',tdeco{:})

% xlabel(strcat({'Time '}, time_units),ldeco{:})
ylabel(strcat({'|E(t)| '}, ef_units),ldeco{:})

ax2 = subplot(2,2,[3,4]);

[hAx,hLine1,hLine2] = ...
    plotyy(sol.x,sol.y(5,:), ...
    sol.x,sol.y(6,:));


set(hLine1,'LineWidth',linewidth)
set(hAx(1),'ytick',linspace(0,max(ylim(hAx(1))),3));
set(hLine2,'LineWidth',linewidth)
set(hAx(2),'ytick',linspace(0,max(ylim(hAx(2))),3));
rhsYLabel = get(hAx(2),'YLabel');
set(rhsYLabel,'Position',get(rhsYLabel,'Position') + [0.1 0 0])

set(hAx,ldeco{:});

xlabel(strcat({'Time '}, time_units),ldeco{:})
ylabel(hAx(1),'\rho (no units)', ldeco{:}) % left y-axis 
ylabel(hAx(2),strcat({'n_r '},n_units),ldeco{:}) % right y-axis

% ax2 = subplot(2,2,[3,4]);
% set(ax2,ldeco{:});
% 
% plot(sol.x,sol.y(5,:), ...
%     '-', 'Color','k','LineWidth',linewidth)
% % title('QD Occupation Prob',tdeco{:})
% xlabel(strcat({'Time '}, time_units),ldeco{:})
% ylabel('\rho(t) (no units)', ldeco{:})
% 
% ax3 = subplot(2,2,4);
% set(ax2,ldeco{:});
% 
% plot(sol.x,sol.y(6,:), ...
%     '-', 'Color','k','LineWidth',linewidth)
% title('Carrier Density',tdeco{:})
% xlabel(strcat({'Time '}, time_units),ldeco{:})
% ylabel(strcat({'n_r(t) '},n_units),ldeco{:})


% Page settings & print
set(fig,'PaperType','a4')
set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','normalized');
set(fig,'PaperPosition', [0 0 1 1]);

% prints into current directory
print(fig, ...
    'turn_on_time_series', ...
    '-dpdf')