%% Plot
%% branch plot

branchplot = figure;

% Plot each hopf_branch
namesHopfBranches = fieldnames(hopf_branches);
for i = 1:numel(namesHopfBranches)
    % Add each hopf_branch
    if ~any(strcmp('error', ... 
            fieldnames(hopf_branches.(namesHopfBranches{i}))))
        % Only plot hopf_branches that DO NOT have errors
        plot_branch(hopf_branches.(namesHopfBranches{i}), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end


% Plot each fold_branch
namesFoldBranches = fieldnames(fold_branches);
for i = 1:numel(namesFoldBranches)
    % Add each hopf_branch
    if ~any(strcmp('error',...
            fieldnames(fold_branches.(namesFoldBranches{i}))))
        % Only plot fold_branches that DO NOT have errors
        plot_branch(fold_branches.(namesFoldBranches{i}), param, ...
            'add_2_gcf', 1, 'color','r');
    end
    
end


% Plot init_branch
plot_branch(branch_stst, param, ...
            'add_2_gcf', 1, 'color','g', ...
            'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
        
% Set + Print to pdf
set(branchplot,'PaperType','a4')
set(branchplot,'PaperOrientation','landscape');
set(branchplot,'PaperUnits','normalized');
set(branchplot,'PaperPosition', [0 0 1 1]);
branchPlotFileName = [master_options.datadir_specific,'BranchPlot.pdf'];
print(branchplot,branchPlotFileName,'-dpdf')


%% Plot omega vs continued param with nunst
[~, omega1plot] = plot_branch(branch_stst,param, ...
    'nunst_color',nunst_branch_stst);

% Set + Print to pdf
set(omega1plot,'PaperType','a4')
set(omega1plot,'PaperOrientation','landscape');
set(omega1plot,'PaperUnits','normalized');
set(omega1plot,'PaperPosition', [0 0 1 1]);
omega1PlotFileName = [master_options.datadir_specific,'Omega1Plot.pdf'];
print(omega1plot,omega1PlotFileName,'-dpdf')

% omega2
[~, omega2plot] = plot_branch(branch_stst, param, ...
    'axes_indParam', [ param.feed_phase.index, param.omega2.index ], ...
    'nunst_color', nunst_branch_stst);

% Set + Print to pdf
set(omega2plot,'PaperType','a4')
set(omega2plot,'PaperOrientation','landscape');
set(omega2plot,'PaperUnits','normalized');
set(omega2plot,'PaperPosition', [0 0 1 1]);
omega2PlotFileName = [master_options.datadir_specific,'Omega2Plot.pdf'];
print(omega2plot,omega2PlotFileName,'-dpdf')

% both omegas
[~, omegaBOTHplot] = plot_branch(branch_stst,param, ...
    'nunst_color',nunst_branch_stst); % Put the first omega plot on there
[~, omegaBOTHplot] = plot_branch(branch_stst, param, ...
    'axes_indParam', [ param.feed_phase.index, param.omega2.index ], ...
    'add_2_gcf', 1, ...
    'nunst_color', nunst_branch_stst); % Now the second
title('Both Omega vs Feedback Phase')
ylabel(['Omega 1 and 2', ' (1/tau_{sp})'])

set(omegaBOTHplot,'PaperType','a4')
set(omegaBOTHplot,'PaperOrientation','landscape');
set(omegaBOTHplot,'PaperUnits','normalized');
set(omegaBOTHplot,'PaperPosition', [0 0 1 1]);
omegaBOTHPlotFileName = [master_options.datadir_specific,'OmegaBOTHPlot.pdf'];
print(omegaBOTHplot,omegaBOTHPlotFileName,'-dpdf')

%% Print to combined PDF

unix(['pdftk ', ...
    ['''',branchPlotFileName,''''],' ', ...
    ['''',omega1PlotFileName,''''], ' ', ...
    ['''',omega2PlotFileName,''''], ' ', ...
    ['''',omegaBOTHPlotFileName,''''], ' ', ...
    'cat output ', ...
    ['''',master_options.datadir_specific,''''], 'BranchOmegaPlotCombi.pdf']);






