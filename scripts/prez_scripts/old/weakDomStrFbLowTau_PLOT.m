%% Weak Dom Plot

clear;

loader(...
    'datadir_specific', ...
    '../data_bimodal-qd-micropillars/weakDomStrFbLowTau/', ...
    'overwrite', 0)


%% Big plotting one:


branchplot = figure;

% V
for i = 5:8 % numel(folds_horBrn006)
    % Add each fold
    if isa(folds_horBrn006(i).error,'double') && folds_horBrn006(i).error == 0
        [~,~,lineHandFold1] = plot_branch(folds_horBrn006(i), param, ...
            'add_2_gcf', 1, 'color','r', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end

% U
for i = 1:numel(folds_vrtBrn00)
    % Add each fold
    if isa(folds_vrtBrn00(i).error,'double') && folds_vrtBrn00(i).error == 0
        plot_branch(folds_vrtBrn00(i), param, ...
            'add_2_gcf', 1, 'color','r', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end

% Hopfs
for i = 1:numel(hopfs_vrtBrn00)
    % Add each hopf
    if i == 1 % problem case
        pt1 = hopfs_vrtBrn00(i);
        pt1.point(314:384) = [];
        pt1.point(end-70:end) = [];
        plot_branch(pt1, param, ...
            'add_2_gcf', 1, 'color','b');
    elseif isa(hopfs_vrtBrn00(i).error,'double') ...
            && hopfs_vrtBrn00(i).error == 0 ...
            && i~= 1
        pt2 = hopfs_vrtBrn00(i);
        pt2.point(end-20:end) = [];
        pt2.point(1:35) = [];
        [~,~,lineHandHopf] = plot_branch(pt2, param, ...
            'add_2_gcf', 1, 'color','b', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end

% 
% plot_branch(horizBranch025Stabl, param, ...
%     'add_2_gcf', 1, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'nunst_color', horizBranch025Stabl.nunst)
% 
% plot_branch(horizBranch035Stabl, param, ...
%     'add_2_gcf', 1, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'nunst_color', horizBranch035Stabl.nunst)



pruned_hopfs = hopfs_horBrn025Stabl;
pruned_hopfs(2).point(end-19:end) = [];
pruned_hopfs(2).point(1:15) = [];
pruned_hopfs(4).point(1:2) = [];

for i = [2,4] % 1:numel(hopfs_horBrn025Stabl)
    % Add each hopf
    if isa(pruned_hopfs(i).error,'double') && pruned_hopfs(i).error == 0
        plot_branch(pruned_hopfs(i), param, ...
            'add_2_gcf', 1, 'color','b', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end


lgnd = legend([lineHandHopf{1}, lineHandFold1{1}], ...
    'Hopf Bifurcation', ...
    'Fold Bifurcation', ...
    'Location', 'SouthEast');
set(gca,'YLim',[0 0.5])
title('Feedback Amp vs Feedback Phase');



%% Save things

% Set + Print to pdf
set(branchplot,'PaperType','a4')
set(branchplot,'PaperOrientation','landscape');
set(branchplot,'PaperUnits','normalized');
set(branchplot,'PaperPosition', [0 0 1 1]);
branchPlotFileName = [master_options.datadir_specific,'BranchPlot.pdf'];
print(branchplot,branchPlotFileName,'-dpdf')