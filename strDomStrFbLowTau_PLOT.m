%% Strong Dom Plot

clear;

loader(...
    'datadir_specific', ...
    '../data_bimodal-qd-micropillars/strDomStrFbLowTau/', ...
    'overwrite', 0)


%% Plot

branchplot = figure;

for i = [3,4,5,6] % numel(folds_horBrn015)
    % Add each fold
    if isa(folds_horBrn015(i).error,'double') && folds_horBrn015(i).error == 0
        [~,~,lineHandFold1] = plot_branch(folds_horBrn015(i), param, ...
            'add_2_gcf', 1, 'color','r', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end

% TOO OUT OF THE WAY
% for i = [3] % numel(hopfs_horBrn05)
%     % Add each hopf
%     if isa(hopfs_horBrn05(i).error,'double') && hopfs_horBrn05(i).error == 0
%         plot_branch(hopfs_horBrn05(i), param, ...
%             'add_2_gcf', 1, 'color','b', ...
%             'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
%     end
% end

pruned_hopfs = hopfs_horBrn045;
pruned_hopfs(2).point(end-120:end) = [];
pruned_hopfs(3).point(end-100:end) = [];

for i = [1,2,3,4] % 1:numel(hopfs_horBrn045)
    % Add each hopf
    if isa(pruned_hopfs(i).error,'double') && pruned_hopfs(i).error == 0
        [~,~,lineHandHopf] = plot_branch(pruned_hopfs(i), param, ...
            'add_2_gcf', 1, 'color','b', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end

pruned_hopfs = hopfs_vrtBrn905;
pruned_hopfs(3).point(end-10:end) = [];

for i = [3,4] % 1:numel(hopfs_vrtBrn905)
    % Add each hopf
    if isa(pruned_hopfs(i).error,'double') && pruned_hopfs(i).error == 0
        plot_branch(pruned_hopfs(i), param, ...
            'add_2_gcf', 1, 'color','b', ...
            'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
    end
end


pruned_hopfs = hopfs_horBrn04Unst;
pruned_hopfs(1).point(end-4:end) = [];
pruned_hopfs(4).point(end-39:end) = [];
pruned_hopfs(6).point(end-140:end) = [];
pruned_hopfs(8).point(end-249:end) = [];

for i = [1,2,3,4,5,6,7,8] % 1:numel(pruned_hopfs) 
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


% for i = 1
%     % Add each hopf
%     if isa(leftHighHopf(i).error,'double') && leftHighHopf(i).error == 0
%         plot_branch(leftHighHopf(i), param, ...
%             'add_2_gcf', 1, 'color','b', ...
%             'PlotStyle', { 'LineStyle', '-', 'Marker', '.' });
%     end
% end