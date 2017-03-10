%% Weak Dom Plot

clear;

loader(...
    'datadir_specific', ...
    '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/weakDomStrFbLowTau/', ...
    'overwrite', 0)


%% Big plotting one:


branchplot = figure;

% V
for i = 5:8 % numel(folds_horBrn006)
    % Add each fold
    if isa(folds_horBrn006(i).error,'double') && folds_horBrn006(i).error == 0
        plot_branch(folds_horBrn006(i), param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

% U
for i = 1:numel(folds_vrtBrn00)
    % Add each fold
    if isa(folds_vrtBrn00(i).error,'double') && folds_vrtBrn00(i).error == 0
        plot_branch(folds_vrtBrn00(i), param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

% Hopfs
for i = 1:numel(hopfs_vrtBrn00)
    % Add each hopf
    if i == 1 % problem case
        pt1 = hopfs_vrtBrn00(i);
        pt1.point(314:384) = [];
        plot_branch(pt1, param, ...
            'add_2_gcf', 1, 'color','c');
    elseif isa(hopfs_vrtBrn00(i).error,'double') ...
            && hopfs_vrtBrn00(i).error == 0 ...
            && i~= 1
        plot_branch(hopfs_vrtBrn00(i), param, ...
            'add_2_gcf', 1, 'color','c');
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


for i = [2,3,4] % 1:numel(hopfs_horBrn025Stabl)
    % Add each hopf
    if isa(hopfs_horBrn025Stabl(i).error,'double') && hopfs_horBrn025Stabl(i).error == 0
        plot_branch(hopfs_horBrn025Stabl(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end
