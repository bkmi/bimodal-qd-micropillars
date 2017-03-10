%% Strong Dom Plot

clear;

loader(...
    'datadir_specific', ...
    '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/strDomStrFbLowTau/', ...
    'overwrite', 0)


%% Plot

branchplot = figure;

for i = [3,4,5,6] % numel(folds_horBrn015)
    % Add each fold
    if isa(folds_horBrn015(i).error,'double') && folds_horBrn015(i).error == 0
        plot_branch(folds_horBrn015(i), param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

% TOO OUT OF THE WAY
% for i = [3] % numel(hopfs_horBrn05)
%     % Add each hopf
%     if isa(hopfs_horBrn05(i).error,'double') && hopfs_horBrn05(i).error == 0
%         plot_branch(hopfs_horBrn05(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end

for i = [1,2,3,4 ] % 1:numel(hopfs_horBrn045)
    % Add each hopf
    if isa(hopfs_horBrn045(i).error,'double') && hopfs_horBrn045(i).error == 0
        plot_branch(hopfs_horBrn045(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

for i = 1:numel(hopfs_vrtBrn905)
    % Add each hopf
    if isa(hopfs_vrtBrn905(i).error,'double') && hopfs_vrtBrn905(i).error == 0
        plot_branch(hopfs_vrtBrn905(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

for i = 1:numel(hopfs_horBrn04Unst)
    % Add each hopf
    if isa(hopfs_horBrn04Unst(i).error,'double') && hopfs_horBrn04Unst(i).error == 0
        plot_branch(hopfs_horBrn04Unst(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

for i = 1
    % Add each hopf
    if isa(leftHighHopf(i).error,'double') && leftHighHopf(i).error == 0
        plot_branch(leftHighHopf(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end