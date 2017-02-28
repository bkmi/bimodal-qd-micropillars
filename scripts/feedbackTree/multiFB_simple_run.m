%% For doing a new scan at a different fb amplitude
%% For strongDom fields AKA sweepup

%% load
treeLoader('datadir_specific', ...
    '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/zeroPhaseOffsetJ=560uA/strongDomAlpha=0/tau_fb=0.83ns/');


%% more stst sweeps
stdyBrnchFB = steadyBranchMultiCreator( ...
    funcs,  ...
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6], ...
    param.feed_phase.index, ...
    300, ...
    param, ...
    'sweep',1);

save([master_options.datadir_specific,'stdyBrnchFB'],'stdyBrnchFB')


%% Probe folds

hopfTest = bifurFoldHopfMultiCreator( ...
    funcs, stdyBrnchFB(2), stdyBrnchFB(2).param, stdyBrnchFB(2).indHopf(2), 200);


%% For weakDom fields aka no sweep
% treeLoader('datadir_specific', ...
%     '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/zeroPhaseOffsetJ=560uA/weakDomAlpha=0/tau_fb=0.83ns/');



%% Branch plot

branchplot = figure;

% Plot each hopf_branch
for i = 1:numel(hopfTest)
    % Add each hopf_branch
    if isa(hopfTest.error,'double') && hopfTest.error == 0
        % Only plot hopf_branches that DO NOT have errors
        plot_branch(hopfTest(1), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end


% % Plot each hopf_branch
% namesHopfBranches = fieldnames(hopf_branches);
% for i = 1:numel(namesHopfBranches)
%     % Add each hopf_branch
%     if ~any(strcmp('error', ... 
%             fieldnames(hopf_branches.(namesHopfBranches{i}))))
%         % Only plot hopf_branches that DO NOT have errors
%         plot_branch(hopf_branches.(namesHopfBranches{i}), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end


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
for i = 1:numel(stdyBrnchFB)
    
    plot_branch(stdyBrnchFB(i), stdyBrnchFB(i).param, ...
                'add_2_gcf', 1, 'color','g', ...
                'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
end