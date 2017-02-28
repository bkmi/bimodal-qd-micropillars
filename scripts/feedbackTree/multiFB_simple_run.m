%% For doing a new scan at a different fb amplitude
% For strongDom fields
% Set our relevant names
% appendToName = '_fb=0_1';
% master_options.save = 0;
% 
% param_fb0_01 = updateParams(param, 'feed_ampli', 0.1);
% 
% sweepSoln_fb0_01 = sweeper(param.J.index, [60e-6, param.J.value], param_fb0_01, master_options);
% dde23_soln_fb0_01 = sweepSoln_fb0_01(end).timeSeries;
% 
% save([master_options.datadir_specific,'dde23_fb0_01.mat'],'sweepSoln_fb0_01','dde23_soln_fb0_01')
% 
% 
% %%
% [branch_stst_fb0_01, nunst_branch_stst_fb0_01, ind_fold_fb0_01, ind_hopf_fb0_01] = ... 
%     init_branch(funcs, ...
%     dde23_soln.y(:,end), param.feed_phase.index, 300, param_fb0_01, ...
%     'max_step',[param.feed_phase.index, (1)*pi/32], ...
%     'minimal_real_part', -0.5, ...
%     'reverse', 1, ...
%     master_options);









% Other steady state branches
feedbackArray = [0.1, 0.2];
steadyBranches_FB = struct('method',struct,'parameter',struct,'point',struct, ...
    'param', struct, 'nunst', 0, 'indFold', 0, 'indHopf', 0, 'sweep', struct);
steadyBranches_FB = repmat(steadyBranches_FB, [numel(feedbackArray), 1]);


for i = 1:numel(feedbackArray)
    % Temporary names and data
    steadyBranch = struct;
    steadyBranchParam = updateParams(param, 'feed_ampli', feedbackArray(i));
    
    % Sweep Up
    steadyBranchSweep = sweeper(steadyBranchParam.J.index, ...
        [60e-6, steadyBranchParam.J.value], steadyBranchParam, ...
        master_options, ...
        'save', 0);

    
    try
        [steadyBranchSTST, nunst_steadyBranchSTST, ...
            ind_foldSteadyBranchSTST, ind_hopfSteadyBranchSTST] = ... 
            init_branch(funcs, ...
            dde23_soln.y(:,end), steadyBranchParam.feed_phase.index, ...
            400, steadyBranchParam, ...
            'max_step',[steadyBranchParam.feed_phase.index, (1)*pi/32], ...
            'minimal_real_part', -0.5, ...
            'reverse', 1, ...
            master_options, ...
            'save', 0);
    catch ME
        rethrow(ME)
    end
    
    steadyBranches_FB(i) = steadyBranchSTST;
    steadyBranches_FB(i).param = steadyBranchParam;
    steadyBranches_FB(i).nunst = nunst_steadyBranchSTST;
    steadyBranches_FB(i).indFold = ind_foldSteadyBranchSTST;
    steadyBranches_FB(i).indHopf = ind_hopfSteadyBranchSTST;
    steadyBranches_FB(i).sweep = steadyBranchSweep;
    
end

save([master_options.datadir_specific,'steadyBranches_FB'],'steadyBranches_FB')
