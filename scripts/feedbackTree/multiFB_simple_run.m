%% load strongDom, tau = 0.83
treeLoader('datadir_specific', ...
    '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/zeroPhaseOffsetJ=560uA/strongDomAlpha=0/tau_fb=0.83ns/');

% CHECK OUT THIS FB VALUE:
% 

%% more stst sweeps
stdyBrnchFB_regSpace = steadyBranchMultiCreator( ...
    funcs,  ...
    [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], ...
    param.feed_phase.index, ...
    275, ...
    param, ...
    'sweep',1);

foldBranches = struct(...
    'method',struct, ...
    'parameter',struct, ...
    'point',struct, ...
    'newFuncs',@()'undefined',...
    'error', 0, ...
    'indError', NaN);
hopfBranches = struct(...
    'method',struct, ...
    'parameter',struct, ...
    'point',struct, ...
    'newFuncs',@()'undefined',...
    'error', 0, ...
    'indError', NaN);


%% 'interesting fb' - small and near fold intersections.

stdyBrnchFB_interesting = steadyBranchMultiCreator( ...
    funcs,  ...
    [0.025, 0.1330, 0.25, 0.3111, 0.3889], ...
    param.feed_phase.index, ...
    225, ...
    param, ...
    'sweep',1);

%% Probe fold

foldNear0133 = bifurFoldHopfMultiCreator( ...
    funcs, stdyBrnchFB_interesting(2), param, stdyBrnchFB_interesting(2).indFold, 100);

% foldTest(1) = br_contn(foldTest(1).newFuncs, foldTest(1), 50);
% foldTest(1) = br_rvers(foldTest(1));
% foldTest(1) = br_contn(foldTest(1).newFuncs, foldTest(1), 50);
% 
% foldTest(2) = br_contn(foldTest(2).newFuncs, foldTest(2), 50);
% foldTest(2) = br_rvers(foldTest(2));
% foldTest(2) = br_contn(foldTest(2).newFuncs, foldTest(2), 50);

%% Probe hopf

% % short
% hopfTest = bifurFoldHopfMultiCreator( ...
%      funcs, stdyBrnchFB_regSpace(2), param, stdyBrnchFB_regSpace(2).indHopf, 35);
% 
% % short
% hopfTest2 = bifurFoldHopfMultiCreator( ...
%     funcs, stdyBrnchFB_regSpace(3), param, stdyBrnchFB_regSpace(3).indHopf, 25);

hopfTestNear025 = bifurFoldHopfMultiCreator( ...
    funcs, stdyBrnchFB_interesting(3), param, stdyBrnchFB_interesting(3).indHopf, 60);

hopfTestNear0311 = bifurFoldHopfMultiCreator( ...
    funcs, stdyBrnchFB_interesting(4), param, stdyBrnchFB_interesting(4).indHopf, 20);

hopfTestNear0389 = bifurFoldHopfMultiCreator( ...
    funcs, stdyBrnchFB_interesting(5), param, stdyBrnchFB_interesting(5).indHopf, 20);



%% Branch plot

% branchplot = figure;
% 
% % Plot each testHopf
% for i = 1:numel(hopfTest)
%     % Add each hopf_branch
%     if isa(hopfTest(i).error,'double') && hopfTest(i).error == 0
%         % Only plot hopf_branches that DO NOT have errors
%         plot_branch(hopfTest(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% % Plot each testHopf
% for i = 1:numel(hopfTest2)
%     % Add each hopf_branch
%     if isa(hopfTest2(i).error,'double') && hopfTest2(i).error == 0
%         % Only plot hopf_branches that DO NOT have errors
%         plot_branch(hopfTest2(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% 
% % Plot each testFold
% for i = 1:numel(foldTest)
%     % Add each fold
%     if isa(foldTest(i).error,'double') && foldTest(i).error == 0
%         plot_branch(foldTest(i), param, ...
%             'add_2_gcf', 1, 'color','r');
%     end
% end
% 
% 
% % Plot stst branches
% for i = 1:numel(stdyBrnchFB_regSpace)
%     
%     plot_branch(stdyBrnchFB_regSpace(i), param, ...
%                 'add_2_gcf', 1, 'color','g', ...
%                 'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
% end






% branchplot = polar([0],[0]);
% 
% % Plot each testHopf
% for i = 1:numel(hopfTest)
%     % Add each hopf_branch
%     if isa(hopfTest(i).error,'double') && hopfTest(i).error == 0
%         % Only plot hopf_branches that DO NOT have errors
%         plot_branch(hopfTest(i), param, ...
%             'add_2_gcf', 1, 'color','c','polar',1);
%     end
% end
% 
% % Plot each testHopf
% for i = 1:numel(hopfTest2)
%     % Add each hopf_branch
%     if isa(hopfTest2(i).error,'double') && hopfTest2(i).error == 0
%         % Only plot hopf_branches that DO NOT have errors
%         plot_branch(hopfTest2(i), param, ...
%             'add_2_gcf', 1, 'color','c','polar',1);
%     end
% end
% 
% 
% % Plot each testFold
% for i = 1:numel(foldTest)
%     % Add each fold
%     if isa(foldTest(i).error,'double') && foldTest(i).error == 0
%         plot_branch(foldTest(i), param, ...
%             'add_2_gcf', 1, 'color','r','polar',1);
%     end
% end
% 
% 
% % Plot stst branches
% for i = 1:numel(stdyBrnchFB)
%     
%     plot_branch(stdyBrnchFB(i), param, ...
%                 'add_2_gcf', 1, 'color','.g', ...
%                 'axes_indParam',[param.feed_phase.index,param.feed_ampli.index], ...
%                 'polar',1);
% end









%% Add the tested bifurcations to the total list

% foldBranches(1) = foldTest(1);
% foldBranches(2) = foldTest(2);
% 
% hopfBranches(1) = hopfTest(2);


%% feed amp dir

stdyBranch_alongFB0_FA = pickAndSwitch(funcs, ...
    stdyBrnchFB_regSpace(1), ...
    param.feed_ampli.index, ...
    250, ...
    param, ...
    'fig', branchplot, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1);


%% SAVE

% % Sure you want to save and overwrite?
% saveSure = ...
%     input('\n\nAre you sure you want to save? \n0 = no \n1 = yes\n\n');
% % Force user to choose: OVERWRITE or not.
% while(1)
%     if saveSure == 1
%         saveit = 1; %Save and OVERWRITE
%         fprintf('Saving\n')
%         break
%     elseif saveSure == 0
%         saveit = 0; %Don't save and don't overwrite.
%         fprintf('Not saving\n')
%         break
%     end
% end
% 
% if saveSure == 1
%     disp('\nSaved\n')
%     save([master_options.datadir_specific,'stdyBrnchFB'],'stdyBrnchFB')
% 
%     save([master_options.datadir_specific,'foldBranches'],'foldBranches')
% 
%     save([master_options.datadir_specific,'hopfBranches'],'hopfBranches')
% end

%% For weakDom, tau = 0.83
% treeLoader('datadir_specific', ...
%     '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/zeroPhaseOffsetJ=560uA/weakDomAlpha=0/tau_fb=0.83ns/');

