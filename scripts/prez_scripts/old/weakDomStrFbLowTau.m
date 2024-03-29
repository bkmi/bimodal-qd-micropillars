%% weak field dom, strong fb 

clear;

feedPhaseMat = [1, 0; 0, 0];
feedAmpMat = [1, 0; 0, 0];


% For feedback phase continuation
setup_params_nonDim_CnstCplRatio(...
    'save',1, ...
    'alpha_par',0, ...
    'feed_ampli',0.1, ...
    'feed_ampliMatrix', feedAmpMat, ...
    'feed_phase',0, ...
    'feed_phaseMatrix', feedPhaseMat, ...
    'clear',0,...
    'J',560e-6,...
    'tau_fb', 0.8, ...
    'datadir_specific', ...
    strcat(data_directory(), '/weakDomStrFbLowTau'))

% sweepSoln = sweeper(param.J.index, [60e-6, param.J.value], param, master_options);
% dde23_soln = sweepSoln(end).timeSeries;

dde23_soln = solver([1e-9;0;1e-9;0;0;0], ...
    [0,10], ...
    param, ...
    master_options, ...
    'plot',1);

step_bound_forHOPF = { ...
    'step',pi/64, ...
    'max_step', ...
    [param.feed_phase.index,(1) * pi/64, ...
    param.feed_ampli.index, (1) * 0.01], ...
    'newton_max_iterations',20, ...
    'max_bound',[param.feed_phase.index,15*pi, param.feed_ampli.index,2], ...
    'min_bound',[param.feed_phase.index,-15*pi, param.feed_ampli.index,-1], ...
    'halting_accuracy',1e-10, ...
    'minimal_accuracy',1e-6 };

step_bound_forFOLD = { ...
    'step',pi/64, ...
    'max_step', ...
    [param.feed_phase.index,(1) * pi/32, ...
    param.feed_ampli.index, (1) * 0.05], ...
    'newton_max_iterations',20, ...
    'max_bound',[param.feed_phase.index,15*pi, param.feed_ampli.index,2], ...
    'min_bound',[param.feed_phase.index,-15*pi, param.feed_ampli.index,-1], ...
    'halting_accuracy',1e-10, ...
    'minimal_accuracy',1e-6 };

save([master_options.datadir_specific,'step_bound_opts.mat'], ...
    'step_bound_forFOLD', 'step_bound_forHOPF');


%%


% opt_inputs = { 'extra_condition', [1],'print_residual_info',[1] };
[branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
    init_branch(funcs, ...
    dde23_soln.y(:,end), param.feed_phase.index, 500, param, ...
    'max_step',[param.feed_phase.index, (1)*pi/32], ...
    'minimal_real_part', -1, ...
    'halting_accuracy', 1e-12, ...
    'minimal_accuracy', 1e-10, ...
    'newton_max_iterations', 15, ...
    'reverse',1, ...
    master_options);
% , ...
%     'opt_inputs', opt_inputs);


%% out of phase feedback
vertBranch00 = pickAndSwitch(funcs, ...
    branch_stst, ...
    param.feed_ampli.index, ...
    250, ...
    param, ...
    'reverse',1, 'nunst_color', nunst_branch_stst, ...
    'point', 421);

vertBranch00 = bif_contnAndStab(funcs,vertBranch00,300,'reverse',1);

plot_branch(vertBranch00, param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'nunst_color', vertBranch00.nunst)

save([master_options.datadir_specific,'vertBranch25.mat'],'vertBranch25');

%% from out of phase feedback
horizBranch006 = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    350, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1, ...
    'point', 100); % directly next to fold, low fa

plot_branch(horizBranch006, param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'nunst_color', horizBranch006.nunst)

save([master_options.datadir_specific,'horizBranch006.mat'],'horizBranch006');

%% verts from horizBranch006 in stable and unstable region.
%These are essentially more comprehensive vertBranch25. I wouldn't use them

horizBranch006.method.point.halting_accuracy = 1e-8;
horizBranch006.method.point.minimal_accuracy = 1e-6;

%stable
vertBranch3Stab = pickAndSwitch(funcs, ...
    horizBranch006, ...
    param.feed_ampli.index, ...
    70, ...
    param, ...
    'reverse',1, ...
    'nunst_color', horizBranch006.nunst, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'point', 321);

plot_branch(vertBranch3Stab, param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'nunst_color', vertBranch3Stab.nunst)

save([master_options.datadir_specific,'vertBranch3Stab.mat'],'vertBranch3Stab');

% unstable
vertBranch3unst = pickAndSwitch(funcs, ...
    horizBranch006, ...
    param.feed_ampli.index, ...
    350, ...
    param, ...
    'reverse',1, ...
    'nunst_color', horizBranch006.nunst, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'point', 339);

plot_branch(vertBranch3unst, param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'nunst_color', vertBranch3unst.nunst)

save([master_options.datadir_specific,'vertBranch3unst.mat'],'vertBranch3unst');


%% folds from horizBranch015

folds_horBrn006 = bifurFoldHopfMultiCreator( ...
    funcs, ...
    horizBranch006, ...
    param, ...
    horizBranch006.indFold, ...
    70, ...
    'step_bound_opt', step_bound_forFOLD);

% extending
for i = 1:numel(folds_horBrn006)
    if mod(i,2) == 0
        %even
        folds_horBrn006(i) = br_contn(folds_horBrn006(i).newFuncs, ...
            folds_horBrn006(i), 100);
    else
        %odd
        folds_horBrn006(i) = br_rvers(folds_horBrn006(i));
        folds_horBrn006(i) = br_contn(folds_horBrn006(i).newFuncs, ...
            folds_horBrn006(i), 100);
    end
end

branchplot = figure;

% plotting
for i = 1:numel(folds_horBrn006)
    % Add each fold
    if isa(folds_horBrn006(i).error,'double') && folds_horBrn006(i).error == 0
        plot_branch(folds_horBrn006(i), param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

save([master_options.datadir_specific,'folds_horBrn006.mat'],'folds_horBrn006');

%% Let's go straight up the feedback phase = 0
vertBranch00 = pickAndSwitch(funcs, ...
    branch_stst, ...
    param.feed_ampli.index, ...
    300, ...
    param, ...
    'reverse',1, 'nunst_color', nunst_branch_stst, ...
    'point', 496);

plot_branch(vertBranch00, param, ...
    'nunst_color', vertBranch00.nunst)

save([master_options.datadir_specific,'vertBranch00.mat'],'vertBranch00');


%% unstable folds + hopf

folds_vrtBrn00 = bifurFoldHopfMultiCreator( ...
    funcs, ...
    vertBranch00, ...
    param, ...
    vertBranch00.indFold, ...
    300, ...
    'step_bound_opt', step_bound_forFOLD);

branchplot = figure;

% plotting
for i = 1:numel(folds_vrtBrn00)
    % Add each fold
    if isa(folds_vrtBrn00(i).error,'double') && folds_vrtBrn00(i).error == 0
        plot_branch(folds_vrtBrn00(i), param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

save([master_options.datadir_specific,'folds_vrtBrn00.mat'],'folds_vrtBrn00');

% hopfs
hopfs_vrtBrn00 = bifurFoldHopfMultiCreator( ...
    funcs, ...
    vertBranch00, ...
    param, ...
    vertBranch00.indHopf, ...
    300, ...
    'step_bound_opt', step_bound_forHOPF);

branchplot = figure;

for i = 1:numel(hopfs_vrtBrn00)
    % Add each hopf
    if isa(hopfs_vrtBrn00(i).error,'double') && hopfs_vrtBrn00(i).error == 0
        plot_branch(hopfs_vrtBrn00(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

save([master_options.datadir_specific,'hopfs_vrtBrn00.mat'],'hopfs_vrtBrn00');


%% Horiz right below U fold heights

horizBranch035Stabl = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    350, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1, ...
    'point', 220); % directly next to fold, low fa

plot_branch(horizBranch035Stabl, param, ...
    'nunst_color', horizBranch035Stabl.nunst)

save([master_options.datadir_specific,'horizBranch035Stabl.mat'], ...
    'horizBranch035Stabl');

% There are only unstable hopfs birthed here
horizBranch035unstab = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    75, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1, ...
    'point', 213); % directly next to fold, low fa

plot_branch(horizBranch035unstab, param, ...
    'nunst_color', horizBranch035unstab.nunst)

save([master_options.datadir_specific,'horizBranch035unstab.mat'], ...
    'horizBranch035unstab');


%% more horiz below U shape
horizBranch025Stabl = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    250, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1,  ...
    'nunst_color', vertBranch00.nunst, ...
    'max_step', [23, 0.0982, param.omega1.index 0.2], ...
    'point', 252); % directly next to fold, low fa

plot_branch(horizBranch025Stabl, param, ...
    'nunst_color', horizBranch025Stabl.nunst)

save([master_options.datadir_specific,'horizBranch025Stabl.mat'], ...
    'horizBranch025Stabl');



horizBranch02Stabl = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    250, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1,  ...
    'nunst_color', vertBranch00.nunst, ...
    'max_step', [23, 0.0982, param.omega1.index 0.2], ...
    'point', 269); % directly next to fold, low fa

plot_branch(horizBranch02Stabl, param, ...
    'nunst_color', horizBranch02Stabl.nunst)

save([master_options.datadir_specific,'horizBranch02Stabl.mat'], ...
    'horizBranch02Stabl');


horizBranch03Stabl = pickAndSwitch(funcs, ...
    vertBranch00, ...
    param.feed_phase.index, ...
    250, ...
    param, ...
    'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
    'reverse',1,  ...
    'nunst_color', vertBranch00.nunst, ...
    'max_step', [23, 0.0982, param.omega1.index 0.2], ...
    'point', 236); % directly next to fold, low fa

plot_branch(horizBranch03Stabl, param, ...
    'nunst_color', horizBranch03Stabl.nunst)

save([master_options.datadir_specific,'horizBranch03Stabl.mat'], ...
    'horizBranch03Stabl');


%% hopfs from horizBranch025Stabl, horizBranch03Stabl
hopfs_horBrn025Stabl = bifurFoldHopfMultiCreator( ...
    funcs, ...
    horizBranch025Stabl, ...
    param, ...
    horizBranch025Stabl.indHopf, ...
    60, ...
    'step_bound_opt', step_bound_forHOPF);

branchplot = figure;

for i = 1:numel(hopfs_horBrn025Stabl)
    % Add each hopf
    if isa(hopfs_horBrn025Stabl(i).error,'double') && hopfs_horBrn025Stabl(i).error == 0
        plot_branch(hopfs_horBrn025Stabl(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

save([master_options.datadir_specific,'hopfs_horBrn025Stabl.mat'], ...
    'hopfs_horBrn025Stabl');


hopfs_horBrn03Stabl = bifurFoldHopfMultiCreator( ...
    funcs, ...
    horizBranch03Stabl, ...
    param, ...
    horizBranch03Stabl.indHopf, ...
    60, ...
    'step_bound_opt', step_bound_forHOPF);

branchplot = figure;

for i = 1:numel(hopfs_horBrn03Stabl)
    % Add each hopf
    if isa(hopfs_horBrn03Stabl(i).error,'double') && hopfs_horBrn03Stabl(i).error == 0
        plot_branch(hopfs_horBrn03Stabl(i), param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

save([master_options.datadir_specific,'hopfs_horBrn03Stabl.mat'], ...
    'hopfs_horBrn03Stabl');



















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












% %% dde for a slightly higher feedback, but in-phase feedback
% 
% param02FB = updateParams(param, 'feed_ampli', 0.2);
% 
% dde23_02FB = solver([1e-9;0;1e-9;0;0;0], ...
%     [0,10], ...
%     param02FB, ...
%     master_options, ...
%     'plot',1, ...
%     'save', 0);
% 
% save([master_options.datadir_specific,'solns02FB.mat'], ...
%     'param02FB', dde23_02FB);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %% No hopf, only fold
% horizBranch031 = pickAndSwitch(funcs, ...
%     vertBranch00, ...
%     param.feed_phase.index, ...
%     350, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 246);
% 
% plot_branch(horizBranch031, param, ...
%     'nunst_color', horizBranch031.nunst)
% 
% %% some hopf, indHopf(2:3) are from stable to unstable, lines on side hopfs
% horizBranch04 = pickAndSwitch(funcs, ...
%     vertBranch00, ...
%     param.feed_phase.index, ...
%     200, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 216);
% 
% plot_branch(horizBranch04, param, ...
%     'nunst_color', horizBranch04.nunst)
% 
% save([master_options.datadir_specific,'horizBranch04.mat'],'horizBranch04');
% 
% 
% %% some hopf, lines on side hopfs
% horizBranch045 = pickAndSwitch(funcs, ...
%     vertBranch00, ...
%     param.feed_phase.index, ...
%     200, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 198);
% 
% plot_branch(horizBranch045, param, ...
%     'nunst_color', horizBranch045.nunst)
% 
% save([master_options.datadir_specific,'horizBranch045.mat'],'horizBranch045');
% 
% 
% %% hopf from horizBranch045, lines on side hopfs
% 
% hopfs_horBrn045 = bifurFoldHopfMultiCreator( ...
%     funcs, ...
%     horizBranch045, ...
%     param, ...
%     horizBranch045.indHopf, ...
%     120, ...
%     'step_bound_opt', step_bound_forHOPF);
% 
% hopfs_horBrn045(2) = br_contn(funcs,hopfs_horBrn045(2), 250);
% hopfs_horBrn045(3) = br_rvers(hopfs_horBrn045(3));
% hopfs_horBrn045(3) = br_contn(funcs,hopfs_horBrn045(3), 250);
% 
% branchplot = figure;
% 
% for i = 1:numel(hopfs_horBrn045)
%     % Add each hopf
%     if isa(hopfs_horBrn045(i).error,'double') && hopfs_horBrn045(i).error == 0
%         plot_branch(hopfs_horBrn045(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% save([master_options.datadir_specific,'hopfs_horBrn045.mat'],'hopfs_horBrn045');
% 
% 
% %% some hopf
% horizBranch05 = pickAndSwitch(funcs, ...
%     vertBranch00, ...
%     param.feed_phase.index, ...
%     200, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 181);
% 
% plot_branch(horizBranch05, param, ...
%     'nunst_color', horizBranch05.nunst)
% 
% save([master_options.datadir_specific,'horizBranch05.mat'],'horizBranch05');
% 
% %% hopf from horizBranch05, lines on side hopfs
% 
% hopfs_horBrn05 = bifurFoldHopfMultiCreator( ...
%     funcs, ...
%     horizBranch05, ...
%     param, ...
%     horizBranch05.indHopf, ...
%     60);
% 
% 
% branchplot = figure;
% 
% for i = 1:numel(hopfs_horBrn05)
%     % Add each hopf
%     if isa(hopfs_horBrn05(i).error,'double') && hopfs_horBrn05(i).error == 0
%         plot_branch(hopfs_horBrn05(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% save([master_options.datadir_specific,'hopfs_horBrn05.mat'],'hopfs_horBrn05');
% 
% %% try and find hopfs vertically
% % Maybe try going above 1??
% 
% vertBranch31 = pickAndSwitch(funcs, ...
%     horizBranch05, ...
%     param.feed_ampli.index, ...
%     200, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 175);
% 
% plot_branch(vertBranch31, param, ...
%     'nunst_color', vertBranch31.nunst)
% 
% save([master_options.datadir_specific,'vertBranch31.mat'],'vertBranch31');
% 
% %% THIS ONE GOES PLACES! LARGE FEEDPHASE VERT CUT
% vertBranch905 = pickAndSwitch(funcs, ...
%     horizBranch05, ...
%     param.feed_ampli.index, ...
%     350, ...
%     param, ...
%     'axes_indParam', [param.feed_phase.index, param.feed_ampli.index], ...
%     'reverse',1, ...
%     'point', 114);
% 
% plot_branch(vertBranch905, param, ...
%     'nunst_color', vertBranch905.nunst, 'twoOmegaNunst', vertBranch905.nunst)
% 
% save([master_options.datadir_specific,'vertBranch905.mat'],'vertBranch905');
% 
% 
% %% long slopeing lines
% 
% hopfs_vrtBrn905 = bifurFoldHopfMultiCreator( ...
%     funcs, ...
%     vertBranch905, ...
%     param, ...
%     vertBranch905.indHopf([4,7:end]), ... 
%     100, ...
%     'step_bound_opt', step_bound_forHOPF);
% 
% hopfs_vrtBrn905(3) = br_contn(funcs, hopfs_vrtBrn905(3), 150);
% hopfs_vrtBrn905(4) = br_contn(funcs, hopfs_vrtBrn905(4), 75);
% 
% branchplot = figure;
% 
% for i = 1:numel(hopfs_vrtBrn905)
%     % Add each hopf
%     if isa(hopfs_vrtBrn905(i).error,'double') && hopfs_vrtBrn905(i).error == 0
%         plot_branch(hopfs_vrtBrn905(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% save([master_options.datadir_specific,'hopfs_vrtBrn905.mat'],'hopfs_vrtBrn905');
% 
% 
% %% feed amp cont
% horizBranch04Unst = pickAndSwitch(funcs, ...
%     vertBranch905, ...
%     param.feed_phase.index, ...
%     400, ...
%     param, ...
%     'nunst_color', vertBranch905.nunst, ...
%     'reverse',1, ...
%     'point', 312);
% 
% horizBranch025Unst = pickAndSwitch(funcs, ...
%     vertBranch905, ...
%     param.feed_phase.index, ...
%     400, ...
%     param, ...
%     'nunst_color', vertBranch905.nunst, ...
%     'reverse',1, ...
%     'point', 368);
% 
% 
% %% hopf from the weird ones above
% 
% hopfs_horBrn04Unst = bifurFoldHopfMultiCreator( ...
%     funcs, ...
%     horizBranch04Unst, ...
%     param, ...
%     horizBranch04Unst.indHopf, ... 
%     100, ...
%     'step_bound_opt', step_bound_forHOPF);
% 
% hopfs_horBrn04Unst(2) = br_rvers(hopfs_horBrn04Unst(2));
% hopfs_horBrn04Unst(2) = br_contn(funcs,hopfs_horBrn04Unst(2), 100);
% hopfs_horBrn04Unst(6) = br_contn(funcs, hopfs_horBrn04Unst(6), 100);
% hopfs_horBrn04Unst(7) = br_contn(funcs, hopfs_horBrn04Unst(7), 100);
% hopfs_horBrn04Unst(8) = br_contn(funcs, hopfs_horBrn04Unst(8), 100);
% 
% branchplot = figure;
% 
% for i = 1:numel(hopfs_horBrn04Unst)
%     % Add each hopf
%     if isa(hopfs_horBrn04Unst(i).error,'double') && hopfs_horBrn04Unst(i).error == 0
%         plot_branch(hopfs_horBrn04Unst(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% save([master_options.datadir_specific,'hopfs_horBrn04Unst.mat'],'hopfs_horBrn04Unst');
% 
% 
% %% Like above "THIS ONES GOES PLACES" but on the other side!
% vertBranchNeg15 = pickAndSwitch(funcs, ...
%     horizBranch05, ...
%     param.feed_ampli.index, ...
%     50, ...
%     param, ...
%     'nunst_color', horizBranch05.nunst, ...
%     'reverse',1, ...
%     'point', 366);
% 
% % vertBranchNeg15 = br_rvers(vertBranchNeg15);
% % vertBranchNeg15 = br_contn(funcs, vertBranchNeg15, 30);
% % [vertBranchNeg15.nunst,~,~,vertBranchNeg15.point] =  ...
% %     GetRotStability(vertBranchNeg15, funcs, 2);
% 
% plot_branch(vertBranchNeg15, param, ...
%     'nunst_color', vertBranchNeg15.nunst, 'twoOmegaNunst', vertBranchNeg15.nunst)
% 
% save([master_options.datadir_specific,'vertBranchNeg15.mat'],'vertBranchNeg15');
% 
% 
% %% hopf, fold plot
% 
% branchplot = figure;
% 
% for i = [3,4,5,6] % numel(folds_horBrn015)
%     % Add each fold
%     if isa(folds_horBrn006(i).error,'double') && folds_horBrn006(i).error == 0
%         plot_branch(folds_horBrn006(i), param, ...
%             'add_2_gcf', 1, 'color','r');
%     end
% end
% 
% % TOO OUT OF THE WAY
% % for i = [3] % numel(hopfs_horBrn05)
% %     % Add each hopf
% %     if isa(hopfs_horBrn05(i).error,'double') && hopfs_horBrn05(i).error == 0
% %         plot_branch(hopfs_horBrn05(i), param, ...
% %             'add_2_gcf', 1, 'color','c');
% %     end
% % end
% 
% for i = [1,2,3,4 ] % 1:numel(hopfs_horBrn045)
%     % Add each hopf
%     if isa(hopfs_horBrn045(i).error,'double') && hopfs_horBrn045(i).error == 0
%         plot_branch(hopfs_horBrn045(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% for i = 1:numel(hopfs_vrtBrn905)
%     % Add each hopf
%     if isa(hopfs_vrtBrn905(i).error,'double') && hopfs_vrtBrn905(i).error == 0
%         plot_branch(hopfs_vrtBrn905(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% for i = 1:numel(hopfs_horBrn04Unst)
%     % Add each hopf
%     if isa(hopfs_horBrn04Unst(i).error,'double') && hopfs_horBrn04Unst(i).error == 0
%         plot_branch(hopfs_horBrn04Unst(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end
% 
% for i = 1
%     % Add each hopf
%     if isa(leftHighHopf(i).error,'double') && leftHighHopf(i).error == 0
%         plot_branch(leftHighHopf(i), param, ...
%             'add_2_gcf', 1, 'color','c');
%     end
% end