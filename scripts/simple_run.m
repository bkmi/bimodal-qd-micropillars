%% Create initial branch and extend bifurcations from there. Basic setup.

clear;

feedPhaseMat = [1, 1; 1, 1];
feedAmpMat = [1, 0; 0, 1];

% Setup parameters, save them?
% setup_params_nonDim_CnstCplRatio('save',0,'feed_ampli',0.373, 'alpha_par',0,'clear',0)
% setup_params_nonDim_CnstCplRatio('save',0,'feed_ampli',0.1, 'alpha_par',0,'clear',0)
% setup_params_nonDim_CnstCplRatio('save',1,'feed_ampli',0.373, 'feed_phase',0,'alpha_par',0,'clear',0)

% For feedback phase continuation
setup_params_nonDim_CnstCplRatio(...
    'save',1, ...
    'alpha_par',0, ...
    'feed_ampli',0.1, ...
    'feed_ampliMatrix', feedAmpMat, ...
    'feed_phase',0, ... 'feed_phaseMatrix', feedPhaseMat, ...
    'clear',0)

% For injection current continuation
% setup_params_nonDim_CnstCplRatio(...
%     'save',1, ...
%     'alpha_par',0, ...
%     'J', 90e-6, ...
%     'feed_ampli',0, ...
%     'feed_ampliMatrix', feedAmpMat, ...
%     'feed_phase',0, ... 'feed_phaseMatrix', feedPhaseMat, ...
%     'clear',0)

% For beta continuation
% setup_params_nonDim_CnstCplRatio(...
%     'save',0, ...
%     'alpha_par',0, ...
%     'J', 90e-6, ...
%     'beta',0.0056,...
%     'feed_ampli',0.01, ...
%     'feed_ampliMatrix', feedAmpMat, ...
%     'feed_phase',0, ... 'feed_phaseMatrix', feedPhaseMat, ...
%     'clear',0)




% Create and save turn on time series
% 'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot)
dde23_soln = solver([1e-9;0;1e-9;0;0;0], ...
    [0,10], ...
    param, ...
    master_options, ...
    'plot',1, ...
    'dde23_options',ddeset('RelTol',10^-8,'OutputFcn', @odeplot) );

%% Create initial branch
% For feedback phase continuation
[branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
    init_branch(funcs, ...
    dde23_soln.y(:,end), param.feed_phase.index, 250, param, ...
    'max_step',[param.feed_phase.index, (0.8)*pi/32], ...
    'minimal_real_part', -1, ...
    'reverse',1, ...
    master_options);

% For alpha continuaiton
% [branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
%     init_branch(funcs, ...
%     dde23_soln.y(:,end), param.alpha_par.index, 200, param, ...
%     'minimal_real_part', -1, master_options, ...
%     'step_bound_opt', ...
%     {'step', 1e-8, ... 
%     'max_step',[param.alpha_par.index, 1e-8], ...
%     'newton_max_iterations',20, ...
%     'max_bound',[param.alpha_par.index, 2],...
%     'min_bound', [param.alpha_par.index, -1] ...
%     },...
%     master_options);

% For feedback amplitude continuation
% [branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
%     init_branch(funcs, ...
%     dde23_soln.y(:,end), param.feed_ampli.index, 50, param, ...
%     'max_step',[param.feed_ampli.index, 0.01], 'reverse', 1, ...
%     master_options);

% For injection current continuation
% [branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
%     init_branch(funcs, ...
%     dde23_soln.y(:,end), param.J.index, 200, param, ...
%     'step_bound_opt', ...
%     {'step', 1e-6, ... 
%     'max_step',[param.J.index, 1e-6], ...
%     'newton_max_iterations',15, ...
%     'max_bound',[param.J.index, 250e-6],...
%     'min_bound', [param.J.index, 0e-6] ...
%     },...
%     master_options);

% For beta continuation
% [branch_stst, nunst_branch_stst, ind_fold, ind_hopf] = ... 
%     init_branch(funcs, ...
%     dde23_soln.y(:,end), param.beta.index, 200, param, ...
%     'step_bound_opt', ...
%     {'step', 0.05e-3, ... 
%     'max_step',[param.beta.index, 0.05e-3], ...
%     'newton_max_iterations',15, ...
%     'max_bound',[param.beta.index, 14e-3],...
%     'min_bound', [param.beta.index, 1e-3] ...
%     },...
%     master_options);

% [TEST_nunst,~,~,pts] = GetStability(branch_stst, ...
%     'exclude_trivial',true, ...
%     'locate_trivial',@(p)[0,0],... % Two rotating waves
%     'funcs',funcs,...
%     'recompute',true);

%% Create structs for fold_branches
% Fold
fold_branches = struct;
for i = 1:length(ind_fold)
    fold_active_branch_name = ...
        strcat('f',num2str(i),'branch');

    try
        fbranch = ...
            bifurContin_FoldHopf( ...
            funcs, ... 
            branch_stst, ...
            ind_fold(i), ...
            [param.feed_phase.index, param.feed_ampli.index], ...
            20, ...
            param,...
            'plot_prog', 1, ...
            master_options,...
            'save',0);

        fold_branches.(fold_active_branch_name) = fbranch;
    catch ME
        switch ME.identifier
            case 'br_contn:start'
                warning(ME.message);
                warning(strcat('During branch=',fold_active_branch_name));
                fold_branches.(fold_active_branch_name).error = ME;
                fold_branches.(fold_active_branch_name).fold_active_ind = ...
                    i;
                fold_branches.(fold_active_branch_name).fold_active_branch_name = ...
                    fold_active_branch_name;
            otherwise
                rethrow(ME)
        end
    end
end


%% Create structs for hopf_branches
% Hopf
hopf_branches = struct;
for i = 1:length(ind_hopf)
    hopf_active_branch_name = ...
        strcat('h',num2str(i),'branch');

    try
        hbranch = ...
            bifurContin_FoldHopf( ...
            funcs, ... 
            branch_stst, ...
            ind_hopf(i), ...
            [param.feed_phase.index, param.feed_ampli.index], ...
            20, ...
            param,...
            'plot_prog', 1, ...
            master_options,...
            'save',0);

        hopf_branches.(hopf_active_branch_name) = hbranch;
    catch ME
        switch ME.identifier
            case 'br_contn:start'
                warning(ME.message);
                warning(strcat('During branch=',hopf_active_branch_name));
                hopf_branches.(hopf_active_branch_name).error = ME;
                hopf_branches.(hopf_active_branch_name).hopf_active_ind = ...
                    i;
                hopf_branches.(hopf_active_branch_name).hopf_active_branch_name = ...
                    hopf_active_branch_name;
            otherwise
                rethrow(ME)
        end
    end
end

% Save hopf and fold branches
save([master_options.datadir_specific,'hopf_branches'],'hopf_branches')
save([master_options.datadir_specific,'fold_branches'],'fold_branches')
