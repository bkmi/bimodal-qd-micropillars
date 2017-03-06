function [ out_branch, varargout ] = bifurContin_FoldHopf( ...
    funcs, ...
    init_branch, ...
    ind_bifurcation, ...
    ind_contin_param, ...
    branch_length, ...
    param_struct, ...
    varargin )
%This function will continue a fold or hopf bifurcation. Just provide the
%index of the bifurcation point and the directions to continue. The first
%parameter in ind_contin_param is the direction 'dir' continued.
%
%   Input:
%       funcs, ...
%       init_branch, ...
%       ind_bifurcation, ...
%       ind_contin_param, ...
%       branch_length, ...
%       param_struct, ...
%       varargin
%
%   Options:
%       'step_bound_opt' = {'step',0.1
%                           'max_step', ...
%                           [ind_contin_param(1),0.1, ...
%                           ind_contin_param(2),1]
%                           'newton_max_iterations',15
%                           'max_bound',...
%                           [ind_contin_param(1),10, ...
%                           ind_contin_param(2),10]
%                           'min_bound',...
%                           [ind_contin_param(1),-10, ...
%                           ind_contin_param(2),-10], ...
%                           'halting_accuracy', 1e-11, ...
%                           'minimal_accuracy', 1e-9 }
%           Calling this will overwrite the default step_bound_opt. You can
%           also call each of these parameters by name to overwrite a
%           single parameter. If you call a cell array then that cell array
%           takes priority.
%
%       'plot_prog' = 0, 1
%           Defaults to 1. 0 means don't plot progress, 1 means plot
%           progress.
%
%       'save_name' = 'branch_name'
%           bifurContin_FoldHopf will save the branch as 'branch_name' in 
%           a datadir_specific given by master_options. If you have choosen
%           'save' = 1 without a name given here AND this is not the first
%           time you've run solver then it will overwrite. Defaults to
%           'branch_FoldHopf'
%
%       'opt_inputs' = { 'extra_condition', [1], ...
%                        'print_residual_info',[0] }
%           By default, opt_inputs is set as above. YOU CANNOT CALL IT, I
%           was too lazy to write this. It can easily be written in if you
%           ever need to change it, but I have never needed to change it.
%
%   master_options:
%       'save' = 0, 1
%           By default, this is set to 0. When 'save' = 0, the function
%           does not try to save anything. When 'save' = 1, the function 
%           tries to save ________.
%       'datadir_specific' = '../data_bimodal-qd-micropillars/'
%           By default, this is set as above.
%       'dimensional' = 0, 1
%           By default, this is set to 0. When 'dimensional' = 0, the
%           function applies a non-dimensionalized system. When
%           'dimensional' = 1, the function applies a dimensionalized
%           system.

%% Defaults + InputParser + Organize Behavior

p = inputParser;

% General option defaults
p.addParameter('plot_prog', 1)
p.addParameter('step_bound_opt', 0)
p.addParameter('save_name', 'branch_FoldHopf')

% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% first parse
p.KeepUnmatched = true;
p.PartialMatching = false;
parse(p,varargin{:})

% Create step_bound_opt
if ~any(strcmp('step_bound_opt',p.UsingDefaults))
    % If the user input step_bound_opt
    step_bound_opt = p.Results.step_bound_opt;
else
    % If the user didn't input step_bound_opt
    
    if all( ind_contin_param == ...
            [ param_struct.feed_phase.index, ...
            param_struct.feed_ampli.index ] )
        % Create defaults for feed_phase, feed_ampli
        % Notice ind_contin_param(1) == feed_phase
        % and ind_contin_param(2) == feed_ampli
        p.addParameter('step',pi/64)
        p.addParameter('max_step', ...
            [ind_contin_param(1),pi/64, ...
            ind_contin_param(2),0.01])
        p.addParameter('newton_max_iterations',15)
        p.addParameter('max_bound',...
            [ind_contin_param(1),15*pi, ...
            ind_contin_param(2),2])
        p.addParameter('min_bound',...
            [ind_contin_param(1),-15*pi, ...
            ind_contin_param(2),-1])
        p.addParameter('halting_accuracy',1e-10)
        p.addParameter('minimal_accuracy',1e-8)

    elseif all( ind_contin_param == ...
            [ param_struct.feed_ampli.index, ...
            param_struct.feed_phase.index ] )
        % Create defaults for feed_ampli, feed_phase
        % Notice ind_contin_param(1) == feed_ampli
        % and ind_contin_param(2) == feed_phase
        p.addParameter('step',0.01)
        p.addParameter('max_step', ...
            [ind_contin_param(1),0.01, ...
            ind_contin_param(2),pi/64])
        p.addParameter('newton_max_iterations',15)
        p.addParameter('max_bound',...
            [ind_contin_param(1),2, ...
            ind_contin_param(2),15*pi])
        p.addParameter('min_bound',...
            [ind_contin_param(1),-2, ...
            ind_contin_param(2),-15*pi])
        p.addParameter('halting_accuracy',1e-10)
        p.addParameter('minimal_accuracy',1e-8)

    else
        % Create defaults for others
        p.addParameter('step',0.1)
        p.addParameter('max_step', ...
            [ind_contin_param(1),0.1, ...
            ind_contin_param(2),1])
        p.addParameter('newton_max_iterations',15)
        p.addParameter('max_bound',...
            [ind_contin_param(1),10, ...
            ind_contin_param(2),10])
        p.addParameter('min_bound',...
            [ind_contin_param(1),-10, ...
            ind_contin_param(2),-10])
        p.addParameter('halting_accuracy',1e-10)
        p.addParameter('minimal_accuracy',1e-8)

    end
    
    % Second parse to finalize options
    p.KeepUnmatched = false;
    parse(p,varargin{:})

    % Create step_bound_opt
    step_bound_flags = {};
    step_bound_vals  = {};

    if any(strcmp('step',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'step';
        step_bound_vals{end+1} = p.Results.step;
    end

    if any(strcmp('max_step',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'max_step';
        step_bound_vals{end+1} = p.Results.max_step;
    end

    if any(strcmp('newton_max_iterations',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'newton_max_iterations';
        step_bound_vals{end+1} = p.Results.newton_max_iterations;
    end

    if any(strcmp('max_bound',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'max_bound';
        step_bound_vals{end+1} = p.Results.max_bound;
    end

    if any(strcmp('min_bound',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'min_bound';
        step_bound_vals{end+1} = p.Results.min_bound;
    end
    
    if any(strcmp('halting_accuracy',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'halting_accuracy';
        step_bound_vals{end+1} = p.Results.halting_accuracy;
    end
    
    if any(strcmp('minimal_accuracy',fieldnames(p.Results)))
        step_bound_flags{end+1} = 'minimal_accuracy';
        step_bound_vals{end+1} = p.Results.minimal_accuracy;
    end
    
    % Zip like python
    stacked = [step_bound_flags(:),step_bound_vals(:)].';
    zipped = stacked(:).';
    
    step_bound_opt = zipped;
end

% Create rotational options, ind_contin_param
opt_inputs = {'extra_condition',1,'print_residual_info',0};
ind_contin_param_w_omega = [ind_contin_param, ...
    param_struct.omega1.index, ...
    param_struct.omega2.index];

% Get nunst_init_branch, determine which kind of bifurcation this is.
nunst_init_branch = GetRotStability(init_branch,funcs, 2);
ind_fold = find(abs(diff(nunst_init_branch))==1);
ind_hopf = find(abs(diff(nunst_init_branch))==2);
if any(ind_bifurcation == ind_fold)
    bifurcation_type = 'fold';
elseif any(ind_bifurcation == ind_hopf)
    bifurcation_type = 'hopf';
else
    error('ind_bifurcation is neither fold nor hopf.')
end


% Make final options
options = p.Results;

% Set save to 1 when the user called 'save_name'
if ~any(strcmp('save_name',p.UsingDefaults))
    options.save = 1;
end


% Prepare plotting
if options.plot_prog == 1
    figure
    out_branch.method.continuation.plot=1;
end


%% Fold Continuation // CALLS m_SetupRWFold in rot folder
if strcmp('fold',bifurcation_type)
    [foldfuncs,out_branch,~]=m_SetupRWFold( ...
        funcs, init_branch, ind_bifurcation,...
        'contpar',ind_contin_param_w_omega, ...
        'dir',ind_contin_param_w_omega(1), ...
        opt_inputs{:},...
        step_bound_opt{:}); 
    
    out_branch = br_contn(foldfuncs,out_branch,branch_length);
    out_branch = br_rvers(out_branch);
    out_branch = br_contn(foldfuncs,out_branch,branch_length);
    
    % Provide foldfuncs if the user gives two outputs AND it's a fold
    if nargout == 2 
        varargout{1} = foldfuncs;
    end
end


%% Hopf Continuation

if strcmp('hopf',bifurcation_type)
    [out_branch,~]=SetupRWHopf( ...
        funcs, init_branch, ind_bifurcation, ...
        'contpar',ind_contin_param_w_omega, ...
        'dir',ind_contin_param_w_omega(1), ...
        opt_inputs{:},...
        step_bound_opt{:});
    
    [out_branch,~,~,~] = br_contn(funcs,out_branch,branch_length);
    out_branch = br_rvers(out_branch);
    out_branch = br_contn(funcs,out_branch,branch_length);
    
    % Provide nothing if the user gives two outputs AND it's a hopf
    if nargout == 2 
        varargout{1} = @()'hopf';
    end
end


%% Save
% Save, if necessary
datadir_specific = options.datadir_specific;

if options.save == 1
    % Where will it save?
    fprintf(strcat('\n\n Saving in subfolder:\n', datadir_specific,'\n'))
end

if options.save == 1 && ...
        ~exist(strcat(datadir_specific,options.save_name,'.mat'),'file')
    save(strcat(datadir_specific,options.save_name),...
        'out_branch')
elseif options.save == 1 && ...
        exist(strcat(datadir_specific,options.save_name,'.mat'),'file')
    warning('That file already exists. Overwriting.')
    save(strcat(datadir_specific,options.save_name),...
        'out_branch')
end


end

