function [ steadyBranches ] = steadyBranchMultiCreator( ...
    funcs, bifparArray, ind_contin_param, branch_length, param, ...
    varargin)
%Creates multiple stst branches based on a list of input parameters.
%   Input:
%       funcs, ...
%       bifparArray, ...
%       ind_contin_param, ...
%       branch_length, ...
%       param
%
%   Options:
%       'sweep' = 0, 1
%           When sweep = 1 sweep up to the relevant current. strongDom
%           Otherwise do a turn on time series at current. weakDom
%       'step_bound_opt' = {'step', 0.1, ... 
%                           'max_step',[ind_contin_param, 0.1], ...
%                           'newton_max_iterations',10, ...
%                           'minimal_real_part', -0.5
%                           'max_bound',[ind_contin_param, 1],...
%                           'min_bound', [ind_contin_param, -1] };
%           Calling this will overwrite the default step_bound_opt. You can
%           also call each of these parameters by name to overwrite a
%           single parameter. If you call a cell array then that cell array
%           takes priority.
%       'minimal_real_part' = -0.5
%           For finding eigenvalues.
%
%   master_options: SAVING CURRENTLY NOT SUPPORTED.
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
%% options
p = inputParser;

% General option defaults
p.addParameter('step_bound_opt', 0)
p.addParameter('sweep',0)
p.addParameter('minimal_real_part',-0.5)

% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% first parse to set 'par' variable
p.KeepUnmatched = true;
parse(p,varargin{:})

% Second parse
% Create step_bound_opt, prepare rotational options    
if ~any(strcmp('step_bound_opt',p.UsingDefaults))
    % If the user input step_bound_opt
    step_bound_opt = p.Results.step_bound_opt;
else
    % If the user didn't input step_bound_opt
    % Create defaults for feed_phase
    if ind_contin_param == param.feed_phase.index
        p.addParameter('step',pi/64)
        p.addParameter('max_step',[ind_contin_param,2*pi/64])
        p.addParameter('newton_max_iterations',10)
        p.addParameter('max_bound',[ind_contin_param, 30*pi])
        p.addParameter('min_bound',[ind_contin_param, -30*pi])
        
        %{
        step_bound_opt_PHASE = { 'step',2*pi/64,...
            'max_step',[ind_feed_phase,(0.25)*2*pi/64], ...
            'newton_max_iterations',10, 'max_bound',[ind_feed_phase,6*pi] };
        %}

    % Create defaults for feed_amplitude
    elseif ind_contin_param == param.feed_ampli.index
        p.addParameter('step',0.003)
        p.addParameter('max_step',[ind_contin_param,0.003])
        p.addParameter('newton_max_iterations',10)
        p.addParameter('max_bound',[ind_contin_param,0.9999])
        p.addParameter('min_bound', [ind_contin_param,0.001])

        %{
        step_bound_opt_AMP = { 'step',0.003,'max_step',[ind_feed_ampli,0.003], ...
            'newton_max_iterations',10, ...
            'max_bound',[ind_feed_ampli,0.9999],...
            'min_bound', [ind_feed_ampli,0.001] };
        %}

    % Create defaults for all others
    else
        p.addParameter('step',0.01*par(ind_contin_param))
        p.addParameter('max_step',...
            [ind_contin_param,0.01*par(ind_contin_param)])
        p.addParameter('newton_max_iterations',10)
        p.addParameter('max_bound',...
            [ind_contin_param,2.0*par(ind_contin_param)])
        p.addParameter('min_bound',...
            [ind_contin_param,-2.0*par(ind_contin_param)])

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
    
    % Zip like python
    stacked = [step_bound_flags(:),step_bound_vals(:)].';
    zipped = stacked(:).';
    
    step_bound_opt = zipped;
end

% Make final options
options = p.Results;


%% setup

steadyBranches = struct('method',struct,'parameter',struct,'point',struct, ...
    'param', struct, 'nunst', 0, 'indFold', 0, 'indHopf', 0, 'timeDomain', struct);
steadyBranches = repmat(steadyBranches, [numel(bifparArray), 1]);


for i = 1:numel(bifparArray)
    % Temporary names and data
    steadyBranchParam = updateParams(param, 'feed_ampli', bifparArray(i));
    
    if options.sweep == 1;
        % Sweep Up
        steadyBranchNumerics = sweeper(steadyBranchParam.J.index, ...
            [60e-6, steadyBranchParam.J.value], steadyBranchParam, ...
            'save', 0);
        dde23_soln = steadyBranchNumerics(end).timeSeries;
    else
        % Turn On at the value
        dde23_soln = solver([1e-9;0;1e-9;0;0;0], ...
            [0,10], ...
            steadyBranchParam, ...
            'plot',0, ...
            'save',0);
    end
        

    
    try
        [steadyBranchSTST, nunst_steadyBranchSTST, ...
            ind_foldSteadyBranchSTST, ind_hopfSteadyBranchSTST] = ... 
            init_branch(funcs, ...
            dde23_soln.y(:,end), ...
            ind_contin_param, ...
            branch_length, ...
            steadyBranchParam, ...
            'minimal_real_part', options.minimal_real_part, ...
            'reverse', 1, ...
            'step_bound_opt', step_bound_opt, ...
            'save', 0);
    catch ME
        rethrow(ME)
    end
    
    steadyBranches(i).method = steadyBranchSTST.method;
    steadyBranches(i).parameter = steadyBranchSTST.parameter;
    steadyBranches(i).point = steadyBranchSTST.point;
    steadyBranches(i).param = steadyBranchParam;
    steadyBranches(i).nunst = nunst_steadyBranchSTST;
    steadyBranches(i).indFold = ind_foldSteadyBranchSTST;
    steadyBranches(i).indHopf = ind_hopfSteadyBranchSTST;
    steadyBranches(i).timeDomain = steadyBranchNumerics;
    
end




end

