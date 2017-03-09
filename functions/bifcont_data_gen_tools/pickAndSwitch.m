function [ outBranch ] = pickAndSwitch( ...
    funcs, branch, ind_contin_param, branch_length, param, ...
    varargin)
%Pick a point on a steady bifurcation branch and choose a new cont param.
%   Input:
%       funcs, ...
%       branch, ...
%       ind_contin_param, ...
%       branch_length, ...
%       param
%
%   Options:
%       'point'
%           Chose this option to avoid picking a point with the gui.
%           Instead it just uses this point number.
%       'fig' = fig_number
%           locate_along_branch will add/move the marker on the figure
%           given by fig_number.
%       'nunst_color' = nunst_branch OR {nunst_branch, int_max}
%           Color the dots based on their nunst value. Overrides 'color'.
%           If you give a cell array with an int_max then the plotter will
%           plot as if int_max is the highest possible nunst. This is
%           useful if you want to plot multiple branches on the same plot.
%       'axes_indParam' = [ 0, 0 ]
%           Calling this sets the axes along different parameters than
%           given by branch.parameter.free. The axes are determined by the
%           index of the parameter you enter. The x-axis goes with the
%           first and the y-axis goes with the second.
%
%   Options from init_branch:
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
%
%       'par_overwrite' = branch.point(1).parameter OR param.values
%           Calling this flag overrides the values given in the
%           param_struct. This is particularly useful if you would like to
%           make a new branch starting at a point along a DDEBIF
%           bifurcation continuation.
% 
%           THE ORDER/INDICES DETERMINED IN param_struct MUST BE THE SAME 
%           AS IN YOUR 'par_overwrite' VALUES OR THERE WILL BE MASSIVE 
%           ERROR.
%
%       'reverse' = 0, 1
%           Default to 0. 0 means the continuation will not reverse. 1
%           means the continuation will reverse.
%
%       'plot_prog' = 0, 1
%           Defaults to 1. 0 means don't plot progress, 1 means plot
%           progress.
%
%       'minimal_real_part' = -0.5
%           Defaults to -0.5. This is for getting stability.
%
%       'opt_inputs' = { 'extra_condition', [1], ...
%                        'print_residual_info',[0] }
%           By default, opt_inputs is set as above. YOU CANNOT CALL IT, I
%           was too lazy to write this. It can easily be written in if you
%           ever need to change it, but I have never needed to change it.
%% options
p = inputParser;
p.KeepUnmatched = true;

% General option defaults
p.addParameter('fig',gcf)
p.addParameter('axes_indParam',[0,0])
p.addParameter('nunst_color',0)
p.addParameter('point', 0)

% Make options
parse(p,varargin{:})
options = p.Results;


% For non-default fig, axes, and nunst values...
locateOpts = struct;

if ~any(strcmp('fig',p.UsingDefaults))
    % If fig IS NOT default
    locateOpts.fig = options.fig;
end

if ~any(strcmp('axes_indParam',p.UsingDefaults))
    % If axes_indParam IS NOT default
    locateOpts.axes_indParam = options.axes_indParam;
end

if ~any(strcmp('nunst_color',p.UsingDefaults))
    % If nunst_color IS NOT default
    locateOpts.nunst_color = options.nunst_color;
end


%% pick
if any(strcmp('point',p.UsingDefaults))
    % when point IS default
    [ptNum, ~] = locate_along_branch( branch, param, locateOpts );
else
    ptNum = options.point;
end


%% switch

% param = updateParams(param, 'par_overwrite', branch.point(ptNum).parameter);

[outBranch, nunstOutBranch, indFold, indHopf] = init_branch( funcs, ...
    branch.point(ptNum).x, ...
    ind_contin_param, ...
    branch_length, ...
    param, ...
    p.Unmatched, ...
    'save', 0, ...
    'par_overwrite', branch.point(ptNum).parameter);

outBranch.nunst = nunstOutBranch;
outBranch.indFold = indFold;
outBranch.indHopf = indHopf;
outBranch.timeDomain = struct;

end

