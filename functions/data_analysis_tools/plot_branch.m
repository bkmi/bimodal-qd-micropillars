function [ ...
    ptBranInd_extrema, ...
    figHandle, ...
    lineHandles ] = plot_branch( branch, ...
    param, ...
    varargin )
%Plot a branch (by default) along its two free continuation parameters.
%   Considering we are using rotational symmetry stst will also work since
%   omega counts as a free parameter. The function plots two continuation 
%   parameters like this by default:
%   
%   branch.parameter.free(1:2).
%   The x-axis get branch.parameter.free(1).
%   The y-axis get branch.parameter.free(2).
%   
%   Input:
%       branch
%       param_struct
%       varargin
%
%   Output:
%       ptBranInd_extrema
%           ptBranInd_extrema is a struct containing the point index
%           along the branch of extrema in the x and y directions. It is
%           organized as follows:
%               ptBranInd_extrema.(xdir).max = 0
%               ptBranInd_extrema.(xdir).min = 0
%               ptBranInd_extrema.(ydir).max = 0
%               ptBranInd_extrema.(ydir).min = 0
%           Where (xdir) and (ydir) are the var_names given by param_struct
%
%   Options:
%       'axes_indParam' = [ 0, 0 ]
%           Calling this sets the axes along different parameters than
%           given by branch.parameter.free. The axes are determined by the
%           index of the parameter you enter. The x-axis goes with the
%           first and the y-axis goes with the second.
%       'add_2_gcf' = 0, 1
%           'add_2_gcf' = 1 -> Adds plot to current figure.
%           'add_2_gcf' = 0 or Not Called -> New figure.
%       'color' = [ 0, 0, 0 ] OR 'b', 'y', etc.
%           The branch will be the given color.
%       'nunst_color' = nunst_branch OR {nunst_branch, int_max}
%           Color the dots based on their nunst value. Overrides 'color'.
%           If you give a cell array with an int_max then the plotter will
%           plot as if int_max is the highest possible nunst. This is
%           useful if you want to plot multiple branches on the same plot.
%       'PlotStyle' = { 'LineStyle', '-', 'Marker', '.' }
%           Input a cell array which will be passed to the plotter. Usual
%           plot commands apply.
%       'polar' = 0, 1
%           With 1, the value for x turns into phi, the value for y turns
%           into rho. Polar can only handle options.color = '--r', etc...
%       'YTickLabel' = {'Off'}
%           List the names you want on the label.

%% Defaults + inputParser + Organize behavior

p = inputParser;

% General option defaults
p.addParameter('axes_indParam',...
    [branch.parameter.free(1),branch.parameter.free(2)])
p.addParameter('add_2_gcf', 0)
p.addParameter('color','b')
p.addParameter('nunst_color',[])
p.addParameter('PlotStyle', { 'LineStyle', 'none', 'Marker', '.' })
p.addParameter('polar',0)
p.addParameter('twoOmegaNunst',0)
p.addParameter('YTickLabel',{'Off'})


% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% Parse, set options
parse(p,varargin{:})
options = p.Results;


% Get hold status and save hold status
held_prior = ishold;

% line Handles
lineHandles = {};

% Handle add_2_gcf
if options.add_2_gcf==1
    % User called add_2_gcf
    figHandle = figure(gcf);
    hold on
elseif options.add_2_gcf==0
    % Default or user called add_2_gcf
    figHandle = figure;
    set(figHandle,'PaperType','a4')
    set(figHandle,'PaperOrientation','landscape');
    set(figHandle,'PaperUnits','normalized');
    set(figHandle,'PaperPosition', [0 0 1 1]);
    clf;
elseif isa(options.add_2_gcf,'struct')
    % This part is where we parse the object handle
    error('This is not supported yet.')
else
    error(['add_2_gcf can either equal 1, 0,', ... 
        'or a obj_struct. Otherwise, do not call it.'])
end

% Handle color settings, particularly nunst_color
if ~any(strcmp('nunst_color',p.UsingDefaults))
    % User called nunst_color
    if isa(options.nunst_color,'double')
        % Behavior based on giving a nunst_branch only
        if length(options.nunst_color) ~= length(branch.point)
            warning(['Len: nunst_branch =/= Len: points. ',...
                'Some points missing.'])
            % Calculate everything from nunst_branch in case the user 
            % wants to avoid plotting some points.
        end
    elseif isa(options.nunst_color,'cell')
        % Behavior is based on giving a cell array with nunst_branch and
        % int_max
        if length(options.nunst_color{1}) ~= length(branch.point)
            warning(['Len: nunst_branch =/= Len: points. ',...
                'Some points missing.'])
            % Calculate everything from nunst_branch in case the user 
            % wants to avoid plotting some points.
        end
    end
end


%% Extrema

if ~any(strcmp('twoOmegaNunst',p.UsingDefaults))
    % Prepare figure and vals
    x_param_vals = arrayfun(@(p)p.parameter(options.axes_indParam(1)), ... 
        branch.point); %Get fold continued parameter values for xdir
    y_param_vals1 = arrayfun(@(p)p.parameter(param.omega1.index), ...
        branch.point); %Get fold continued parameter values for ydir
    y_param_vals2 = arrayfun(@(p)p.parameter(param.omega2.index), ...
        branch.point); %Get fold continued parameter values for ydir
    
    y_param_vals = y_param_vals1;
else
    % Prepare figure and vals
    x_param_vals = arrayfun(@(p)p.parameter(options.axes_indParam(1)), ... 
        branch.point); %Get fold continued parameter values for xdir
    y_param_vals = arrayfun(@(p)p.parameter(options.axes_indParam(2)), ...
        branch.point); %Get fold continued parameter values for ydir
end

% Report extrema
ptBranInd_extrema = struct;
ptBranInd_extrema.(param.var_names{options.axes_indParam(1)}) = struct;
ptBranInd_extrema.(param.var_names{options.axes_indParam(2)}) = struct;
% x-dir
[max_val_x, max_ind_x] = max(x_param_vals);
[min_val_x, min_ind_x] = min(x_param_vals);
ptBranInd_extrema.(param.var_names{options.axes_indParam(1)}).Val = ... 
    [min_val_x, max_val_x];
ptBranInd_extrema.(param.var_names{options.axes_indParam(1)}).Ind = ... 
    [min_ind_x, max_ind_x];
% y-dir
[max_val_y, max_ind_y] = max(y_param_vals);
[min_val_y, min_ind_y] = min(y_param_vals);
ptBranInd_extrema.(param.var_names{options.axes_indParam(2)}).Val = ... 
    [min_val_y, max_val_y];
ptBranInd_extrema.(param.var_names{options.axes_indParam(2)}).Ind = ... 
    [min_ind_y, max_ind_y];


%% Plot


if options.polar ~= 1
    % For non-polar
    
    % Plot points
    if ~any(strcmp('twoOmegaNunst',p.UsingDefaults))
        
        % top plot, omega1
        subplot(1,2,1)
        if isa(options.twoOmegaNunst,'double')
            % Behavior based on giving a nunst_branch only
            nunst_pts = options.twoOmegaNunst;
            max_nunst = max(nunst_pts);
        elseif isa(options.twoOmegaNunst,'cell')
            % Behavior is based on giving a cell array with nunst_branch and
            % int_max
            nunst_pts = options.twoOmegaNunst{1};
            max_nunst = options.twoOmegaNunst{2};
        end
        
        sel=@(x,i)x(nunst_pts==i);
        colors = lines(max_nunst+1);

        hold on
        for i=0:max(nunst_pts)
            lineHandles{end+1} = plot(sel(x_param_vals,i),sel(y_param_vals1,i), ...
                'Color',colors(i+1,:), options.PlotStyle{:} );
        end
        hold off
        
        % Add title, axes
        title(strcat(param.plot_names(param.omega1.index), ...
            '-vs-', param.plot_names(options.axes_indParam(1)), ...
                  ' - Bifurcation Continuation'))
        xlabel([param.plot_names(options.axes_indParam(1)), ... 
            param.units(options.axes_indParam(1))])
        ylabel([param.plot_names(param.omega1.index), ...
            param.units(param.omega1.index)])
        
        %bottom plot, omega2
        subplot(1,2,2)
        if isa(options.twoOmegaNunst,'double')
            % Behavior based on giving a nunst_branch only
            nunst_pts = options.twoOmegaNunst;
            max_nunst = max(nunst_pts);
        elseif isa(options.twoOmegaNunst,'cell')
            % Behavior is based on giving a cell array with nunst_branch and
            % int_max
            nunst_pts = options.twoOmegaNunst{1};
            max_nunst = options.twoOmegaNunst{2};
        end
        
        sel=@(x,i)x(nunst_pts==i);
        colors = lines(max_nunst+1);

        hold on
        for i=0:max(nunst_pts)
            lineHandles{end+1} = plot(sel(x_param_vals,i),sel(y_param_vals2,i), ...
                'Color',colors(i+1,:), options.PlotStyle{:} );
        end
        hold off
        
        % Add title, axes
        title(strcat(param.plot_names(param.omega2.index), ...
            '-vs-', param.plot_names(options.axes_indParam(1)), ...
                  ' - Bifurcation Continuation'))
        xlabel([param.plot_names(options.axes_indParam(1)), ... 
            param.units(options.axes_indParam(1))])
        ylabel([param.plot_names(param.omega2.index), ...
            param.units(param.omega2.index)])
    
        colormap(colors)
        colorbar('YTickLabel','Off')
        
    elseif ~any(strcmp('nunst_color',p.UsingDefaults))

        if isa(options.nunst_color,'double')
            % Behavior based on giving a nunst_branch only
            nunst_pts = options.nunst_color;
            max_nunst = max(nunst_pts);
        elseif isa(options.nunst_color,'cell')
            % Behavior is based on giving a cell array with nunst_branch and
            % int_max
            nunst_pts = options.nunst_color{1};
            max_nunst = options.nunst_color{2};
        end


        sel=@(x,i)x(nunst_pts==i);
        colors = lines(max_nunst+1);

        hold on
        for i=0:max(nunst_pts)
            lineHandles{end+1} = plot(sel(x_param_vals,i),sel(y_param_vals,i), ...
                'Color',colors(i+1,:), options.PlotStyle{:} );
        end
        hold off

        %{
        for i=unique(nunst_pts)
            unique_nunst_vals = num2str(i);
        end
        %legend(unique_nunst_vals)
        %}

        colormap(colors)
        colorbar('YTickLabel',options.YTickLabel)
        %{
        colorbar('YTickLabel', ...
            [{(max_nunst/10)-1}, ...
            num2cell((max_nunst/10)-1 + max_nunst/10:max_nunst/10:max_nunst+1)])
        %}
    else
        lineHandles{end+1} = plot(x_param_vals,y_param_vals, ...
            'Color',options.color, options.PlotStyle{:} );
    end
elseif options.polar == 1
    % For polar
    
    polar(x_param_vals,y_param_vals, ...
        options.color );
else
    error('Polar == 1 for a polar plot, any other number is non-polar, you picked something else')
end

if any(strcmp('twoOmegaNunst',p.UsingDefaults))
    % WHEN twoOmegaNunst is NOT SELECTED
    
    % Add title, axes
    title(strcat(param.plot_names(options.axes_indParam(2)), ...
        '-vs-', param.plot_names(options.axes_indParam(1)), ...
              ' - Bifurcation Continuation'))
    xlabel([param.plot_names(options.axes_indParam(1)), ... 
        param.units(options.axes_indParam(1))])
    ylabel([param.plot_names(options.axes_indParam(2)), ...
        param.units(options.axes_indParam(2))])
end


% Return hold to what it was before running this function.
if held_prior==1
    hold on
elseif held_prior==0
    hold off
end

end

