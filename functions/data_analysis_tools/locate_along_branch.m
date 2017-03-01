function [ branchpt, figHandle ] = locate_along_branch( branch, param, varargin )
%The function places a cursor at a point along the plot of a branch
%continuation. Once you choose 'y', the function returns the point number
%on that branch.
%
%   Input:
%       branch, ...
%       param
%
%   Output:
%       point_number
%
%   Options:
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

p = inputParser;

% General option defaults
p.addParameter('fig',gcf)
p.addParameter('axes_indParam',[0,0])
p.addParameter('nunst_color',0)

% Make options
parse(p,varargin{:})
options = p.Results;

% Select the correct paramter axes
if any(strcmp('axes_indParam',p.UsingDefaults))
    axes_indParam = branch.parameter.free;
else
    axes_indParam = options.axes_indParam;
end

% Select correct plot
if any(strcmp('fig',p.UsingDefaults))
    % No fig given
    if any(strcmp('nunst_color',p.UsingDefaults))
        % no nunst branch given
        [~,figHandle] = ...
            plot_branch(branch, param,'axes_indParam',axes_indParam);
    else
        % nunst branch given
        [~,figHandle] = ...
            plot_branch(branch, param,'axes_indParam',axes_indParam,...
            'nunst_color',options.nunst_color);
    end
else
    % Fig given
    figHandle = figure(options.fig);
    
    if any(strcmp('nunst_color',p.UsingDefaults))
        % no nunst branch given
        plot_branch(branch, param, ...
            'add_2_gcf', 1, 'color', 'b', ...
            'axes_indParam', axes_indParam);
    else
        % nunst branch given
        plot_branch(branch, param, ...
            'add_2_gcf', 1, 'nunst_color',options.nunst_color, ...
            'axes_indParam', axes_indParam);
    end
end




% Init the plotting
hold on

% Begin user control, create cursor default
curPoint = 1;

% init the cursor
x = branch.point(curPoint).parameter(axes_indParam(1));
y = branch.point(curPoint).parameter(axes_indParam(2));
cursor = plot(x, y, ...
    'linestyle', 'none', ...
    'Marker','o', ...
    'Color', 'r', ...
    'MarkerSize', 5);

while(1)
    % Update the cursor
    delete(cursor)
    x = branch.point(curPoint).parameter(axes_indParam(1));
    y = branch.point(curPoint).parameter(axes_indParam(2));
    cursor = plot(x, y, ...
        'linestyle', 'none', ...
        'Marker','o', ...
        'Color', 'r', ...
        'MarkerSize', 5);
    
    % User control
    disp(['Current point is ', num2str(curPoint), '.']);
    prompt = 'Modify pt value with a number. Select with "y" \n:';
    userInpt = input( prompt, 's' );
    disp('')
    
    % Update curPoint
    branchpt = curPoint; % For returning the value
    curPoint = curPoint + str2double(userInpt);
    
    if strcmp(userInpt,'y');
        % Break while loop if the user selects 'y'
        break
    end
    
end

% End of plotting
hold off


end

