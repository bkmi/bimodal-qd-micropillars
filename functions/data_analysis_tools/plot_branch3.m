function [  ] = plot_branch3( branch, param, varargin )
%Plot a branch along its first three free continuation parameters.
%   The x-axis gets branch.parameter.free(1).
%   The y-axis gets branch.parameter.free(2).
%   The z-axis gets branch.parameter.free(3).
%% Defaults + inputParser + Organize behavior

p = inputParser;

% General option defaults

% What param/value goes on each axis
if numel(branch.parameter.free) >= 3
    % When there are at least three free parameters. Get the first three.
    p.addParameter('axes_indParam',...
        {branch.parameter.free(1), ...
        branch.parameter.free(2), ...
        branch.parameter.free(3)} )
else 
    % Put strong field intensity on the z axis
    p.addParameter('axes_indParam',...
        {branch.parameter.free(1), ...
        branch.parameter.free(2), ...
        'x1'} )
end

p.addParameter('add_2_gcf', 0)

p.addParameter('color','b')
p.addParameter('nunst_color',[])
p.addParameter('PlotStyle', { 'LineStyle', 'none', 'Marker', '.' })

% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% Parse, set options
parse(p,varargin{:})
options = p.Results;

% Save hold status
held_prior = ishold;

% Handle add_2_gcf
if options.add_2_gcf==1
    % User called add_2_gcf
    figHandle = figure(gcf);
    hold on
elseif options.add_2_gcf==0
    % Default or user called add_2_gcf
    figHandle = figure;
    hold on
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

%% Retrieve data from Branch

for i = 1:numel(options.axes_indParam)
    % Each row is along an axis: x, y, z.
    % Each column is the value of that parameter or variable
    if isa(options.axes_indParam{i},'double')
        % A double corresponds to a parameter index
        plotVal(i,:) = arrayfun(@(p)p.parameter(options.axes_indParam{i}), ... 
            branch.point);
    elseif isa(options.axes_indParam{i},'char')
        % A string with 'x___' gives the index of a dynamic variable
        string = options.axes_indParam{i};
        indVar = str2double(string(2:end));
        plotVal(i,:) = arrayfun(@(p)p.x(indVar), ... 
            branch.point);
    end
end

% x = arrayfun(@(p)p.parameter(options.axes_indParam{1}), ... 
%     branch.point);
% y = arrayfun(@(p)p.parameter(options.axes_indParam{2}), ...
%     branch.point);
% z = arrayfun(@(p)p.parameter(options.axes_indParam{3}), ...
%     branch.point);


% YOU NEED TO UPDATE THIS THING SO IT DOESN'T ASSUME YOU ARE USING
% NUNST_COLOR ALL THE TIME.




%% Plot

if ~any(strcmp('nunst_color',p.UsingDefaults))
    % When nunst_color is used
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
        plot3(sel(plotVal(1,:),i),sel(plotVal(2,:),i),sel(plotVal(3,:),i), ...
            'Color',colors(i+1,:), options.PlotStyle{:});
    end
    hold off

    colormap(colors)
    colorbar()
    colorbar('YTickLabel',options.YTickLabel)
    
else
    % General plot
    plot3(plotVal(1,:),plotVal(2,:),plotVal(3,:), ...
        'Color',options.color, options.PlotStyle{:})
end


% Add title, axes
for i = 1:numel(options.axes_indParam)
    if isa(options.axes_indParam{i},'double')
        % simply a double gives a parameter value
        if i == 1
            xlabel([param.plot_names(options.axes_indParam{i}), ... 
                param.units(options.axes_indParam{i})])
        elseif i == 2
            ylabel([param.plot_names(options.axes_indParam{i}), ... 
                param.units(options.axes_indParam{i})])
        elseif i == 3
            zlabel([param.plot_names(options.axes_indParam{i}), ... 
                param.units(options.axes_indParam{i})])
        end
    elseif isa(options.axes_indParam{i},'char')
        % A string with 'x___' gives the index of a dynamic variable
        if i == 1
            xlabel(options.axes_indParam{i})
        elseif i == 2
            ylabel(options.axes_indParam{i})
        elseif i == 3
            zlabel(options.axes_indParam{i})
        end
    end
end



%% reset hold
if held_prior == 1
    hold on
elseif held_prior == 0
    hold off
end


end

