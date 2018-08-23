function [  ] = plot_RPO3( branch, param, varargin )
%A hack which plots the minimum and maximum of an RPO as two branches
%   It works specifically when for strong mode on the z axis
x1_norm = @(x)sqrt(sum(x(1:2,:).^2,1));
bmax = arrayfun(@(p)max(x1_norm(p.profile)), branch.point);
for i = 1:numel(bmax)
    branch.point(i).x(1) = bmax(i);
end
plot_branch3(branch, param, varargin{:})

bmin = arrayfun(@(p)min(x1_norm(p.profile)), branch.point);
for i = 1:numel(bmin)
    branch.point(i).x(1) = bmin(i);
end
plot_branch3( ...
    branch, ...
    param, ...
    'add_2_gcf', 1, ...
    'PlotStyle', {'LineStyle', '--', 'Marker', 'o'}, ...
    varargin{:})

end

