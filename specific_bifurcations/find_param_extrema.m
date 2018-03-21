function [ ptmax, ptmin ] = find_param_extrema( branch, par_index )
%Return the point numbers for extrema of a parameter in the branch

max = branch.point(1).parameter(par_index);
min = branch.point(1).parameter(par_index);
ptmax = 1;
ptmin = 1;

for i = 2:numel(branch.point)
    if branch.point(i).parameter(par_index) > max
        max = branch.point(i).parameter(par_index);
        ptmax = i;
    end
    
    if branch.point(i).parameter(par_index) < min
        min = branch.point(i).parameter(par_index);
        ptmin = i;
    end
end

end

