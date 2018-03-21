function [ inds ] = branch_pts_near_param( branch, par_index, par_value, tol )
%I.e. search for points in a branch near fbA = 0.3
inds = [];
for i = 1:numel(branch.point)
    v = branch.point(i).parameter(par_index);
    if abs(par_value - v) < tol
        inds(end+1) = i;
    end
end
end

