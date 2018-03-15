function [ branch ] = branch_stability( funcs, branch )
%Gets stability and updates the branch
[branch.nunst,~,~,branch.point] = GetRotStability( ...
    branch, ...
    funcs, 1);
branch.indFold = find(abs(diff(branch.nunst))==1);
branch.indHopf = find(abs(diff(branch.nunst))==2);
% branch.error = 0;

end

