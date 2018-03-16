function [ branch ] = bcont( funcs, branch, numForward, numRvers )
%Continue a branch forward and backward
[branch,~,~,~] = br_contn(funcs, branch, numForward);
if numRvers > 0
    branch = br_rvers(branch);
    [branch,~,~,~] = br_contn(funcs,branch,100);
end

end

