function [ tree ] = clean_bif_tree( tree )
%Removes branches with errors from a set of bifurcation trees.

for i = 1:numel(tree)
    for j = 1:numel(tree{i})
        if tree{i}(j).error == 1
            tree{i}(j) = [];
        end
    end
end

end

