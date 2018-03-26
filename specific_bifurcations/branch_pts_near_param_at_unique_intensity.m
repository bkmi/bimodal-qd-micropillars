function [ inds ] = branch_pts_near_param_at_unique_intensity(branch, par_index, par_value, partol, intentol)
%As the title, avoids returning degenerate cases.

% ind = branch_pts_near_param(branch, par_index, par_value, partol);
%     
% intens = [];
% for j = 1:numel(ind)
%     intens(end + 1) = branch.point(ind(j)).x(1);
% end
% 
% num_decimals_to_round = int_num_dec_round;
% 
% rounder = 10.^num_decimals_to_round;
% intens = round(rounder * intens)/rounder;
% [~,unqiue_inds,~] = unique(intens);
% 
% inds = ind(unqiue_inds);


%%
ind = branch_pts_near_param(branch, par_index, par_value, partol);
    
intens = [];
inds = [];
for j = 1:numel(ind)
    inten = branch.point(ind(j)).x(1);
    if isempty(intens) || (min(abs(round_to(inten, 8) - round_to(intens, 8))) > intentol)
        intens(end + 1) = branch.point(ind(j)).x(1);
        inds(end + 1) = ind(j);
    end
end

end

function [ rnded ] = round_to(val, decimal_place)
    rounder = 10.^decimal_place;
    rnded = round(rounder * val)/rounder;
end

