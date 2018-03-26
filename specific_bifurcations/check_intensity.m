function [ inten ] = check_intensity( ecm_cell )
%Returns an array with intensity values at the second point in that cell
shape = size(ecm_cell);
inten = zeros(shape);

for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(ecm_cell{i,j})
            inten(i,j) = ecm_cell{i,j}.point(2).x(1);
        end
    end
end

end

