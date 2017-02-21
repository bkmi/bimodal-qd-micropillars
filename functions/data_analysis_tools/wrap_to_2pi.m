function [ out_branch, out_struct ] = wrap_to_2pi( branch, ...
    ind_wrap,...
    wrapped_branch_name, container_struct)
%Wraps the parameter in the index to 2pi. The first output is the wrapped
%branch itself. The second output is a struct containing the wrapped
%branch.
%
%Calling wrapped_branch_name defines the name of the branch in the
%out_struct. Calling container_struct will set out_struct to
%container_struct plus the out_branch. (out_branch according to
%wrapped_branch_name.)
%
%   Input:
%       branch, ...
%       ind_wrap, ...
%       ('wrapped_branch_name'), ...
%       (container_struct)
%
%   Output:
%       out_branch, ...
%       container_struct
%


switch nargin
    case 4
        if ~isa(wrapped_branch_name,'char')
            % Only accept a 'name'
            error('wrapped_branch_name must be a char.')
        end
        
        if ~isa(container_struct,'struct')
            % Only accept a struct
            error('container_struct must be a struct.')
        end
        out_struct = container_struct;
        
        if any(strcmp(wrapped_branch_name,fieldnames(out_struct)))
            % Check if the function will overwrite a field
            error('That fieldname is already taken.')
        end
        
    case 3
        if ~isa(wrapped_branch_name,'char')
            % Only accept a 'name'
            error('wrapped_branch_name must be a char.')
        end
        
        % Create default out_struct
        out_struct = struct;
        
    otherwise
        % Create default wrapped_branch_name and out_struct.
        wrapped_branch_name = 'wrap_branch';
        out_struct = struct;
end


% Wrap points in the branch
for i = 1:length(branch.point)
    branch.point(i).parameter(ind_wrap) = ...
        mod(branch.point(i).parameter(ind_wrap),2*pi);
end

out_branch = branch;


% Add out_branch to struct.
out_struct.(wrapped_branch_name) = out_branch;


end

