function [ directory ] = specific_bif_data_dir(  )
%Returns string of main directory.
if ispc
    filepath = mfilename('fullpath');
    parts = strsplit(filepath, '\');
    part = parts(1:end-1);
    directory = strjoin(part, '\');
    directory = strcat(directory, '\', 'specific_bif_plotting\');
else
    filepath = mfilename('fullpath');
    parts = strsplit(filepath, '/');
    part = parts(1:end-1);
    directory = strjoin(part, '/');
    directory = strcat(directory, '/', 'specific_bif_plotting/');
end
end

