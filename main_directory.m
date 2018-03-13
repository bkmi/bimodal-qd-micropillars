function [ directory ] = main_directory(  )
%Returns string of main directory.
if ispc
    filepath = mfilename('fullpath');
    parts = strsplit(filepath, '\');
    part = parts(1:end-1);
    directory = strjoin(part, '\');
    directory = strcat(directory, '\');
else
    filepath = mfilename('fullpath');
    parts = strsplit(filepath, '/');
    part = parts(1:end-1);
    directory = strjoin(part, '/');
    directory = strcat(directory, '/');
end
end

