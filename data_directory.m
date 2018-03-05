function [ directory ] = data_directory(  )
%Returns string of main directory.
filepath = mfilename('fullpath');
parts = strsplit(filepath, '\');
part = parts(1:end-1);
directory = strjoin(part, '\');
directory = strcat(directory, '\', 'output_data\');
end

