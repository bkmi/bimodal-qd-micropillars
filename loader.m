function [ fileData ] = loader( varargin )
%Load data from the data directory.
%   The loader will load all save files from a folder.
%   If you call it with a single output that output will be a struct with
%   all of the saved data.
%
%   Options:
%       'datadir_parent' = str (location of parent datadir)
%           This is unnecessary if you called datadir_specific. By default
%           the directory '../data_qd-micropillar-laser-ddebif/' is chosen.
%       'datadir_specific' = str (location of the datadir you want to load)
%           This dir will be tried before trying datadir_parent.
%       'exclude' = {str, str, ...} (cell array of files not to include)
%           These files will not be loaded.
%       'include' = {str, str, ...} (cell array of files include)
%           These files will be loaded. Nothing else will be loaded from
%           the folder.
%       'overwrite' = 0, 1
%           When overwrite is 0, the function will ask if the user wants to
%           overwrite files. When overwrite is 1 the function will
%           automatically overwrite.
%       

%% Options
p = inputParser;

% General option defaults
p.addParameter('include', cell(0));
p.addParameter('exclude', cell(0));
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','')
p.addParameter('overwrite',0)

% Set options
parse(p,varargin{:})

% Make final options
options = p.Results;


% Add trailing / to datadir_specific and  datadir_parent
if ~any(strcmp('datadir_specific',p.UsingDefaults)) ...
    && options.datadir_specific(end) ~= '/'
    options.datadir_specific(end+1) = '/';
end

if ~any(strcmp('datadir_parent',p.UsingDefaults)) ...
        && options.datadir_parent(end) ~= '/'
    options.datadir_parent(end+1) = '/';
end


% Organize behavior from options
% Check for datadir_parent
if ~any(strcmp('datadir_parent',p.UsingDefaults)) ...
        && ~isdir(options.datadir_parent)
    % When datadir_parent is non-default and doesn't exist.
    warning(strcat(options.datadir_parent, ...
        ': Is not a dir. Using default parent directory.'))
elseif any(strcmp('datadir_parent',p.UsingDefaults)) ...
        && ~isdir(options.datadir_parent)
    % When datadir_parent is default and doesn't exist.
    error('default parent dir does not exist. update your code.')
    % Ideally you would let the user choose the parent directory here, but
    % I am too lasy to write that code.
end
% Check for datadir_specific
if ~any(strcmp('datadir_specific',p.UsingDefaults)) ...
        && ~isdir(options.datadir_specific)
    % When the datadir_specific is non-default and doesn't exist.
    warning(strcat(options.datadir_specific, ...
        ': Is not a dir. Using regular behavior.'))
end


% Select data folder location
if ~isdir(options.datadir_specific) ...
        || any(strcmp('datadir_specific',p.UsingDefaults))
    % If options.datadir_specific doesn't point to a directory
    directory = dir(options.datadir_parent);
    sub_all_list = {directory.name};
    which_are_folders = [directory.isdir];
    subfolder_list = sub_all_list(which_are_folders);
    [selection,exit_box] = listdlg('PromptString','Select a folder:',...
                                            'SelectionMode','single',...
                                            'ListString',subfolder_list);
    datadir_subfolder = strcat(subfolder_list{selection},'/');
    if exit_box==0
        error('You did not choose a folder. Nothing was loaded')
    end
    loadDir = strcat(options.datadir_parent,datadir_subfolder);
else
    loadDir = options.datadir_specific;
end


%% Prepare matlab + Load
% Prepare matlab
addpath_setup()


% Load
datadir = dir(loadDir);
fileIndex = find(~[datadir.isdir]);
% Load only 'included' things
if ~any(strcmp('include',p.UsingDefaults))
    % When include is non-default. Load included things.
    names = {datadir.name};
    fileNames = names(fileIndex);
    for i = 1:length(options.include)
        include_fileName = options.include{i};
        if any(strcmp(include_fileName,fileNames))
            try
                data.(include_fileName(1:end-4)) = ...
                    load(strcat(loadDir,include_fileName));
                disp(include_fileName);
            catch ME
                switch ME.identifier
                    case 'MATLAB:load:numColumnsNotSame'
                        warning(strcat(include_fileName, ... 
                            [': Failed to load. ', ...
                            'MATLAB:load:numColumnsNotSame']))
                    otherwise
                        rethrow(ME)
                end
            end
        elseif ~any(strcmp(include_fileName,fileNames))
            warning(strcat(include_fileName, ...
                ': No matches to that filename.'))
        end
    end
else
    % Otherwise, Load everything.
    for i = 1:length(fileIndex)
        fileName = datadir(fileIndex(i)).name;
        if ~any(strcmp(fileName, options.exclude))
            try
                data.(fileName(1:end-4)) = ... 
                    load(strcat(loadDir,fileName));
            catch ME
                switch ME.identifier
                    case 'MATLAB:load:numColumnsNotSame'
                        warning(strcat(fileName,...
                            [': Failed to load. ',...
                            'MATLAB:load:numColumnsNotSame']))
                    otherwise
                        rethrow(ME)
                end
            end
        else
            warning(strcat(fileName,': Was excluded.'))
        end
    end
end


% Display data list aka files names.
fprintf('\n\n\n\nFiles loaded:\n')
disp(data)
fprintf('\n\n')

if any(strcmp(p.UsingDefaults, 'overwrite'))
    % Would the user like to overwrite files?
    loaded_overwrite_opt = ...
        input(['\n\nWARNING: You have loaded this data.',...
        '\nShould scripts overwrite saved results? \n0 = no \n1 = yes \n\n']);
    % Force user to choose: OVERWRITE or not.
    while(1)
        if loaded_overwrite_opt == 1
            saveit = 1; %Save and OVERWRITE
            fprintf('\nScripts will OVERWRITE saved results!!!\n\n')
            break
        elseif loaded_overwrite_opt == 0
            saveit = 0; %Don't save and don't overwrite.
            fprintf('\nScripts will NOT save results!!!\n\n')
            break
        end
    end
elseif options.overwrite == 0 ...
        && ~any(strcmp(p.UsingDefaults, 'overwrite'))
    saveit = 0; %Don't save and don't overwrite.
    fprintf('\nScripts will NOT save results!!!\n\n')
elseif options.overwrite == 1
    % Overwrite, but don't ask!!
    saveit = 1;
    fprintf('\nScripts will OVERWRITE saved results!!!\n\n')
end


% Set master_options based on above selection
if any(strcmp('master_options',fieldnames(data)))
    % If there is a master_options, load s.t. saving is turned off.
    data.master_options.master_options.save = saveit;
else
    % If there is no master_options, make one.
    data.master_options = struct;
    data.master_options.master_options = struct;
    data.master_options.master_options.save = saveit;
    data.master_options.master_options.datadir_specific = loadDir;
end


% Did they call with an output?
switch nargout
    case 0 % No output arguments, load into workspace.
    % Bring it all into the main workspace
    savefileNames = fieldnames(data);
    for i=1:length(fieldnames(data))
        varNames = fieldnames(data.(savefileNames{i}));
        for j=1:length(varNames)
            assignin('base', varNames{j}, ...
                data.(savefileNames{i}).(varNames{j}))
        end
    end
    case 1 % They called with an output arg, don't load into workspace
        fileData = data;
    otherwise
        error(['Use no output to load into workspace, ',...
            'one output to load into argument'])
end


end

