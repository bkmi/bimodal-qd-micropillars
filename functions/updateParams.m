function [ newParam ] = updateParams( param, varargin )
%Given a param structure and parameters to change, updateParams will
%produce a new parameter structure with the update made throughout the
%structure.
%
%   Inputs:
%       param, ...
%       varargin
%       
%   Options:
%       'par_name', par_value
%           This combination allows the user to call any param, including
%           'feed_ampliMatrix' and 'feed_phaseMatrix', and replace it with
%           a new value. Cannot be used in conjunction with
%           'par_overwrite'.
%
%       'par_overwrite', branch.point(1).parameter OR param.values
%           Calling this flag overrides the values given in the
%           param_struct. Cannot be used in conjunction with 'par_name',
%           par_value.
%           
%           REMEMBER: FEEDBACK MATRICES ARE NOT CONTAINED
%           IN par_overwrite!!
%
%           THE ORDER/INDICES DETERMINED IN param_struct MUST BE THE SAME 
%           AS IN YOUR 'par_overwrite' VALUES OR THERE WILL BE MASSIVE 
%           ERROR.
%
%   Outputs:
%       newParam
%% Collect input/defaults

% Create parser
p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;

% Add parameters
p.addParameter('par_overwrite', param.values)
p.addParameter('feed_ampliMatrix', [0, 0; 0, 0])
p.addParameter('feed_phaseMatrix', [0, 0; 0, 0])

% parse
parse(p,varargin{:})
newParam = param;

%% Organize based on par_overwrite or par_name, par_value

if ~any(strcmp('par_overwrite', p.UsingDefaults)) ...
        && numel(fieldnames(p.Unmatched)) ~= 0
    % When par_overwrite is chosen AND they called a par_name, par_value.
    error('You cannot use both par_overwrite and par_name, par_value)')
    
elseif ~any(strcmp('par_overwrite', p.UsingDefaults)) 
    % When par_overwrite is given.
    
    % Find the index of changes:
    ind_ParUpdt = find(param.values ~= p.Results.par_overwrite);
    
    for i = 1:numel(ind_ParUpdt)
        % values for loop
        ind_tempPar = ind_ParUpdt(i);
        tempParValue = p.Results.par_overwrite(ind_tempPar);
        tempParName = param.var_names{ind_tempPar};
        
        % fix in newParam
        newParam.values(ind_tempPar) = tempParValue;
        newParam.(tempParName).value = tempParValue;
    end
    
elseif numel(fieldnames(p.Unmatched)) ~= 0
    % When par_name, par_value is given.
    
    % Get all par_names
    parNames = fieldnames(p.Unmatched);
    
    for i = 1:numel(parNames)
        % values for loop
        tempParName = parNames{i};
        ind_tempPar = param.(tempParName).index;
        tempParValue = p.Unmatched.(tempParName);
        
        % fix in newParam
        newParam.values(ind_tempPar) = tempParValue;
        newParam.(tempParName).value = tempParValue;
    end
end

%% Specific for feed_phaseMatrix, feed_ampliMatrix
% Update the coupleing parameters cplPar

if ~any(strcmp('feed_ampliMatrix', p.UsingDefaults))
    newParam.cplPar.feed_ampliMatrix = p.Results.feed_ampliMatrix;
end

if ~any(strcmp('feed_phaseMatrix', p.UsingDefaults))
    newParam.cplPar.feed_phaseMatrix = p.Results.feed_phaseMatrix;
end
    
end

