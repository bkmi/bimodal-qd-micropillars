function [ branch ] = bif_contnAndStab( ...
    funcs, ...
    branch, ...
    numPoints, ...
    varargin )
%Continues the branch, either reversed or not. Outputs the branch with
%stability calculated.
%   Input:
%   funcs, ...
%   branch, ...
%   numPoints, ...
%   varargin
%
%   Options:
%       'reverse' = 0, 1
%           If reverse = 1 then the function will reverse the branch first.
%       'numWaves' = 2
%           For calculating stability, the smallest numWaves eigen values
%           are thrown out.
%% Input parser
p = inputParser;

% Options
p.addParameter('reverse', 0)
p.addParameter('numWaves', 2)

% Parse
parse(p,varargin{:})
options = p.Results;

%% Function

if options.reverse == 1
    branch = br_rvers(branch);
end

% Extend and calc stability
branch = br_contn(funcs, branch, numPoints);
[branch.nunst,~,~,branch.point] =  ...
    GetRotStability(branch, funcs, options.numWaves);

% Update indFold and indHopf
branch.indFold = find(abs(diff(branch.nunst))==1);
branch.indHopf = find(abs(diff(branch.nunst))==2);


end

