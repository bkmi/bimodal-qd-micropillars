function [nunst,dom,triv_defect,points] = GetRotStability( branch, ...
    funcs, numOfWaves )
%Expand on DDEBIF GetStability. The options to account for rotating waves
%are all checked when using this function.
%
%numOfWaves = 1 in single mode case
%numOfWaves = 2 in bi modal case
%
%This is with the follow options marked:
%   GetStability(branch,...
%       'exclude_trivial',true,'locate_trivial',@(p)0,'funcs',funcs);
%
%   Input:
%       branch, ...
%       funcs, ...
%       numOfWaves
%
%   Output:
%       nunst
%       dom
%       triv_defect
%       points

numTrivial = zeros(numOfWaves,1);

[nunst,dom,triv_defect,points] = GetStability(branch, ...
    'exclude_trivial',true,'locate_trivial',@(p)numTrivial,'funcs',funcs);


end

