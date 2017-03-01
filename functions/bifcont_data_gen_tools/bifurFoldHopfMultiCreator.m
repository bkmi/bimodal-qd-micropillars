function [ multiBranch ] = bifurFoldHopfMultiCreator( ...
    funcs, ststBranch, param, indBifurcations, lengthBranch, ...
    varargin)
%Creates a package of bifurcation continuations from each bifurcation point
%   Input:
%
%   Options:
%       'indCont' = [param.feed_phase.index, param.feed_ampli.index]
%   
%% Options
p = inputParser;

% General option defaults
p.addParameter('indCont', [param.feed_phase.index, param.feed_ampli.index])

% Master option defaults
p.addParameter('save',0)
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('datadir_specific','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',0)

% Parse, set options
parse(p,varargin{:})
options = p.Results;


%% Calculate

multiBranch = struct(...
    'method',struct, ...
    'parameter',struct, ...
    'point',struct, ...
    'newFuncs',@()'undefined',...
    'error', 0, ...
    'indError', NaN);
multiBranch = repmat(multiBranch, [numel(indBifurcations), 1]);

for i = 1:length(indBifurcations)
    % Update params to correspond with the values at the bifurcation.
    param = updateParams(param, ...
        'par_overwrite', ststBranch.point(indBifurcations(i)).parameter);
    
    try
        [branch,newFuncs] = ...
            bifurContin_FoldHopf( ...
            funcs, ... 
            ststBranch, ...
            indBifurcations(i), ...
            options.indCont, ...
            lengthBranch, ...
            param,...
            'plot_prog', 1, ...
            'save',0);

        multiBranch(i).method = branch.method;
        multiBranch(i).parameter = branch.parameter;
        multiBranch(i).point = branch.point;
        multiBranch(i).newFuncs = newFuncs;

        
    catch ME
        switch ME.identifier
            case 'br_contn:start'
                warning(ME.message);
                warning( ...
                    strcat('During bifurFoldHopfMultiCreator. i=', ...
                    num2str(i)));
                multiBranch(i).error = ME;
                multiBranch(i).indError = i;
            otherwise
                rethrow(ME)
        end
    end
end

end

