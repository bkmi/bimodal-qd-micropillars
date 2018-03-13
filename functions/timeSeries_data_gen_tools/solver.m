function [ dde23_soln, figHandle ] = solver( hist, timeSpan, ...
    param, ...
    varargin )
%Takes in parameters from the experimental setup and returns a time
%series plot starting at the values given by a "hist" vector. If "hist" is
%zero, except for electric field, then this is a "turn on" timeseries.
%
%This function requires a param structure input for organizational
%purposes. Note: If you would like to create a timeseries with a different
%parameter set than the one saved in the param structure (i.e. from a
%bifurcation continuation in DDEBIF) then use the 'par_overwrite' option
%below!!
%
%THE ORDER/INDICES DETERMINED IN param_struct MUST BE THE SAME AS IN YOUR
%'par_overwrite' VALUES OR THERE WILL BE MASSIVE ERROR.
%
%   Input:
%       hist
%       timeSpan
%       param
%       varargin
%
%   Options:
%       'par_overwrite' = branch.point(1).parameter OR param.values
%           Calling this flag overrides the values given in the
%           param_struct. This is particularly useful if you would like to
%           make a time series starting at a point along a DDEBIF
%           bifurcation continuation.
%
%           THE ORDER/INDICES DETERMINED IN param_struct MUST BE THE SAME 
%           AS IN YOUR 'par_overwrite' VALUES OR THERE WILL BE MASSIVE 
%           ERROR.
%
%       'plot' = 1,0
%           If you choose to plot = 1 then the solver will output the
%           solver plot. These are timeseries with ef, rho, and n.
%       'dde23_options' = ddeset('RelTol',10^-8), ...;
%           Calling this flag will let you change the precision, behavior
%           of the dde23 solver. The default is written above. For more
%           information check out the following link:
%         https://www.mathworks.com/help/matlab/ref/ddeset.html#f81-1031913
%           ddeset('RelTol',10^-8,'OutputFcn', @odeplot)
%       'quiet' = 0, 1
%           0 -> usual output.
%           1 -> no disp output.





% Add path
addpath_setup()


%% Defaults + inputParser + Organize behavior

p = inputParser;

% General option defaults
p.addParameter('par_overwrite',0)
p.addParameter('plot',0)
p.addParameter('dde23_options',ddeset('RelTol',10^-8))
p.addParameter('save_name', 'dde23_soln')
p.addParameter('quiet', 0)

% Master option defaults
% historical reasons
p.addParameter('save', 0)
p.addParameter('datadir_parent', strcat(data_directory(), 'output_data/'))
p.addParameter('datadir_specific', strcat(data_directory(), 'output_data/'))
p.addParameter('dimensional',0)

parse(p,varargin{:})
options = p.Results;

% Set par equal to the relevant parameter set.
if options.par_overwrite == 0
    par = param.values;
elseif length(options.par_overwrite) ~= length(param.values)
    error([ 'The length of param_struct.values does not agree with ', ...
        'the length of your par_overwrite input.'])
else
    par = options.par_overwrite;
end

% Set the coupling matricies
feed_phaseMatrix = param.cplPar.feed_phaseMatrix; 
feed_ampliMatrix = param.cplPar.feed_ampliMatrix;

%define solver ready func, dimensionless
sys_4solver = @(x)nonDim_bimodalSystem_CnstCplRatio(...
    x(1,1,:)+1i*x(2,1,:), x(1,2,:)+1i*x(2,2,:),... %EF1
    x(3,1,:)+1i*x(4,1,:), x(3,2,:)+1i*x(4,2,:),... %EF2
    x(5,1,:),x(6,1,:),...
    feed_phaseMatrix, feed_ampliMatrix, ...
    par(1),par(2),par(3),par(4),par(5),par(6),par(7),par(8),par(9), ...
    par(10),par(11),par(12),par(13),par(14),par(15),par(16), ...
    par(17),par(18),par(19),par(20),par(21),par(22),par(23), ...
    par(24),par(25),par(26),par(27),par(28),par(29));

if options.quiet == 0
    fprintf('Solver is using dimensionless units for output.\n')
end

dde23_ef_units = '(\epsilon_{ss} \epsilon_{tilda})^{-1/2}';
dde23_time_units = '(\tau_{sp})';
dde23_n_units = '(S^{in} \tau_{sp})^{-1}';

%% Solver + Plotter
% Setup/use dde23 solver
if options.quiet == 0
    if isa(hist,'function_handle')
        % If hist is a function
        fprintf('History vector is a function. \n')

    elseif isa(hist,'double')
        % If hist is a vector
        fprintf(strcat('History vector: \nhist=',mat2str(hist),'\n'))
    else
        disp('I do not know what type hist is.')
    end
end
    
lags = par(param.tau_fb.index); % lags == feedback time
dde23_soln = dde23(@(t,y,z)sys_4solver([y,z]),...
    lags,hist,timeSpan, options.dde23_options);


% Plot, if necessary
if options.plot == 1
    figHandle = plot_solver( dde23_soln, par(param.J.index), ...
        hist, ...
        dde23_time_units, dde23_ef_units, dde23_n_units );
else
    figHandle = 0;
end

end