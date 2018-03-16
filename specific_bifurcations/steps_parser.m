function [ p ] = steps_parser( defaults, varargin )
%input parser for steps_() functions
% options:
% step, newton_max_iterations, max_bound, min_bound, halting_accuracy, 
% minimal_accuracy

p = inputParser;

% General option defaults
p.addParameter('step', defaults.step)
p.addParameter('max_step', defaults.max_step)
p.addParameter('newton_max_iterations', defaults.newton_max_iterations)
p.addParameter('max_bound', defaults.max_bound)
p.addParameter('min_bound', defaults.min_bound)
p.addParameter('halting_accuracy', defaults.halting_accuracy)
p.addParameter('minimal_accuracy', defaults.minimal_accuracy)

% Set options
parse(p,varargin{:})

end