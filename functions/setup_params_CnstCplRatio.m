function [ param ] = setup_params_CnstCplRatio( varargin )
%Setup the parameters for a bimodal electric field with coherent feedback.
%
%   By default this function will clear the base workspace, sets scripts in
%   DO NOT SAVE mode, and add all parameter values and indices to the
%   base workspace.
%
%   The default feedback parameters are as follows: 
%       feed_phase = 0
%       feed_ampli = 0.032 
%       tau_fb     = 0.8
%
%   Input:
%       varargin
%
%   Output:
%       param (struct)
%
%   Options:
%       'populate_wrkspc' = 0, 1
%           By default, this is set to 1. When 'populate_wrkspc' = 0, the
%           function will NOT add parameter values and indices to the
%           base workspace. When 'populate_wrkspc' = 1, the function will
%           add each parameter, param structure, and index to the
%           workspace.
%       'datadir_parent' = '../data_bimodal-qd-micropillars/'
%           By default, this is set as above. Determine the parent data 
%           saving directory.
%       'dimensional' = 1, 0
%           By default, this is set to 1. When 'dimensional' = 0, the
%           function applies a non-dimensionalized system. When
%           'dimensional' = 1, the function applies a dimensionalized
%           system.
%       'clear' = 0, 1
%           By default, this is set to 1. When 'clear' = 0, the function
%           leaves the workspace untouched. When 'clear' = 1, the base
%           workspace is cleared.
%       'save' = 0, 1
%           By default, this is set to 0. When 'save' = 0, the function
%           does not try to create a directory and does not try to save
%           anything. When 'save' = 1, the function tries to create a
%           directory and save master_options, parameters,
%           parameters_index, rotation_settings.
%

%% Options/Setup
p = inputParser;

% General option defaults
p.addParameter('datadir_parent','../data_bimodal-qd-micropillars/')
p.addParameter('dimensional',1)
p.addParameter('clear',1)
p.addParameter('save',0)
p.addParameter('populate_wrkspc', 1)

% Make options
parse(p,varargin{:})
options = p.Results;

% Add path, prepare matlab
addpath_setup()
if options.clear == 1
    evalin('base','clear all') % Clear workspace
end


%% Define default params/index


%define indicies
ind_kappa_s = 1;
ind_kappa_w = 2;
ind_mu_s = 3; 
ind_mu_w = 4;
ind_epsi_ss = 5; 
ind_epsi_ww = 6;
ind_epsi_sw = 7;
ind_epsi_ws = 8;
ind_beta = 9;
ind_J_p = 10;
ind_eta = 11;
ind_tau_r = 12;
ind_S_in = 13;
ind_V = 14;
ind_Z_QD = 15;
ind_n_bg = 16;
ind_tau_sp = 17;
ind_T_2 = 18;
ind_A = 19;
ind_hbar_omega = 20;
ind_epsi_tilda = 21;
ind_J = 22;
ind_feed_phase = 23;
ind_feed_ampli = 24;
ind_tau_fb = 25;
ind_feed_phaseMatrix = 26;
ind_feed_ampliMatrix = 27;
ind_epsi0 = 28;
ind_hbar = 29;
ind_e0 = 30;
ind_alpha_par = 31;

% Define params in kg, m, s aka SI UNITS.
% Physical constants
epsi0 = 8.85e-12;               % F m^-1 == J V^-2 m^-1
hbar = 1.05e-34;                % Js
e0 = 1.6e-19;                   % C
c0 = 3.0e8;                     % m s^-1
% Params from table 1 in Redlich
default_kappa_s    = 0.039*(1/1e-12);	% ps^-1 -> 1/s
default_kappa_w    = 0.041*(1/1e-12);	% ps^-1 -> 1/s
default_mu_s 	   = 3.70*(1e-9*e0);	% nm * e0 -> m*C
default_mu_w 	   = 3.75*(1e-9*e0);	% nm * e0 -> m*C
default_epsi_ss    = 70e-10;            % m^2 A^-1 V-^1
default_epsi_ww    = 50e-10;            % m^2 A^-1 V^-1
default_epsi_sw    = 160e-10;           % m^2 A^-1 V^-1
default_epsi_ws    = 150e-10;           % m^2 A^-1 V^-1
default_beta	   = 5.6e-3;
default_J_p 	   = 42.5*1e-6;         % microAmps -> Amps
default_eta 	   = 1.28e-3;
default_tau_r 	   = 150*(1e-12);       % ps -> s
default_S_in 	   = 10^(-16)*(1/1e-12);% m^2 ps^-1 -> m^2 s^-1
default_V          = 6.3*(1e-6)^3;      % micro m^3 -> m^3
default_Z_QD 	   = 110;
default_n_bg 	   = 3.34;
default_tau_sp 	   = 1*(1e-9);          % ns -> s
default_T_2 	   = 0.33*(1e-12);  	% ps -> s
default_A          = 3.14*(1e-6)^2; 	% micro m^2 -> m^2

% These two params cannot be altered by the user
hbar_omega = 1.38*(1.6e-19); 	% eV -> J
%epsi_tilda = epsi0*default_n_bg*c0; % Defined later!!
default_J          = 2.5*90*(1e-6);     % microAmps -> Amps
                                        % 2.5 * Threshold.
                                        % Threshold from Redlich paper,
                                        % *2.5 from ott10

% Feedback setup
default_feed_phase = 0;
default_feed_ampli = 0.373;
default_tau_fb 	   = 0.8;
% (Physical constants go here in index)
% matrix = [ss, sw; ws, ww]
default_feed_phaseMatrix = [1,0;0,1];
default_feed_ampliMatrix = [1,0;0,1];

% Further params
default_alpha_par  = 0;                 % Linewidth Enhancement Factor
default_omega      = 0;                 % Rotational parameter for DDEBIF


% Unit system
param_units = 'SI units used in definding parameters.';

if options.dimensional == 1 %dimensional units used in solver
    calc_units = ['Dimensional units used in calculations ',...
        'i.e. solver and bifurcations.'];
elseif options.dimensional == 0 %non-dimensional units used in solver
    calc_units = ['Non-dimensional units used in calculations ',...
        'i.e. solver and bifurcations.'];
else
    calc_units = ['No information was given regarding the use of ',...
        'dimensional or non-dimensional units in options.'];
    warning(['No information was given regarding the use of ',...
        'dimensional or non-dimensional units in options. ',...
        '\n This is recorded.']);
end





%% inputParser
% This section allows the user to set parameter values by using the syntax:
% 'parameter', value
p = inputParser;

% Add parameters
addParameter(p,'kappa_s', default_kappa_s);
addParameter(p,'kappa_w', default_kappa_w);
addParameter(p,'mu_s', default_mu_s);
addParameter(p,'mu_w', default_mu_w);
addParameter(p,'epsi_ss', default_epsi_ss);
addParameter(p,'epsi_ww', default_epsi_ww);
addParameter(p,'epsi_sw', default_epsi_sw);
addParameter(p,'epsi_ws', default_epsi_ws);
addParameter(p,'beta', default_beta);
addParameter(p,'J_p', default_J_p);
addParameter(p,'eta', default_eta);
addParameter(p,'tau_r', default_tau_r);
addParameter(p,'S_in', default_S_in);
addParameter(p,'V', default_V);
addParameter(p,'Z_QD', default_Z_QD);
addParameter(p,'n_bg', default_n_bg);
addParameter(p,'tau_sp', default_tau_sp);
addParameter(p,'T_2', default_T_2);
addParameter(p,'A', default_A);
%addParameter(p,'hbar_omega', hbar_omega);
%addParameter(p,'epsi_tilda', epsi_tilda);
addParameter(p,'J', default_J);
addParameter(p,'feed_phase', default_feed_phase);
addParameter(p,'feed_ampli', default_feed_ampli);
addParameter(p,'tau_fb', default_tau_fb);
addParameter(p,'feed_phaseMatrix', default_feed_phaseMatrix);
addParameter(p,'feed_ampliMatrix', default_feed_ampliMatrix);
%addParameter(p,'epsi0', epsi0);
%addParameter(p,'hbar', hbar);
%addParameter(p,'e0', e0);
%addParameter(p,'c0', c0);
addParameter(p,'alpha_par', default_alpha_par);
addParameter(p,'omega', default_omega);

% Add default options so parse doesn't freak out on us.
addParameter(p,'datadir_parent', 'option');
addParameter(p,'dimensional', 'option');
addParameter(p,'clear', 'option');
addParameter(p,'save', 'option');
addParameter(p,'populate_wrkspc','option');

% Parse inputs/options.
parse(p,options)

% Raise flags for values which are not allowed to be changed
if any(strcmp('epsi0',p.UsingDefaults))
    error('You cannot alter the value of epsi0')
elseif any(strcmp('hbar',p.UsingDefaults))
    error('You cannot alter the value of hbar')
elseif any(strcmp('e0',p.UsingDefaults))
    error('You cannot alter the value of e0')
elseif any(strcmp('c0',p.UsingDefaults))
    error('You cannot alter the value of c0')
elseif any(strcmp('hbar_omega',p.UsingDefaults))
    error('You cannot alter the value of hbar_omega')    
elseif any(strcmp('epsi_tilda',p.UsingDefaults))
    error('You cannot alter the value of epsi_tilda')
end


%% Update values, Create parameter struct and arrays
% Define params from user input + defaults
% Params from table 1 in Redlich
kappa_s    = p.Results.kappa_s;
kappa_w    = p.Results.kappa_w;
mu_s 	   = p.Results.mu_s;
mu_w 	   = p.Results.mu_w;
epsi_ss    = p.Results.epsi_ss;
epsi_ww    = p.Results.epsi_ww;
epsi_sw    = p.Results.epsi_sw;
epsi_ws    = p.Results.epsi_ws;
beta	   = p.Results.beta;
J_p 	   = p.Results.J_p;
eta 	   = p.Results.eta;
tau_r 	   = p.Results.tau_r;
S_in 	   = p.Results.S_in;
V          = p.Results.V;
Z_QD 	   = p.Results.Z_QD;
n_bg 	   = p.Results.n_bg;
tau_sp 	   = p.Results.tau_sp;
T_2 	   = p.Results.T_2;
A          = p.Results.A;
% hbar_omega = 1.38*(1.6e-19); %Defined above
epsi_tilda = epsi0*n_bg*c0;
J          = p.Results.J;       
% Feedback setup
feed_phase = p.Results.feed_phase;
feed_ampli = p.Results.feed_ampli;
tau_fb 	   = p.Results.tau_fb;
% (Physical constants go here in index)
% Further params
feed_phaseMatrix = p.Results.feed_phaseMatrix;
feed_ampliMatrix = p.Results.feed_ampliMatrix;
alpha_par  = p.Results.alpha_par;
omega      = p.Results.omega;

% Create parameter array called par
par_names = { ...
    'kappa_s', 'kappa_w', 'mu_s', 'mu_w', ...
    'epsi_ss', 'epsi_ww', 'epsi_sw', 'epsi_ws', ...
    'beta', 'J_p', 'eta', ...
    'tau_r', 'S_in', 'V', 'Z_QD', 'n_bg', ...
    'tau_sp', 'T_2', 'A', 'hbar_omega', 'epsi_tilda', 'J', ...
    'feed_phase', 'feed_ampli', 'tau_fb', ...
    'feed_phaseMatrix', 'feed_ampliMatrix', ...
    'epsi0', 'hbar', 'e0', ...
    'alpha_par' };
    
if options.dimensional == 1
    % Units for a dimensional run
    par_units = { ...
        '1/s', '1/s', 'mC', 'mC', ... 
        'm^2/(AV)', 'm^2/(AV)', 'm^2/(AV)', 'm^2/(AV)', ... 
        '(no units)', 'A', '(no units)', ... 
        's', 'm^2/s', 'm^3', '(no units)', '(no units)', ...
        's', 's', 'm^2', 'J', 'J/(V^2 s)', 'A', ...
        '(radians)', '(no units)', 's', ...
        '(no units)', '(no units)', ...
        'J/(V^2 m)', 'Js', 'C', ...
        '(no units)' };
elseif options.dimensional == 0
    % Units for a non-dimensional run
    par_units = { ...
        '1/s', '1/s', 'mC', 'mC', ... 
        'm^2/(AV)', 'm^2/(AV)', 'm^2/(AV)', 'm^2/(AV)', ... 
        '(no units)', 'A', '(no units)', ... 
        's', 'm^2/s', 'm^3', '(no units)', '(no units)', ...
        's', 's', 'm^2', 'J', 'J/(V^2 s)', 'A', ...
        '(radians)', '(no units)', 'tau_{sp}', ...
        '(no units)', '(no units)', ...
        'J/(V^2 m)', 'Js', 'C', ...
        '(no units)' };
end

par_plot_names = { ...
    '\kappa_s', '\kappa_w', '\mu_s', '\mu_w', ... 
    '\epsilon_{ss}', '\epsilon_{ww}', '\epsilon_{sw}', '\epsilon_{ws}', ...
    '\beta', 'J_p', '\eta', ...
    '\tau_r', 'S^{in}', 'V', 'Z^{QD}', 'n_{bg}', ...
    '\tau_{sp}', 'T_2', 'A', 'hbar\omega', '\epsilon_0n_{bg}c_{0}', 'J',...
    'Feedback Phase', 'Feedback Amp', '\tau_{fb}', ...
    'Phase Matrix', 'FB Amplitude Matrix', ...
    '\epislon0', 'hbar', 'e_0', ...
    '\alpha' };

par = { kappa_s, kappa_w, mu_s, mu_w, epsi_ss, epsi_ww, epsi_sw, epsi_ws, ...
    beta, J_p, eta, ...
    tau_r, S_in, V, Z_QD, n_bg, ...
    tau_sp, T_2, A, hbar_omega, epsi_tilda, J, ...
    feed_phase, feed_ampli, tau_fb, ...
    feed_phaseMatrix, feed_ampliMatrix, ...
    epsi0, hbar, e0, ...
    alpha_par };

% Append rotational parameters for DDEBIF tool
ind_omega=numel(par)+1;
par{ind_omega} = omega;
if options.dimensional == 1
    par_units(ind_omega) = {'1/s'}; % dimensional units
elseif options.dimensional == 0
    par_units(ind_omega) = {'1/tau_sp'}; % non-dimensional units
else
    error('Your dimensionality choice does not make sense!');
end
par_names(ind_omega) = {'omega'};
par_plot_names(ind_omega) = {'\omega'};

% Create a param struct with all of this information.
param.values = par;
param.var_names = par_names;
param.units = par_units;
param.plot_names = par_plot_names;
param.unit_system = [param_units, ' ', calc_units];

for i=1:length(param.var_names)
    param.(param.var_names{i}) = struct;
    param.(param.var_names{i}).index = i;
    param.(param.var_names{i}).value = param.values(i);
    param.(param.var_names{i}).var_name = param.var_names{i};    
    param.(param.var_names{i}).units = param.units{i};
    param.(param.var_names{i}).plot_name = param.plot_names{i};
end


%% Rotational functionality for DDEBIF 

%% I am really worried about this part
% I don't think that this system is correct. I believe I will need to use
% multiple omegas and a more complicated rotational symmetry. Something
% like this:
% 
%     x = [re(E_strong);
%            im(E_strong):
%            re(E_weak);
%            im(E_weak):
%            rho;
%            n ]
% 
%     A =
%         [0,-1,0, 0,0,0;
%          1, 0,0, 0,0,0;
%          0, 0,0,-1,0,0;
%          0, 0,1, 0,0,0;
%          0, 0,0, 0,0,0;
%          0, 0,0, 0,0,0 ];
% 
%     expA=@(phi, theta)
%         [cos(phi),-sin(phi),0         ,0          ,0,0; ...
%          sin(phi), cos(phi),0         ,0          ,0,0; ...
%          0       ,0        ,cos(theta),-sin(theta),0,0; ...
%          0       ,0        ,sin(theta),cos(theta) ,0,0; ...
%          0       ,0        ,0         ,0          ,0,0; ...
%          0       ,0        ,0         ,0          ,0,0 ];

% create rotation matrix
A_rot = ...
    [0,-1,0, 0,0,0; ...
     1 ,0,0, 0,0,0; ...
     0 ,0,0,-1,0,0; ...
     0 ,0,1, 0,0,0; ...
     0 ,0,0, 0,0,0; ...
     0 ,0,0, 0,0,0 ];
expA_rot=@(phi) ...
    [cos(phi),-sin(phi),0       ,0        ,0,0; ...
     sin(phi), cos(phi),0       ,0        ,0,0; ...
     0       ,0        ,cos(phi),-sin(phi),0,0; ...
     0       ,0        ,sin(phi), cos(phi),0,0; ...
     0       ,0        ,0       ,0        ,0,0; ...
     0       ,0        ,0       ,0        ,0,0 ];
    

% Choose system based on dimensional choice
if options.dimensional == 1
    %define rhs ready ddebif, dimensional units
    rhs = @(x,p)bimodalSystem_CnstCplRatio(...
    x(1,1,:)+1i*x(2,1,:), ... Es
    x(1,2,:)+1i*x(2,2,:), ... EsTau 
    x(3,1,:)+1i*x(4,1,:), ... Ew
    x(3,2,:)+1i*x(4,2,:), ... EwTau
    x(5,1,:), ... rho
    x(6,2,:), ... n
    p{:});

elseif options.dimensional == 0
    
    error('Not written yet')
    
    %define rhs ready func, dimensionless
%     rhs = @(x,p)qd_1ef_phaseAmp_nondim(x(1,1,:)+1i*x(2,1,:),...
%         x(1,2,:)+1i*x(2,2,:),x(3,1,:),x(3,2,:),x(4,1,:),x(4,2,:),...
%         p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),...
%         p(12),p(13),p(14),p(15),p(16),p(17),p(18),p(19),p(20),p(21),...
%         p(22),p(23),p(24),p(25),p(26),p(27),p(28),p(29),p(30));
else
    error('Your dimensionality choice does not make sense!')
end

% Prepare 'funcs' for DDEBIF
funcs = set_rotfuncs('sys_rhs',rhs, ... 
    'rotation',A_rot,'exp_rotation',expA_rot, ... 
    'sys_tau',@()ind_tau_fb,'x_vectorized',true);
opt_inputs = {'extra_condition',1,'print_residual_info',0};

%% Saving Section

% Determine name of folder with relevant parameters.
% mono mode (single) electric field
mode_report = 'BiModal_';

% report dimensionality
if options.dimensional == 1
    dimension_report = 'dim_'; % dimensional units
elseif options.dimensional == 0
    dimension_report = ''; % non-dimensional units
else
    error('Your dimensionality choice does not make sense!');
end

% report current IN AMPS
current_report = strcat('J=',num2str(J,'%1.1e'),'_');

% report feedback params
feed_tau_report = strcat('tau=',num2str(tau_fb),'_');
feed_amp_report = strcat('amp=',num2str(feed_ampli));

% report alpha (line width enchancement) if NONZERO
if alpha_par == 0
    alpha_par_report = '';
elseif alpha_par ~= 0
    alpha_par_report = strcat('_alpha=',num2str(alpha_par));
end

% Folder shall be named below:
datadir_specific = strcat(options.datadir_parent, ... 
    mode_report,dimension_report,current_report, ...
    'FEED_',feed_tau_report,feed_amp_report, ...
    alpha_par_report, '/');


% Make user confirm overwrite
if options.save == 1
    if exist(datadir_specific,'dir') == 7
        while(1)
            overwrite = input( ... 
                [ '\n\nWARNING: Directory already exists, ' ...
                'would you like to overwrite data? ', ... 
                '\n0 = no \n1 = yes \n\n'] );
            if overwrite == 1 % overwrite
                fprintf(...
                    '\n\nYou have chosen to OVERWRITE the files!\n\n')
                break
            elseif overwrite == 0 % DO NOT overwrite
                fprintf([ '\n\nYou have chosen not to ', ...
                    'overwrite the files. ', ...
                    '\nThe program WILL NOT SAVE.\n\n\n' ])
                options.save = 0; %turn off saving
                break
            end
        end
    else
        mkdir(datadir_specific)
    end
end


% Create 'master_options' structure
master_options = struct;
master_options.save = options.save;
master_options.datadir_parent = options.datadir_parent;
master_options.datadir_specific = datadir_specific;
master_options.dimensional = options.dimensional;


% Save
if options.save == 1
    % Where will it save?
    fprintf(strcat('Subfolder:\n', datadir_specific,'\n'))

    % Save parameter index
    save(strcat(datadir_specific,'parameters_index.mat'),...
        'ind_kappa_s','ind_kappa_w','ind_mu_s', ...
        'ind_mu_w','ind_epsi_ss','ind_epsi_ww','ind_epsi_sw',...
        'ind_epsi_ws','ind_beta','ind_J_p','ind_eta','ind_tau_r',...
        'ind_S_in','ind_V','ind_Z_QD','ind_n_bg','ind_tau_sp',...
        'ind_T_2','ind_A','ind_hbar_omega','ind_epsi_tilda','ind_J',...
        'ind_feed_phase','ind_feed_ampli','ind_tau_fb','ind_epsi0',...
        'ind_hbar','ind_e0','ind_c0','ind_alpha_par','ind_omega')
    
    % Save parameters
    save(strcat(datadir_specific,'parameters.mat'),...
        'kappa_s','kappa_w','mu_s', ...
        'mu_w','epsi_ss','epsi_ww','epsi_sw',...
        'epsi_ws','beta','J_p','eta','tau_r',...
        'S_in','V','Z_QD','n_bg','tau_sp',...
        'T_2','A','hbar_omega','epsi_tilda','J',...
        'feed_phase','feed_ampli','tau_fb','epsi0',...
        'hbar','e0','c0','alpha_par','omega', ...
        'param')
    
    % Save rotational settings
    save(strcat(datadir_specific,'rotation_settings.mat'), ... 
        'rhs','funcs','opt_inputs','A_rot','expA_rot')
    
    % Save options
    save(strcat(datadir_specific,'master_options.mat'),...
        'master_options')
    
end


%% Populate Base Workspace

if options.populate_wrkspc == 1
    
    assignin('base','param', param);
    
    % parameters
    assignin('base','kappa_s', kappa_s);
    assignin('base','kappa_w', kappa_w);
    assignin('base','mu_s', mu_s);
    assignin('base','mu_w', mu_w);
    assignin('base','epsi_ss', epsi_ss);
    assignin('base','epsi_ww', epsi_ww);
    assignin('base','epsi_sw', epsi_sw);
    assignin('base','epsi_ws', epsi_ws);
    assignin('base','beta', beta);
    assignin('base','J_p', J_p);
    assignin('base','eta', eta);
    assignin('base','tau_r', tau_r);
    assignin('base','S_in', S_in);
    assignin('base','V', V);
    assignin('base','Z_QD', Z_QD);
    assignin('base','n_bg', n_bg);
    assignin('base','tau_sp', tau_sp);
    assignin('base','T_2', T_2);
    assignin('base','A', A);
    assignin('base','hbar_omega', hbar_omega);
    assignin('base','epsi_tilda', epsi_tilda);
    assignin('base','J', J);
    assignin('base','feed_phase', feed_phase);
    assignin('base','feed_ampli', feed_ampli);
    assignin('base','tau_fb', tau_fb);
    assignin('base','feed_phaseMatrix', feed_phaseMatrix);
    assignin('base','feed_ampliMatrix', feed_ampliMatrix);
    assignin('base','epsi0', epsi0);
    assignin('base','hbar', hbar);
    assignin('base','e0', e0);
    assignin('base','alpha_par', alpha_par);
    assignin('base','omega', omega);
    
    % indicies
    assignin('base','ind_kappa_s', ind_kappa_s);
    assignin('base','ind_kappa_w', ind_kappa_w);
    assignin('base','ind_mu_s', ind_mu_s);
    assignin('base','ind_mu_w', ind_mu_w);
    assignin('base','ind_epsi_ss', ind_epsi_ss);
    assignin('base','ind_epsi_ww', ind_epsi_ww);
    assignin('base','ind_epsi_sw', ind_epsi_sw);
    assignin('base','ind_epsi_ws', ind_epsi_ws);
    assignin('base','ind_beta', ind_beta);
    assignin('base','ind_J_p', ind_J_p);
    assignin('base','ind_eta', ind_eta);
    assignin('base','ind_tau_r', ind_tau_r);
    assignin('base','ind_S_in', ind_S_in);
    assignin('base','ind_V', ind_V);
    assignin('base','ind_Z_QD', ind_Z_QD);
    assignin('base','ind_n_bg', ind_n_bg);
    assignin('base','ind_tau_sp', ind_tau_sp);
    assignin('base','ind_T_2', ind_T_2);
    assignin('base','ind_A', ind_A);
    assignin('base','ind_hbar_omega', ind_hbar_omega);
    assignin('base','ind_epsi_tilda', ind_epsi_tilda);
    assignin('base','ind_J', ind_J);
    assignin('base','ind_feed_phase', ind_feed_phase);
    assignin('base','ind_feed_ampli', ind_feed_ampli);
    assignin('base','ind_tau_fb', ind_tau_fb);
    assignin('base','ind_feed_phaseMatrix', ind_feed_phaseMatrix);
    assignin('base','ind_feed_ampliMatrix', ind_feed_ampliMatrix);
    assignin('base','ind_epsi0', ind_epsi0);
    assignin('base','ind_hbar', ind_hbar);
    assignin('base','ind_e0', ind_e0);
    assignin('base','ind_alpha_par', ind_alpha_par);
    assignin('base','ind_omega', ind_omega);
    
    % rotational settings
    assignin('base','rhs', rhs);
    assignin('base','funcs', funcs);
    assignin('base','opt_inputs', opt_inputs);
    assignin('base','A_rot', A_rot);
    assignin('base','expA_rot', expA_rot);
    
    % Options
    assignin('base', 'master_options', master_options)

end


end

