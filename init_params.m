function [ params ] = init_params( varargin )
%Choose parameters for a non-dimensionalized, two mode polarization,
%electric field with coherent feedback.

end

function [ p ] = parser( varargin )
%% Define default params/index
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

% Feedback setup
default_feed_phase = 0;
default_feed_ampli = 0.373;
default_tau_fb 	   = 0.8;
default_feed_phaseMatrix = [1, 1; 1, 1];
default_feed_ampliMatrix = [1, 1; 1, 1];
% (Physical constants go here in index)
% Further params
default_J          = 2.5*90*(1e-6);     % microAmps -> Amps
                                        % 2.5 * Threshold.
                                        % Threshold from Redlich paper,
                                        % *2.5 from ott10
default_alpha_par  = 0;                 % Linewidth Enhancement Factor
default_omega1     = 0;                 % Rotational parameter for DDEBIF
default_omega2     = 0;                 % Rotational parameter for DDEBIF

%% inputParser
% This section allows the user to set parameter values by using the syntax:
% 'parameter', value
p = inputParser;
p.KeepUnmatched = true;
p.PartialMatching = false;

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
addParameter(p,'omega1', default_omega1);
addParameter(p,'omega2', default_omega2);

% Parse inputs/options.
parse(p,varargin{:})



end