%% Self Coupled, Zero Feedback Phase Offset Tree Creator
% We begin by creating a tree structure. This structure organizes our saved
% data and our thinking about the bifurcations.
%
% The (trivial) trunk of the tree is studying bimodal micropillars with
% self coupled feedback.
%
% The first branch is how the feedback phase is varied. This particular
% script considers just one feedback phase. I.e. the strong feedback phase
% and the weak feedback phase behave in the same way (linearly) without any
% offset. Essentially, there is only one feedback phase variable.
%
% The second branches are at the bistability for high injection currents.
% There is the case where the weak field dominates and the case where the
% strong field dominates. We consider bifurcation diagrams beginning at
% with one of these two options. The current where weakDom takes over
% strongDom is somewhere around 180mA. 
%
% Next, there is a branch concerning which feedback time, tau_fb, we use.
% Experimentalists used five feedback times. I will use the same.
%
% That gives us the following directory structure layout:
% 
% zeroPhaseOffset (with J=560 uA)
%   |-weakDom
%   |   -alpha1
%   |   |   -tau_fb1
%   |   |   -tau_fb2
%   |   |   -...
%   |   -alpha2
%   |   |   -tau_fb1
%   |   |   -...
%   |   -...
%   |
%   |-strongDom
%   |   -alpha1
%   |   |   -tau_fb1
%   |   |   -tau_fb2
%   |   |   -...
%   |   -alpha2
%   |   |   -tau_fb1
%   |   |   -...
%   |   -...

%% System Parameters
% Feedback Parameters
feedPhaseMat = [1, 0; 0, 1];
feedAmpMat = [1, 0; 0, 1];
feedTimesShort = [0.83, 1.17, 1.5, 1.83, 2.17, 2.5, 2.83, 3.17]; %Short Cavity 
feedTimesLong  = [2.93,4.67,5.60,7.40,9.33]; %Long Cavity
tau_fbArray = [0.83, 1.5, 2.93, 5.60, 9.33];

% Laser Parameters
J = 560e-6;
alphas = [0, 3];

% directory creation
datadir = '/home/bkmiller/qd-micropillar-laser-project/data_bimodal-qd-micropillars/';
trunkdir = strcat(datadir,'zeroPhaseOffsetJ=',num2str(J*1e6),'uA/');

mkdir(trunkdir);

for i = 1:numel(alphas)
    weakdir = strcat(trunkdir,'weakDomAlpha=',num2str(alphas(i)),'/');
    mkdir(weakdir)
    for j = 1:numel(tau_fbArray)
        tempdir = strcat(weakdir,'tau_fb=',num2str(tau_fbArray(j)),'ns/');
        
        setup_params_nonDim_CnstCplRatio(...
            'save',1, ...
            'J', J,...
            'alpha_par',alphas(i), ...
            'tau_fb',tau_fbArray(j),...
            'feed_ampli',0.15, ...
            'feed_ampliMatrix', feedAmpMat, ...
            'feed_phase',0, ...
            'feed_phaseMatrix', feedPhaseMat, ...
            'clear',0, ...
            'populate_wrkspc', 0, ...
            'datadir_specific',tempdir);
    end
    
    strongdir = strcat(trunkdir,'strongDomAlpha=',num2str(alphas(i)),'/');
    mkdir(strongdir)
    for j = 1:numel(tau_fbArray)
        tempdir = strcat(strongdir,'tau_fb=',num2str(tau_fbArray(j)),'ns/');
        
        setup_params_nonDim_CnstCplRatio(...
            'save',1, ...
            'J', J,...
            'alpha_par',alphas(i), ...
            'tau_fb',tau_fbArray(j),...
            'feed_ampli',0.15, ...
            'feed_ampliMatrix', feedAmpMat, ...
            'feed_phase',0, ...
            'feed_phaseMatrix', feedPhaseMat, ...
            'clear',0, ...
            'populate_wrkspc', 0, ...
            'datadir_specific',tempdir);
    end
end


