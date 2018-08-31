%% Hopf stuck? Perhaps there are degenerate hopf solutions? 

clear;

directory = 'C:/Users/Public/matlab/bimodal-qd-micropillars';
addpath(strcat(directory,'/ddebiftool_multirot/')) % Multi_rot
addpath(strcat(directory,'/dde_biftool_v3.1.1/ddebiftool/')); 
addpath(strcat(directory,'/dde_biftool_v3.1.1/ddebiftool_extra_psol/'));
addpath(strcat(directory,'/dde_biftool_v3.1.1/ddebiftool_utilities/'));
addpath(strcat(directory,'/dde_biftool_v3.1.1/ddebiftool_extra_rotsym/'));

load('conts.mat')

amp = cell(2,1);
for i = 3:4
    amp{i - 2} = famp_cur560_phasePi2_pruned{i};
end

hopf = cell(2,1);
hopf_inds = [1, 3];

for i = 2 % 1:numel(amp)
    [hopf{i}, ~]=SetupHopf( ...
        funcs, ...
        amp{i}, ...
        amp{i}.indHopf(hopf_inds(i)), ...
        'contpar', ...
        [param.feed_phase.index, param.feed_ampli.index, param.omega1.index], ...
        'dir', param.feed_phase.index, ...
        opt_inputs{:},...
        st_bif_cur_amp{:});
    
    [hopf{i},~,~,~] = br_contn(funcs, hopf{i}, 250);
    hopf{i} = br_rvers(hopf{i});
    [hopf{i},~,~,~] = br_contn(funcs, hopf{i}, 250);
end


f = figure;
for i = 1:2
    plot_branch3( amp{i,1}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', amp{i,1}.nunst, ...
        'add_2_gcf', 1)
end

plot_branch3( hopf{2}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'c', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)