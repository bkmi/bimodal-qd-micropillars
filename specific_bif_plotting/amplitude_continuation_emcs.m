%% Bifurcation continuations in the amplitude direction (ECMs)
%  axes = feedback amp, feedback phase, intensity
%  found in 'ecms.m' data in 'ecms.mat'

clear;
load(strcat(specific_bif_data_dir(), ...
    'amplitude_continuation_emcs_pruned.mat'))

% plot all emcs as phase = 0
figure;
shape = size(famp_cur560_phase0);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phase0{i,j})
            plot_branch3( famp_cur560_phase0{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phase0{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end

% plot all emcs as phase = pi/2
figure;
shape = size(famp_cur560_phasePi2);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phasePi2{i,j})
            plot_branch3( famp_cur560_phasePi2{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phasePi2{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end


% plot all emcs as phase = pi/4
figure;
shape = size(famp_cur560_phasePi4);
for i = 1:shape(1)
    for j = 1:shape(2)
        if ~isempty(famp_cur560_phasePi4{i,j})
            plot_branch3( famp_cur560_phasePi4{i,j}, ...
                param, ...
                'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
                'nunst_color', famp_cur560_phasePi4{i,j}.nunst, ...
                'add_2_gcf', 1)
        end
    end
end
