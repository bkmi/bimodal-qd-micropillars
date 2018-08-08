%% Bifurcation continuations in the phase direction
%  axes = feedback amp, feedback phase, intensity
%  found in 'pretty_phases.m' data in 'rephases.mat'

clear;
load(strcat(specific_bif_data_dir(), ...
    'phase_continuations.mat'))

% Plot all relevant continuations at all feedback amplitudes.
f = figure;
for i = 1:numel(rephases)
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations')
view([0,1,0])

% Plot continuations at just above zero feedback amplitude.
f = figure;
for i = 1:3
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations at Feedback Amp = 0')
view([0,1,0])
saveas(f, strcat(data_directory(), 'phases_fbamp_0.png'))

% Plot continuations at feedback amplitude = 0.2
f = figure;
for i = 4:6
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations at Feedback Amp = 0.2')
view([0,1,0])
saveas(f, strcat(data_directory(), 'phases_fbamp_02.png'))

% Plot continuations at feedback amplitude = 0.34
f = figure;
for i = 7:8
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations at Feedback Amp = 0.34')
view([0,1,0])
saveas(f, strcat(data_directory(), 'phases_fbamp_034.png'))

% Plot continuations at feedback amplitude = 0.36
f = figure;
for i = 9
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations at Feedback Amp = 0.36')
view([0,1,0])
saveas(f, strcat(data_directory(), 'phases_fbamp_036.png'))

% Plot continuations at feedback amplitude = 0.38
f = figure;
for i = 10
    plot_branch3( rephases{i}, ...
        param, ...
        'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
        'nunst_color', rephases{i}.nunst, ...
        'add_2_gcf', 1)
end
zlabel('Strong Field Intensity')
title('Phase Continuations at Feedback Amp = 0.38')
view([0,1,0])
saveas(f, strcat(data_directory(), 'phases_fbamp_038.png'))
