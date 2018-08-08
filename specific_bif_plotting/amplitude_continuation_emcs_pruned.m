%% Bifurcation continuations in the amplitude direction (ECMs)
%  axes = feedback amp, feedback phase, intensity
%  found in 'ecms.m' data in 'ecms.mat'

clear;
load(strcat(specific_bif_data_dir(), ...
    'amplitude_continuation_emcs_pruned.mat'))

% plot all emcs as phase = 0
f = figure;
for i = 1:numel(famp_cur560_phase0_pruned)
    if ~isempty(famp_cur560_phase0_pruned{i,1})
        plot_branch3( famp_cur560_phase0_pruned{i,1}, ...
            param, ...
            'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
            'nunst_color', famp_cur560_phase0_pruned{i,1}.nunst, ...
            'add_2_gcf', 1)
    end
end
zlabel('Strong Field Intensity')
title('Amplitude Continuations at \phi_{fb} = 0')
view([1,0,0])
saveas(f, strcat(data_directory(), 'emc_0.png'))

% plot all emcs as phase = pi/2
f = figure;
for i = 1:numel(famp_cur560_phasePi2_pruned)
    if ~isempty(famp_cur560_phasePi2_pruned{i,1})
        plot_branch3( famp_cur560_phasePi2_pruned{i,1}, ...
            param, ...
            'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
            'nunst_color', famp_cur560_phasePi2_pruned{i,1}.nunst, ...
            'add_2_gcf', 1)
    end
end
zlabel('Strong Field Intensity')
title('Amplitude Continuations at \phi_{fb} = \pi/2')
view([1,0,0])
saveas(f, strcat(data_directory(), 'emc_pi2.png'))

% plot all emcs as phase = pi/4
f = figure;
for i = 1:numel(famp_cur560_phasePi4_pruned)
    if ~isempty(famp_cur560_phasePi4_pruned{i,1})
        plot_branch3( famp_cur560_phasePi4_pruned{i,1}, ...
            param, ...
            'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
            'nunst_color', famp_cur560_phasePi4_pruned{i,1}.nunst, ...
            'add_2_gcf', 1)
    end
end
zlabel('Strong Field Intensity')
title('Amplitude Continuations at \phi_{fb} = \pi/4')
view([1,0,0])
saveas(f, strcat(data_directory(), 'emc_pi4.png'))


% %%Pruning
% famp_cur560_phase0_pruned = cell(2,1);
% famp_cur560_phase0_pruned{1,1} = famp_cur560_phase0{1,1};
% famp_cur560_phase0_pruned{2,1} = famp_cur560_phase0{1,4};
% 
% famp_cur560_phasePi2_pruned = cell(4, 1);
% famp_cur560_phasePi2_pruned{1,1} = famp_cur560_phasePi2{1,1};
% famp_cur560_phasePi2_pruned{2,1} = famp_cur560_phasePi2{1,3};
% famp_cur560_phasePi2_pruned{3,1} = famp_cur560_phasePi2{1,4};
% famp_cur560_phasePi2_pruned{4,1} = famp_cur560_phasePi2{1,5};
% 
% famp_cur560_phasePi4_pruned = cell(4,1);
% famp_cur560_phasePi4_pruned{1,1} = famp_cur560_phasePi4{1,1};
% famp_cur560_phasePi4_pruned{2,1} = famp_cur560_phasePi4{1,3};
% famp_cur560_phasePi4_pruned{3,1} = famp_cur560_phasePi4{1,4};
% famp_cur560_phasePi4_pruned{4,1} = famp_cur560_phasePi4{1,5};
% 
% %% change
% temp = famp_cur560_phasePi4_pruned;
% ind = 4;
% 
% temp{ind, 1} = bcont( funcs, ...
%                   temp{ind, 1}, ...
%                   800, ...
%                   800 );
% temp{ind, 1} = branch_stability(funcs, temp{ind, 1});
% 
% figure;
% for i = 1:numel(temp)
%     if ~isempty(temp{i,1})
%         plot_branch3( temp{i,1}, ...
%             param, ...
%             'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
%             'nunst_color', temp{i,1}.nunst, ...
%             'add_2_gcf', 1)
%     end
% end
% 
% %% save
% 
% famp_cur560_phasePi4_pruned = temp;