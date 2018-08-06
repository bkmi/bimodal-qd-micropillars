%% Mode seperation
%  axes = feedback amp, feedback phase, intensity
%  found in 'vfolds.m' (I think) data in 'folds_main.mat'
%  found in 'fold_mode_sep.m' data in 'kappaw45b.mat'
%  found in 'fold_mode_sep.m' data in 'kappaw47b.mat'

% Standard seperation (4.1e10)
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_main.mat'))

figure;
for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'k', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end


hold on

% bigger seperation (4.5e10)
clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_kappaw45.mat'))

for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'b', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

% biggest seperation (4.7e10)

clear;
load(strcat(specific_bif_data_dir(), ...
    'folds_kappaw47.mat'))

for i = 1:numel(folds)
    plot_branch3( folds{i}, ...
    param, ...
    'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
    'color', 'r', ...
    'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
    'add_2_gcf', 1)
end

% %% Find the good bifurcations...
% clear;
% load(strcat(data_directory(), ...
%     '\specific_bifurcations\kappaw47b.mat'))
% % 1, 2, 3, 4, 5, 7 for kappaw45
% % have to get the positive ones for kappa47
% % pos: 1, 2, 3, 7
% % pos shift: 4, 5, 6, 7, 8, 10
% 
% %% shift phase
% for i = 1:numel(folds)
%     for j = 1:numel(folds{i}.point)
%         val = folds{i}.point(j).parameter(param.feed_phase.index);
%         folds{i}.point(j).parameter(param.feed_phase.index) = val + 2 * pi;
%     end
% end
% 
% %%
% hold on
% for i = 10
%     plot_branch3( pfolds{i}, ...
%     param, ...
%     'axes_indParam', {param.feed_phase.index, param.feed_ampli.index, 'x1'}, ...
%     'color', 'r', ...
%     'PlotStyle', { 'LineStyle', '-', 'Marker', '.' }, ...
%     'add_2_gcf', 1)
% end
% 
% %%
% posinds = [];
% for i = 1:numel(folds)
%     if folds{i}.point(1).x(1) > 0
%         posinds(end+1) = i;
%     end
% end
% pfolds = cell(numel(posinds), 1);
% 
% counter = 1;
% for i = posinds
%     pfolds{counter} = folds{i};
%     counter = counter + 1;
% end
% % 1, 
% 
% %%
% sfolds = cell(6,1);
% counter = 1;
% for i = [4, 5, 6, 7, 8, 10]
%     sfolds{counter} = pfolds{i};
%     counter = counter + 1;
% end
% 
% for i = 1:numel(sfolds)
%     sfolds{i} = bcont( foldfuncs, sfolds{i}, ...
%                             4000, 4000 );
% end
% 
% %%
% if 1
%     save(strcat(specific_bif_data_dir(), 'folds_kappaw47o'), ...
%     'funcs', ...
%     'foldfuncs', ...
%     'famp_phase0', ...
%     'phase_famp05', ...
%     'folds', ...
%     'param', ...
%     'st_amp', ...
%     'st_bif', ...
%     'st_phase', ...
%     'opt_inputs')
% end