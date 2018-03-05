function [  ] = addpath_setup(  )
%This function adds the relevant folders to matlab path.

%% My Files
present_working_directory = main_directory();
addpath(present_working_directory)

%funcs
addpath(strcat(present_working_directory,'/functions/'))
addpath(strcat(present_working_directory,'/functions/numerical_systems/'))
addpath(strcat(present_working_directory,'/functions/timeSeries_data_gen_tools/'))
addpath(strcat(present_working_directory,'/functions/bifcont_data_gen_tools/'))
addpath(strcat(present_working_directory,'/functions/data_analysis_tools/'))

%scripts
addpath(strcat(present_working_directory,'/scripts/'))
addpath(strcat(present_working_directory,'/scripts/feedbackTree/'))
addpath(strcat(present_working_directory,'/scripts/simple_run/'))
addpath(strcat(present_working_directory,'/scripts/IOcurveScripts/'))
addpath(strcat(present_working_directory,'/scripts/prez_scripts/'))
addpath(strcat(present_working_directory,'/scripts/bifurcationSelfFB/'))
addpath(strcat(present_working_directory,'/scripts/report/'))


%% DDE-BIFTOOL
addpath(strcat(present_working_directory,'/ddebiftool_multirot/')) % Multi_rot
addpath(strcat(present_working_directory,'/dde_biftool_v3.1.1/ddebiftool/')); 
addpath(strcat(present_working_directory,'/dde_biftool_v3.1.1/ddebiftool_extra_psol/'));
addpath(strcat(present_working_directory,'/dde_biftool_v3.1.1/ddebiftool_utilities/'));
addpath(strcat(present_working_directory,'/dde_biftool_v3.1.1/ddebiftool_extra_rotsym/'));

%% Plotting
addpath(strcat(present_working_directory,'/BrewerMap/'))

end