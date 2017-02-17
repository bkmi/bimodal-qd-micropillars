function [  ] = addpath_setup(  )
%This function adds the relevant folders to matlab path.

% Add multi_rot_sym directly from location on my own personal drive. If you
% aren't me, and are hoping to use mult_rot then you will need to direct
% this add path to the multirot directory in DDE-BIFTOOL.
% addpath(strcat(present_working_directory,'/functions/rot'))
addpath('/home/bkmiller/qd-micropillar-laser-project/ddebiftool_multirot')

present_working_directory = pwd;
addpath(present_working_directory)
addpath(strcat(present_working_directory,'/functions/'))
addpath(strcat(present_working_directory,'/functions/numerical_systems/'))
addpath(strcat(present_working_directory,'/functions/timeSeries_data_gen_tools'))
addpath(strcat(present_working_directory,'/functions/bifcont_data_gen_tools'))
addpath(strcat(present_working_directory,'/functions/data_analysis_tools'))
addpath(strcat(present_working_directory,'/scripts/'))
addpath('~/dde_biftool_v3.1.1/ddebiftool/'); 
addpath('~/dde_biftool_v3.1.1/ddebiftool_extra_psol/');
addpath('~/dde_biftool_v3.1.1/ddebiftool_utilities/');
addpath('~/dde_biftool_v3.1.1/ddebiftool_extra_rotsym');

end

