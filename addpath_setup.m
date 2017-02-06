function [  ] = addpath_setup(  )
%This function adds the relevant folders to matlab path.

present_working_directory = pwd;
addpath(present_working_directory)
addpath(strcat(present_working_directory,'/functions/'))
addpath(strcat(present_working_directory,'/functions/numerical_systems/'))
addpath(strcat(present_working_directory,'/functions/rot'))
addpath(strcat(present_working_directory,'/scripts/'))
addpath('~/dde_biftool_v3.1.1/ddebiftool/'); 
addpath('~/dde_biftool_v3.1.1/ddebiftool_extra_psol/');
addpath('~/dde_biftool_v3.1.1/ddebiftool_utilities/');
addpath('~/dde_biftool_v3.1.1/ddebiftool_extra_rotsym');

end

