%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXECUTABLE

addpath('./Fourier_Optics_Tools');

clear all

% Generate set of parameters
Set_Parameters;

%%%%%

Study_name = '';
Filename_parameters = ['/home/*/' Study_name, '_param.mat'];
Filename_results = ['/home/*/',Study_name, '_results.mat'];

% Save current parameters of the study
save(Filename_parameters);

% Load parameters of the study
load(Filename_parameters);
Parameters = load(Filename_parameters);

% Run the simulator
[Results] = Coronagraph_model(Parameters, Filename_results)
