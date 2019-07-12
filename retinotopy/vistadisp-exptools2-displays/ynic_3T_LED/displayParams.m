function params = displayParams

% Critical parameters
%%%% scanner set up
params.numPixels = [1920 1080];

%%% debug/home/linux set up
%%%params.numPixels = [1680 1050] ; %%%%;
% % params.numPixels = [1024 768];

%%% Scanner set up
params.dimensions =  [40 23];%%%[40 23]; %%%%

%%%% debug/home/linux set up
%%%%params.dimensions = [47 29.5]; %cm  %%


params.distance = [57]; %cm
params.frameRate = 120;
% % params.frameRate = 60;
params.cmapDepth = 8;
params.screenNumber = 0;

% Descriptive parameters
params.computerName = 'niwin402';
params.monitor = 'PROPixx_MRI_Proj';
params.card = 'AMD Radeon';
params.position = '3T scanner 16 channel setup';

