function runLocalRetino(setupname,site, WaitMRITrigger)

cd('C:\Users\Ryan Maloney\Documents\Kirstie\retinotopy')

setupname = [setupname '.mat'];

if ~exist(setupname, 'file')
    disp('File does not exist, moron')
    return
end

dat = load(setupname);

% Set up the parameters
params = setRetinotopyParams(dat.params.experiment, dat.params);

% set response device
params = setRetinotopyDevices(params);
% keyboard
% go
doRetinotopyScan_Kirstie(params, site, WaitMRITrigger);

% Flush any events
FlushEvents;

end