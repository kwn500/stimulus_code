function [event_trigger] = get_event_trigger(trial_data,...
    ref_data, thisTrial, thisBlock)

% specifies trigger value to send to MEG for each event on the basis of the
% combination of stimulus changes between the target and the reference.

if trial_data{thisBlock}(1,thisTrial) == ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) == ref_data.C...
        && trial_data{thisBlock}(3,thisTrial) == ref_data.A
    event_trigger = 200; % no change
    
elseif trial_data{thisBlock}(1,thisTrial) ~= ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) == ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) == ref_data.A
    event_trigger = 201; % orientation change
    
elseif trial_data{thisBlock}(1,thisTrial) == ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) ~= ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) == ref_data.A
    event_trigger = 202; % contrast change
    
elseif trial_data{thisBlock}(1,thisTrial) == ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) == ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) ~= ref_data.A
    event_trigger = 203; % shape change
    
elseif trial_data{thisBlock}(1,thisTrial) ~= ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) ~= ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) == ref_data.A
    event_trigger = 204; % orientation & contrast change
    
elseif trial_data{thisBlock}(1,thisTrial) ~= ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) == ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) ~= ref_data.A
    event_trigger = 205; % orientation & shape change
    
elseif trial_data{thisBlock}(1,thisTrial) == ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) ~= ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) ~= ref_data.A
    event_trigger = 206; % contrast and shape change
    
elseif trial_data{thisBlock}(1,thisTrial) ~= ref_data.phi &&...
        trial_data{thisBlock}(2,thisTrial) ~= ref_data.C &&...
        trial_data{thisBlock}(3,thisTrial) ~= ref_data.A
    event_trigger = 207; % all change.
end
end

%% --------------------------------------------------------------------- %%