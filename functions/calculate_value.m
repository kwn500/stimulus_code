function [change_value, target_value, ii, movement, max] = calculate_value(tTest, reference, feature_task)

% specifies change in feature and direction of change in target radial
% frequency pattern stimulus. 

% specify the maximum extent of change in each stimulus. 
if isequal(feature_task, 'Orientation')
    max = pi/8; movement = 'clockwise';
elseif isequal(feature_task, 'Contrast')
    max = 0.45; movement = 'lower contrast';
elseif isequal(feature_task, 'RAM')
    max = 0.1; movement = 'smoother';
end

ii = rand();

% convert threshold value to radians. 
change_value=10^tTest;    
change_value=min(change_value,max); 

% alter the direction of change on approximaely 50% of trials.
if ii > 0.5
    change_value = -change_value;
    if isequal(feature_task,'Orientation')
        movement = 'anticlockwise';
    elseif isequal(feature_task, 'Contrast')
        movement = 'higher contrast';
    elseif isequal(feature_task, 'RAM')
        movement = 'spikier';
    end
end
target_value = reference - change_value;
end
%% --------------------------------------------------------------------- %%