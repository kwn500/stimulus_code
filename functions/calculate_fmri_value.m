function [o_final_target, o_tracker, o_direction, c_final_target, c_tracker,...
    c_direction, r_final_target, r_tracker, r_direction]...
    = calculate_fmri_value(o_target, c_target, r_target)

% reads in stimulus attributes and calculates direction of change and
% correct response.

o_final_target = o_target; o_tracker = 0; o_direction = 'no_change';
c_final_target = c_target; c_tracker = 0; c_direction = 'no_change';
r_final_target = r_target; r_tracker = 0; r_direction = 'no_change';

C = 0.5;	% Contrast.
A = 0.5;    % Radial modulation amplitude (A < 1).
phi = 0;    % Phi - angular rotation.

% If target orientation is different from the reference value:
if isequal(o_target, phi) == 0
    o_tracker = 1; % Mark this as an orientation change.
    o_target = -o_target; % Alter the direction of change. 
    o_direction = 'anticlockwise'; % Store direction in relevant cell array. 
    o_final_target = phi - o_target; % Set orientation below the reference value.      
end

% repeat the same process for contrast and shape. 
if isequal(c_target, C) == 0
    c_tracker = 1;
    c_target = -c_target;
    c_direction = 'higher contrast';
    c_final_target = C - c_target;
end

if isequal(r_target, A) == 0
    r_tracker = 1;
    r_target = -r_target;
    r_direction = 'spikier';
    r_final_target = A - r_target;
end
end
%% --------------------------------------------------------------------- %%