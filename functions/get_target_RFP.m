function [stim_draw] = get_target_RFP(stim, stim_changes)

% this function takes in a struct containing colour-specific LMS cone 
% modulations and contrasts, and a 1 x 3 logical with numbers corresponding
% to the presence/absence of changes in particular visual features 
% (orientation, contrast, shape). 

% 1 corresponds to no-change (i.e the feature-specific values remain the
% same, as the reference). 
% 0 corresponds to a change in feature-specific values from their reference. 

% this function finds the pre-generated colour-RFP stimulus matching the
% changes specified in the 1x3 logical, and outputs the image matrix to be 
% presented in the psychtoolbox experimental paradigm. 

if isequal(stim_changes,[1 1 1]) % no change. 
    stim_draw = stim.no_change_colour;
elseif isequal(stim_changes, [0 1 1]) % orientation change. 
    stim_draw = stim.o_colour;
elseif isequal(stim_changes, [1 0 1]) % contrast change. 
    stim_draw = stim.c_colour;
elseif isequal(stim_changes, [1 1 0]) % shape change. 
    stim_draw = stim.s_colour;
elseif isequal(stim_changes, [0 0 1]) % orientation & contrast change. 
    stim_draw = stim.oc_colour;
elseif isequal(stim_changes, [0 1 0]) % orientation & shape change. 
    stim_draw = stim.os_colour;
elseif isequal(stim_changes, [1 0 0]) % contrast & shape change. 
    stim_draw = stim.cs_colour;
elseif isequal(stim_changes, [0 0 0]) % orientation, contrast & shape change. 
    stim_draw = stim.ocs_colour;
end
end

