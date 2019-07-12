function [img, instructions, reference] = get_feature_task_info(feature_task, base_directory, combine, toggle)

% if this is a psychophysics experiment (i.e a feature-specific threshold
% staircase):
if isequal(toggle, 'psychophysics')
    
    % select the feature-specific fixation letter (we're using .png files
    % here with the same mid-gray background as psychtoolbox to prevent
    % issues from the previous experiment with the letter size/thickness
    % changing across trials. 
    img = imread(strcat(base_directory,'/fixation_letters', combine,...
        sprintf('%s_affinity.png', lower(feature_task))));
    
    % also specify the feature-specific introduction. 
    line1 = sprintf('%s Discrimination Task', feature_task);
    line4 = '\n\n Press the spacebar to continue';
    
    % if this is an attend-orientation staircase:
    if isequal(lower(feature_task), 'orientation') 
        reference = 0; % specify the default reference value.
        
        % specify the orientation-specific instructions (i.e.
        % clockwise/anticlockwise). 
        line2 = sprintf('\n\n If the second is rotated ANTICLOCKWISE with respect the first, press "N"');
        line3 = sprintf('\n\n If the second is rotated CLOCKWISE with respect to the first, press "U"');
        
    % repeat the same process for attend-contrast and attend-shape staircases.     
    elseif isequal(lower(feature_task), 'contrast') 
        reference = 1;
        line2 = sprintf('\n\n If the second is HIGHER contrast with respect the first, press "U"');
        line3 = sprintf('\n\n If the second is LOWER contrast with respect the first, press "N"');
        
    elseif isequal(lower(feature_task), 'shape')
        reference = 0.5;
        line2 = sprintf('\n\n If the second is SPIKIER (higher RAM) than the first, press "U"');
        line3 = sprintf('\n\n If the second is SMOOTHER (lower RAM) than the first, press "N"');
    end
    
    instructions = [line1, line2, line3, line4]; % combine the instructions in order.

% else, if this is a contrast detection task:
elseif isequal(toggle, 'contrast_detection')
    
    % get the relevant contrast fixation letter. 
    img = imread(strcat(base_directory,'/fixation_letters', combine, 'contrast_affinity.png'));
     
    reference = 1; % specify the default reference value (unused). 
    
    % specify the contrast detection specific instructions and keypresses.
    line1 = 'Contrast Discrimination Task'; 
    line2 = '\n\n Press A if you thought the stimulus was in the first interval';
    line3 = '\n\n Press L if you thought the stimulus was in the second interval';
    line4 = '\n\n Press the spacebar to continue';
    
    instructions = [line1, line2, line3, line4]; % combine the instructions in order.  
end

