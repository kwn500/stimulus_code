%% --------- RFP 2AFC SELECTIVE/NON-SELECTIVE ATTENTION TASK ----------- %%

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

    % ORIENTATION, CONTRAST, RADIAL AMPLITUDE MODULATION & All 3 %  
       
% Participant must judge whether target RFP is different in a cued feature, 
% (selective) or different in any feature with respect to a reference RFP 
% (non-selective attention). Participants must response with 'A' or 'L' 
% button presses indicating presence/absence of change. 

%% ----------------- PRE-ALLOCATING VARIABLES FOR SPEED ---------------- %%   
trial_number = cell(1); conditions = cell(1); pp_orientation = cell(1);   
reference_orientation = cell(1); target_orientation = cell(1); orientation_direction = cell(1);  
orientation_change_tracker = cell(1); pp_contrast = cell(1); reference_contrast = cell(1);
target_contrast = cell(1); contrast_direction = cell(1); contrast_change_tracker = cell(1); 
pp_ram = cell(1); reference_ram = cell(1); target_ram = cell(1); ram_direction = cell(1); 
ram_change_tracker = cell(1); key_pressed = cell(1); response_made = cell(1);    
marked_responses = cell(1); reaction_time = cell(1); change_count = cell(1);
change_type = cell(1); conjunction = cell(1); response_required = cell(1); 
block_conditions = cell(1); orientation_seq = cell(1); contrast_seq = cell(1);  
ram_seq = cell(1); reference_RFPs = cell(1); reference_RFP_pres = cell(1);
task_condition = cell(1); trial_offsets = cell(1); trial_fixations = cell(1);
trial_fixation_pres = cell(1); target_RFPs = cell(1); target_RFP_pres = cell(1);
response_fixations = cell(1); response_fixation_pres = cell(1); 
trial_ends = cell(1); total_trial_starts_block = cell(1); 
total_trial_starts_scan = cell(1); total_trial_ends_block = cell(1);
total_trial_ends_scan = cell(1);

%% -------------------------- GENERAL SETUP ---------------------------- %%
sca; close all; clearvars; KbName('UnifyKeyNames');  jheapcl;
start_dir = pwd;

% CHECK THESE SETTINGS MATCH PIXELS PER DEGREE SETTINGS IN RFP_GENERATOR. 
% Calculate number of pixels per degree to control later stimulus sizes. 
display_width = 52.4;       % Monitor width (cm).
display_horRes = 1920;      % Horizontal resolution of monitor (pixels).
display_viewDist = 57;      % Participant viewing distance (cm).
display_ppd = display_horRes/(2*atand((display_width/2)/display_viewDist)); % Calculates pixels/°.

% Create and present dialog box that gathers user input.
prompt = {'Enter Participant Number:','Enter Scan Number:', 'Enter Trial Type:'};
dlg_title = 'Participant Information';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

% Set participant id, scan number and feature task variables from GUI input. 
pid = str2double(answer{1}); scan_number = str2double(answer{2}); trial_type = strcat(lower(answer{3}),'_trials');

% Toggle trial type; 'equal_trials' vs 'unequal_trials'.
% Equal = equal trials to correctly respond to whether selective or non-selective (less changes in non-selective). 
% Unequal = less trials to correctly respond to in selective than non-selective, but same amount of changes across...
    % all features in both conditions (recommended). 
if isequal(trial_type,'unequal_trials')   
    short_trial_type = 'UT';            % Short trial type variable for data output ('UT' & 'ET');
elseif isequal(trial_type,'equal_trials')
    short_trial_type = 'ET';
end

% Vital experimental information - check this is correct before running!. 
trials_desired = 10;                % Number of trials.
task = 'T';                         % Type of task - 'T' = temporal, 'S' = spatial. 
threshold_repeat = 3;               % Number of threshold runs (ensures you load in correct thresholds). 
condition_task = 'temporal';        % Task condition string (for data output)- 'temporal'/'spatial'. 
block_repeats = 2;                  % Number of block repeats (one scan).
prob_change_1 = 0.8;                % Probability of change in stimulus feature (0.8 = 20% likelihood change) (unequal & equal trials (selective only)).
prob_change_2 = (1 - prob_change_1) / 3;  % Probability of change in stimulus feature for non-selective condition if equal trials desired. 
prob_change_2 = 1 - prob_change_2;
stimulus_presentation_time = 0.2;   % Stimulus presentation time (s). 
response_presentation_time = 0.8;   % Response period duration (s). 
total_trial_length = 1.5;           % Total trial duration (s).  
    
% Sets up directories for data storage and access to functions. 
% root_directory = ('/Users/Kirstie/Documents/radial_frequency_code');    % Top-level directory.
root_directory = ('\Users\PS-RMLab\Documents\Kirstie\updated');    % Top-level directory- Shuttle.

orientation_threshold = 0; contrast_threshold = 0; ram_threshold = 0; 

% Access specific participants unique experimental thresholds directory. 
% file_path = strcat(root_directory,sprintf('/data_output/experimental_thresholds/%d_experimental_thresholds',pid)); cd(file_path);  
if pid < 10
    file_path = strcat(root_directory,sprintf('\\data_output\\experimental_thresholds\\0%d_experimental_thresholds',pid)); cd(file_path); 
else
    file_path = strcat(root_directory,sprintf('\\data_output\\experimental_thresholds\\%d_experimental_thresholds',pid)); cd(file_path);  
end

for i = 1:threshold_repeat    % for one to the number of threshold repeats each participant completed:
    
    if pid < 10
        load(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat', condition_task, i, pid));                    % Load in participants threshold values.
    else
        load(sprintf('thresholds_%s_scan_0%d_pp_%d.mat', condition_task, i, pid));                    % Load in participants threshold values.
    end
    
    % Add all thresholds for each feature to single variable (running total). 
    orientation_threshold = orientation_threshold + orientationT;   contrast_threshold = contrast_threshold + contrastT; ram_threshold = ram_threshold + RAMT;   
end

% set threshold variables to average ot thresholds from number of threshold scan repeats- then set to double threshold. 
orientation_threshold = orientation_threshold / threshold_repeat; orientation_threshold = orientation_threshold * 2;
contrast_threshold = contrast_threshold / threshold_repeat; contrast_threshold = contrast_threshold * 2;
ram_threshold = ram_threshold / threshold_repeat;   ram_threshold = ram_threshold * 2;

% cd(strcat(root_directory,'/supporting_files'));                                                     % Return to relevant function directory. 
base_directory = strcat(root_directory,'\supporting_files');             % Path to relevant functions- Shuttle.
cd(base_directory); 

% data_storage = strcat(root_directory,'/data_output/selective_vs_non-selective');                    % Set data storage directory. 
data_storage = strcat(root_directory,'\data_output\selective_vs_non_selective');                    % Set data storage directory. 

% General experiment 'housekeeping'. 
block_label = {'ORIENTATION' 'CONTRAST' 'SHAPE'};                   % Specify block labels. 
block_tracker = [1 2 3];                                            % Specify block numbers.
fixation_label = {'O' 'C' 'S' 'N'};                                 % Specify fixation letters. 
fixation_colours = {[1 0 0], [0 1 0], [0 0 1], [0.5 0 0.9]};        % Specify fixation colours. 
 
ii  = 1;                                                            % Establish counter variable. 
for i = 1:block_repeats                                            % For number of blocks specified:
    block_conditions{ii} = 1; block_conditions{ii+1} = 2; ii = ii+ 2;   % Add selective and non-selective conditions to block array. 
end

block_tracker = block_tracker(randperm(length(block_tracker)));          % Randomise selective condition order. 
block_conditions = block_conditions(randperm(length(block_conditions))); % Randomise order of blocks within scan. 

space = KbName('space'); escape = KbName('escape');  trigger = KbName('5%');        % Set keynames for potential response keypresses.   

PsychDefaultSetup(2); screens = Screen('Screens');  screenNumber = max(screens);    % Call default PTB settings & draw to max screen.
black = BlackIndex(screenNumber); white = WhiteIndex(screenNumber);                 % Define black (0) & white (1).

%% -------------------------- SCREEN SETUP------------------------------ %%

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, (128/255));    % Open gray screen window. 
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');     % Alpha-blending for anti-aliased lines.
HideCursor; ListenChar (2)                                                     % Hide mouse cursor & disable keypress output to the command window. 

[screenXpixels, screenYpixels] = Screen('WindowSize', window);                 % Window size (pixels). 
[xCenter, yCenter] = RectCenter(windowRect);                                   % Centre of screen (pixels). 

% No longer controls stimulus size (this is hard-coded in Rich's RFP_generator function), but controls position of spatial targets. 
baseRect = [0 0 500 500]; 

xCenter_left = xCenter - (5 * display_ppd) ; xCenter_right = xCenter + (5 * display_ppd) ; % Position RFPs 5° visual angle from centre.

% Set positions of target % reference RFP (reference - left, target - right).  
rect_ref_left = CenterRectOnPointd(baseRect, xCenter_left, yCenter);
rect_target_right = CenterRectOnPointd(baseRect, xCenter_right, yCenter); 

%% ----------------------- INTRODUCTION SCREEN ------------------------- %%

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial');    % Set text parameters.

% Present relevant instructions for selected condition. 
line1 = 'Detection Task'; 
if task == 'T'
    line2 = '\n\n Compare the first and second shapes.';
    line3 = '\n\n If the second stimulus is different from the first, press "A"';
    line4 = '\n If the second stimulus is the same as the first, press "L"'; 
elseif task == 'S'
    line2 = '\n\n Compare the left and right shapes.'; 
    line3 = '\n\n If the right stimulus is different from the left, press "A"';
    line4 = '\n If the right stimulus is the same as the left, press "L"'; 
end
line5 = '\n\n You will be cued to attend to different features in each block';
line6 = '\n Press 5 to continue';  

DrawFormattedText(window, [line1 line2 line3 line4 line5 line6], 'center', 'center', white);
Screen('Flip',window); [secs,KeyCode] = KbWait; % Draw text to screen and wait for keypress.  

% Wait for trigger response to progress. 
while KeyCode(trigger) == 0  % While trigger is not pressed.
    [secs,KeyCode] = KbWait;    % Wait for trigger response. 
    if KeyCode(escape)         % If trigger pressed:
        sca; ListenChar(); ShowCursor;  % Close all, re-enable keypresses to command window and make mouse cursor visible. 
    end
end

%% ----------------- RADIAL FREQUENCY PATTERN SETUP -------------------- %%

% Set parameters for RFP. 
C = 0.5;	% Contrast.  
A = 0.5;    % Radial modulation amplitude (A < 1).  
phi = 0;    % Phi - angular rotation. 
        
distances = RFP_generator(phi,C,A); % Create reference RFP. 
scn = Screen('MakeTexture', window, distances); % Convert reference RFP number matrix to texture.

%% ---------------------- STIMULUS PRESENTATION ------------------------ %%

% Set script to maximum processing priority. 
topPriorityLevel = MaxPriority(window); Priority(topPriorityLevel); 

i = 1; b = 1; c = 1;        % Establish counter variables. 
trial_end_in_scan = 0;      % Establish timing variables. 
trial_start_in_scan= 0;
trial_offset = 0;
total_trial_start_scan = GetSecs();

for j = 1:length(block_conditions)   % For 1 to the number of blocks desired:
    
    b_condition = block_conditions{j};  % Select block condition from randomised array.
    
    if isequal(trial_type, 'unequal_trials')   % If unequal trials have been specified:
        for a = 1:trials_desired                   % For number of trials desired:

            o_randomiser = rand(); c_randomiser = rand(); r_randomiser = rand();    % Establish random numbers for loop iterations. 
            orientation_seq{a} = phi; contrast_seq{a} = C; ram_seq{a} = A;          % Fill cell array position with default reference value. 

            if o_randomiser > prob_change_1                    % If random value is above 0.8 (20% chance).
                orientation_seq{a} = orientation_threshold;         % Add orientation value to orientation trial sequence.   
            end
            if c_randomiser > prob_change_1    
                contrast_seq{a} = contrast_threshold;
            end
            if r_randomiser > prob_change_1
                ram_seq{a} = ram_threshold;
            end
        end 
    elseif isequal(trial_type, 'equal_trials') % Else Iif unequal trials have been specified:
        if b_condition == 2                          % If condition is selective:
            for a = 1:trials_desired                   % For number of trials desired:

                o_randomiser = rand(); c_randomiser = rand(); r_randomiser = rand();    % Establish random numbers for loop iterations. 
                orientation_seq{a} = phi; contrast_seq{a} = C; ram_seq{a} = A;          % Fill cell array position with default reference value. 

                if o_randomiser > prob_change_1                    % If random value is above 0.8 (20% chance).
                    orientation_seq{a} = orientation_threshold;         % Add orientation value to orientation trial sequence.  
                end
                if c_randomiser > prob_change_1   
                    contrast_seq{a} = contrast_threshold;
                end
                if r_randomiser > prob_change_1
                    ram_seq{a} = ram_threshold;
                end
            end 
        elseif b_condition == 1                      % Else if condition is non-selective: 
            for a = 1:trials_desired                   % For number of trials desired:
                
                o_randomiser = rand(); c_randomiser = rand(); r_randomiser = rand();    % Establish random numbers for loop iterations. 
                orientation_seq{a} = phi; contrast_seq{a} = C; ram_seq{a} = A;          % Fill cell array position with default reference value. 
                
                if o_randomiser > prob_change_2                    % If random value is above 0.933 (6.33% chance).
                    orientation_seq{a} = orientation_threshold;         % Add orientation value to orientation trial sequence.  
                end
                if c_randomiser > prob_change_2  
                    contrast_seq{a} = contrast_threshold;
                end
                if r_randomiser > prob_change_2
                    ram_seq{a} = ram_threshold;
                end
            end
        end
    end
   
    % Shuffle trial order within these sequences.
    orientation_seq = orientation_seq(randperm(length(orientation_seq))); 
    contrast_seq = contrast_seq(randperm(length(contrast_seq)));  
    ram_seq = ram_seq(randperm(length(ram_seq)));  

    if isequal(b_condition,1)      % If the condition is non-selective, set specific values for certain parameters. 
        block_condition = 'Non-Selective'; specific_condition = 'All'; fixation = fixation_label{4}; 
        instructions = ('Respond to changes in ORIENTATION or CONTRAST or SHAPE'); 
    elseif isequal(b_condition,2) 	% Else if condition is selective: 
        if b == 4                      % Select each of the selective conditions from their randomised order in turn. 
            b = 1;                          % If more than three conditions are specified, conditions are re-randomised and sequence begins again.
            block_tracker = block_tracker(randperm(length(block_tracker))); % Randomise selective condition order.
        end
        final_condition = block_tracker(b); specific_condition = block_label{final_condition};                       % Select final_condition. 
        fixation = fixation_label{final_condition};             % Set relevant fixation letter/colour. 
        block_condition = 'Selective'; instructions = sprintf('Respond to changes only in %s', specific_condition);  % Set relevant instructions. 
        b = b + 1;                                                                                                   % Increment counter variable.    
    end
    
    % Present condition-specific instruction screen.
    Screen('TextSize', window, 25);  
    DrawFormattedText(window, sprintf(instructions), 'center', 'center', white); 
    Screen('Flip',window); [secs,KeyCode] = KbWait; 
    
    % Wait for spacebar response to progress. 
    while KeyCode(space) == 0  % While spacebar is not pressed.
        [secs,KeyCode] = KbWait;    % Wait for spacebar response. 
        if KeyCode(escape)         % If escape key pressed:
            sca; ListenChar(); ShowCursor;  % Close all, re-enable keypresses to command window and make mouse cursor visible. 
        end
    end
    
    tic; total_trial_start = GetSecs(); % Record time at start of trials loop. 
    trial_onset = 0;                    % Set individual trial start time to 0.
    
    for k = 1:trials_desired           % For each trial to the total number of trials desired.

        ind_trial_start = GetSecs;       % Record time at start of each trial. 
        
        if exist('timing_diff','var') == 0                                     % If the variable 'timing_diff' does not yet exist (we will only do this once):
            timing_diff = total_trial_start_scan - ind_trial_start;                 % Calculate the difference between the total trial start time and individual trial start time. 
            total_trial_start_scan = total_trial_start_scan - timing_diff;          % Set the total scan start time as the total scan start time - this timing difference. 
        end
        
        trial_start_in_block = ind_trial_start - total_trial_start;         % Set trial start time for block as difference between individual trial start time and total block start time.  
        trial_start_in_scan = ind_trial_start - total_trial_start_scan;     % Set trial start time for scan as difference between individual trial start time and total scan start time. 
        
        if k == 1                       % If this is the first trial:
            trial_start_in_block = 0;       % Set start time of first trial to 0 (for cumulative trial time counter). 
            trial_end_in_block = 0;         % Set counter so that cumulative trial end per block is refreshed at start of each block. 
        end
        
        clear distances_target; % Clear the target RFP each loop, so a new one is generated each time. 
        
        Screen('TextSize', window, 15); % Set smaller text size for fixation cross. 
   
        pp_orientation{i} = orientation_threshold;  pp_contrast{i} = contrast_threshold; ...    % Add participant thresholds to array. 
        pp_ram{i} = ram_threshold;     
        ref=[phi C A];  inp=[orientation_seq{k} contrast_seq{k} ram_seq{k}];                    % Vectors of reference and input conditions.                   
        conditions{i} = lower(block_label(b_condition));                                     % Store attention condition of block. 
        
        % Set default response parameters. 
        response = 'no'; RT = 0; keyLabel = 'N/A'; co_occur = 'no'; num_change = 0; type_changes = 'no_change';
        orientation_tracker = 0; contrast_tracker = 0; ram_tracker = 0;  response_needed= 'N'; key_press_counter = 0;
        
        if task == 'T'
            Screen('DrawTexture', window, scn);                                         % Draw reference RFP texture to screen. 
            DrawFormattedText(window, fixation, 'center', 'center', white);   % Draws condition-specific fixation letter in condition colour. 
            
            t0 = GetSecs; Screen('Flip', window); reference_RFP = GetSecs();    % Flip RFP to screen and record time. 
            trial_offset = reference_RFP - ind_trial_start;                     % Calculate time between start of trial loop and RFP flip. 
            reference_RFP = reference_RFP - (ind_trial_start + trial_offset);
            time_diff = reference_RFP - t0;                                     % Calculate time taken to flip screen. 
            
            waitTime = (stimulus_presentation_time - time_diff);                % Specify stimulus presentation time (- time taken to flip screen).        
            endTime = waitTime + reference_RFP;                                 % Specify end time as time just after flip + stimulus presentation time. 
            
            trial_start_ref_RFP = 0;                                            % Set start time of individual trial to 0.     

            while GetSecs < endTime                                                            % While current time is less than the specified end time:
                    DrawFormattedText(window, fixation, 'center', 'center', white);   % Prepare the temporal between-RFP fixation cross. 
            end
        end
        
        t0 = GetSecs; Screen('Flip', window); trial_fixation = GetSecs;             % Flip fixation cross to screen and record time. 
        time_diff = trial_fixation - t0;                                            % Calculate time taken to flip screen. 
        
        waitTime = (stimulus_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen).
        endTime = waitTime + trial_fixation;                                        % Specify end time as time just after flip + stimulus presentation time. 
        trial_fixation = trial_fixation - (ind_trial_start + trial_offset);         % Record time temporal between-RFP fixation cross was presented. 
        
        while GetSecs < endTime    % While current time is less than the specified end time:
            if exist('distances_target','var') == 0    % If the target RFP texture does not exist: 
                
            changes=inp-ref~=0;         % Logical of nonzero locations (where were input and reference values unequal?). 
            num_change=sum(changes);    % Sum these logicals together. 

            if num_change > 0                                      % If some change from reference values occured:
                fmt=[repmat('%s & ',1,num_change-1) '%s'];              % Create a dynamic format string
                type_changes = sprintf(fmt,block_label{changes});       % Build 'type_changes' variable specifying the change. 
                response_needed = 'Y';                                  % Mark as response needed. 
            end

            if num_change > 1          % If more than one change occured:
                co_occur = 'yes';           % Mark as feature co_occurence. 
            end
            
            o_target = orientation_seq{k}; c_target = contrast_seq{k}; r_target = ram_seq{k};   % Save target values as variables for function below. 
            
            % Call 'calculate_fmri_value' function to calculate feature-value for test RFP and categorise this change (e.g. clockwise, lower contrast...).
            [o_final_target, o_tracker, o_direction, c_final_target, c_tracker, c_direction, r_final_target, r_tracker, r_direction]...
            = calculate_psychophysics_value(o_target, c_target, r_target);
            
            % Add variables to relevant cell arrays.
            orientation_direction{i} = o_direction; orientation_seq{k} = o_final_target;
            contrast_direction{i} = c_direction; contrast_seq{k} = c_final_target;
            ram_direction{i} = r_direction; ram_seq{k} = r_final_target;a;
            target_orientation{i} = orientation_seq{k}; target_contrast{i} = contrast_seq{k}; 
            target_ram{i} = ram_seq{k};  orientation_change_tracker{i} = o_tracker; 
            contrast_change_tracker{i} = c_tracker; ram_change_tracker{i} = r_tracker;
            
            % Create a target RFP with change(s) in the relevant feature dimension and create target RFP texture. 
            distances_target = RFP_generator(orientation_seq{k},contrast_seq{k},ram_seq{k});
            sn = Screen('MakeTexture', window, distances_target);                    
            end
        end
        
       if task == 'T'                                                             % If temporal task:
            Screen('DrawTexture', window, sn);                                          % Draw target RFP to the centre of the screen. 
        elseif task == 'S'                                                         % If spatial task:
            Screen('DrawTexture', window, scn, [], rect_ref_left);                      % Draw reference RFP to left of the screen.  
            Screen('DrawTexture', window, sn, [], rect_target_right);               % Draw target RFP to right of the screen. 
        end
        DrawFormattedText(window, fixation, 'center', 'center', white);   % Draw relevant fixation letter & colour on top of RFPs.  
        
        t0 = GetSecs; Screen('Flip', window); target_RFP = GetSecs;                 % Flip fixation cross and RFPs to screen and record time. 
        time_diff = target_RFP - t0;                                                % Calculate time taken to flip screen. 
        
        waitTime = (stimulus_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen).
        endTime = waitTime + target_RFP;                                            % Specify end time as time just after flip + stimulus presentation time. 

        if task == 'S'                                                             % If spatial task:
            trial_offset = target_RFP - ind_trial_start;                                % Calculate time between start of trial loop and RFP flip. 
        end
   
        target_RFP = target_RFP - (ind_trial_start + trial_offset);                 % Record time target RFP fixation stimulus was presented. 
                                  
        while GetSecs < endTime                                                    % While current time is less than the specified end time:
            DrawFormattedText(window, fixation, 'center', 'center', [1 1 1]);       % Prepare the response period fixation cross. 
        end
        
        t0 = GetSecs; Screen('Flip', window); response_fixation = GetSecs;          % Flip response fixation to screen and record time. 
        time_diff = response_fixation - t0;                                         % Calculate time taken to flip screen. 
        
        waitTime = (response_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen). 
        endTime = waitTime + response_fixation;                                     % Specify end time as time just after flip + stimulus presentation time.
    
        while GetSecs < endTime                                                    % While current time is less than the specified end time:
        [keyisDown, secs, KeyCode] = KbCheck;                                       % Wait for keypress. 
            if key_press_counter == 0
                if keyisDown                                                           % If a key has been pressed. 
                    key_press_counter = 1;
                    keyLabel = KbName(KeyCode);                                             % Store keypress name corresponding to keyboard letter.              
                    RT = secs - response_fixation;                                          % Record response reaction time. 
                    if isequal(keyLabel,'ESCAPE')                                      % If participant pressed the escape key. 
                        sca; ListenChar(); ShowCursor;                                      % Close the screen and re-enable keyboard and mouse input. 
                    end  
                end
            end    
        end
        
        response_fixation = response_fixation - (ind_trial_start + trial_offset);   % Record time response fixation was presented. 
        endTime = total_trial_length + ind_trial_start;           % Specify end time as total trial length + trial start time.  
        
        while GetSecs < endTime                                  % While current time is less than the specified end time:
        
            if isequal(keyLabel,'l') == 0 && isequal(keyLabel,'a') == 0     % If an invalid response was made:
                response = 'no_response'; response_made{i} = response;...       % Add default values to relevant cell arrays.   
                    key_pressed{i} = keyLabel; marked_responses{i} = 'No Response';            
            end
        
            if isequal(keyLabel,'a')         % If 'A' response was made:
                response = 'change';        % Mark as change response. 
            elseif isequal(keyLabel,'l')      % Else if 'L' response was made:
                response = 'no_change';     % Mark as no change response. 
            end
            
            response_made{i} = response; key_pressed{i} = keyLabel;  % Add values to relevant cell arrays. 
            
            % Classify response as hit, miss, correct rejection or false alarm.
            [marked_response] = d_prime_classification_psychophysics(response, type_changes, specific_condition, b_condition);
            marked_responses{i} = marked_response;   % Add value to relevant cell array. 
        end
        
        trial_end = GetSecs();          % Record time at end of each trial.
        trial_end = trial_end - ind_trial_start;                % Calculate total length of trial - trial start time. 
        trial_end_in_block = trial_end_in_block + trial_end;      % Calculate total length of trials so far within a block.
        trial_end_in_scan = trial_end_in_scan + trial_end;        % Calculate total length of trials so far within a scan.

        if task == 'T'                                                         % If temporal task:
            ref_RFP_presentation = trial_fixation - trial_start_ref_RFP;            % Calculate presentation period of reference RFP (period it was visible).  
            reference_RFPs{i} = trial_start_ref_RFP * 1000;                           % Convert values to milliseconds and add to relevant cell arrays. 
            reference_RFP_pres{i} = (trial_fixation - reference_RFP) * 1000;
        end

        % Calculate presentation periods for inter-trial fixation, target RFP and response fixation.
        trial_fixation_presentation = target_RFP - trial_fixation;             
        target_RFP_presentation = response_fixation - target_RFP;
        response_fixation_presentation = trial_end - response_fixation;

        % Convert values to milliseconds and add to relevant cell arrays. 
        task_condition{i} = lower(specific_condition); trial_offsets{i} = trial_offset * 1000; trial_fixations{i} = trial_fixation * 1000; 
        trial_fixation_pres{i} = trial_fixation_presentation * 1000; target_RFPs{i} = target_RFP * 1000;
        target_RFP_pres{i} = target_RFP_presentation * 1000; response_fixations{i} = response_fixation * 1000;
        response_fixation_pres{i} = response_fixation_presentation * 1000; trial_ends{i} = trial_end * 1000;
        total_trial_ends_block{i} = trial_end_in_block * 1000; total_trial_ends_scan{i} = trial_end_in_scan * 1000;
        total_trial_starts_block{i} = trial_start_in_block * 1000; total_trial_starts_scan{i} = trial_start_in_scan * 1000;
        
        % Store trial variables in relevant arrays
        trial_number{i} = i; reference_orientation{i} = phi; reference_contrast{i} = C; reference_ram{i} = A; response_made{i} = response;    
        reaction_time{i} = RT * 1000; change_count{i} = num_change; change_type{i} = lower(type_changes); conjunction{i} = co_occur;
        response_required{i} = response_needed; i = i + 1;
        
        cd(data_storage);   % Change directory to data storage directory (specified above). 
        
        if pid < 10
            pid_directory = (sprintf('participant_0%d', pid)); % Specify unique participant directory. 
        else
            pid_directory = (sprintf('participant_%d', pid)); % Specify unique participant directory.
        end

        if scan_number == 1    % If this is the first scan for the participant:
            if exist(pid_directory,'dir') == 0 % If the unique participant directory does not already exist:
                mkdir(pid_directory); cd(pid_directory);     % Create a directory for the participant and index into this directory. 
            else
                cd(pid_directory);
            end
        elseif scan_number > 1 % If this is not the first scan for the participant:
                cd(pid_directory);    % Index into their unique participant directory. 
        end
        
        if pid < 10
            save(sprintf('%s_s_vs_ns_workspace_pp_0%d_scan_0%d', short_trial_type,pid, scan_number));
        else
            save(sprintf('%s_s_vs_ns_workspace_pp_%d_scan_0%d', short_trial_type,pid, scan_number));
        end
        
%         cd(strcat(root_directory,'/supporting_files'));                                                     % Return to relevant function directory.   
        cd(strcat(root_directory,'\supporting_files'));     % Shuttle.   
    end
    jheapcl;
end
toc    % Calculate elapsed time for all trials.

sca; ListenChar(); ShowCursor; % Close screen and re-enable keypresses and mouse input. 

%% -------------------------- DATA OUTPUT ------------------------------ %%

cd(data_storage);   % Change directory to data storage directory (specified above). 
cd(pid_directory);

% Store all response/feature data in one transposed matrix. 
all_data = [trial_number; task_condition; pp_orientation; reference_orientation; target_orientation; orientation_direction;...
    orientation_change_tracker; pp_contrast; reference_contrast; target_contrast; contrast_direction; contrast_change_tracker;...
    pp_ram; reference_ram; target_ram; ram_direction; ram_change_tracker; key_pressed; response_made; marked_responses; reaction_time;...
    change_count; change_type; conjunction; response_required; task_condition]';            
        
% Convert response/feature data to table with headings.
data_table = cell2table(all_data, 'VariableNames', {'Trial' 'Condition' 'O_Threshold' 'O_Reference' 'O_Target' 'O_Direction' 'O_Tracker'...
    'C_Threshold' 'C_Reference' 'C_Target' 'C_Direction' 'C_Tracker' 'S_Threshold' 'S_Reference' 'S_Target' 'S_Direction' 'S_Tracker' 'Key_Press'...
    'Response' 'Detection' 'RT' 'Change_Count' 'Change_Type' 'Conjunction' 'Response_Required' 'Condition_2'});

if task == 'T'         %If temporal task:
    
    % Store all timing data in one transposed matrix (temporal stores reference RFP and inter-trial fixation times). 
    all_times = [trial_number; task_condition; total_trial_starts_scan; total_trial_starts_block; trial_offsets; reference_RFPs; reference_RFP_pres; trial_fixations; trial_fixation_pres;...
        target_RFPs; target_RFP_pres; response_fixations; response_fixation_pres; trial_ends; total_trial_ends_block; total_trial_ends_scan]';   
    
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames', {'Trial_Number' 'Condition' 'C_Trial_Start_In_Scan' 'C_Trial_Start_In_Block' 'Trial_Offset' 'Trial_Start'...
                                                        'Reference_RFP_Presentation' 'Trial_Fixation' 'Trial_Fixation_Presentation' 'Target_RFP'...
                                                        'Target_RFP_Presentation' 'Response_Fixation' 'Response_Fixation_Presentation' 'Trial_End' 'C_Trial_End_In_Block' 'C_End_In_Scan'});  
                     
elseif task == 'S'     % If spatial task:

    % Store all timing data in one transposed matrix (spatial does not store reference RFP or inter-trial fixation times). 
    all_times = [trial_number; task_condition; total_trial_starts_scan; total_trial_starts_block; trial_offsets; target_RFPs; target_RFP_pres;...
        response_fixations; response_fixation_pres; trial_ends;  total_trial_ends_block; total_trial_ends_scan]';
     
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames',{'Trial_Number' 'Condition' 'C_Trial_Start_In_Scan' 'C_Trial_Start_In_Block' 'Trial_Offset' 'Trial_Start'...
                                                        'RFP_Presentation' 'Response_Fixation' 'Response_Fixation_Presentation' 'Trial_End'...
                                                        'C_Trial_End_In_Block' 'C_End_In_Scan'});
end

if pid < 10
    writetable(data_table, sprintf('s_vs_non_s_%s_%s_pp_0%d_scan_0%d_%s.csv', condition_task, short_trial_type, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));             % Save response data as .csv file.
    writetable(time_table, sprintf('s_vs_non_s_timing_log_%s_%s_pp_0%d_scan_0%d_%s.csv', condition_task, short_trial_type, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));  % Save timing data as .csv.
    save(sprintf('%s_s_vs_ns_workspace_pp_0%d_scan_0%d_%s', short_trial_type,pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));
else
    writetable(data_table, sprintf('s_vs_non_s_%s_%s_pp_%d_scan_0%d_%s.csv', condition_task, short_trial_type, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));             % Save response data as .csv file.
    writetable(time_table, sprintf('s_vs_non_s_timing_log_%s_%s_pp_%d_scan_0%d_%s.csv', condition_task, short_trial_type, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));  % Save timing data as .csv.
    save(sprintf('%s_s_vs_ns_workspace_pp_%d_scan_0%d_%s', short_trial_type,pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));
end
cd(start_dir);
%% --------------------------------------------------------------------- %%