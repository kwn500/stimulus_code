                                %% ---------------- RFP 2AFC FMRI DISCRIMINATION TASK ------------------ %%
                                    % ORIENTATION, RADIAL AMPLITUDE MODULATION & PASSIVE %  
                                    
% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

% Participant must judge whether target RFP is different in a cued feature with respect to a reference RFP. Participant must respond with '1'/'2' or 
% '3'/'4' button presses indicating presence/absence of change. Participants will be cued to attend to one feature per block, or a passive viewing 
% condition. 

%% ---------------- PRE-ALLOCATE CELL ARRAYS FOR SPEED ----------------- %%  
sca; close all; clearvars; KbName('UnifyKeyNames'); %%jheapcl;
start_dir = pwd;    % Store start location as directory to return to at end. 

trial_number = cell(1); condition = cell(1); pp_orientation = cell(1); reference_orientation = cell(1); target_orientation = cell(1);
orientation_direction = cell(1); orientation_change_tracker = cell(1); pp_contrast = cell(1); reference_contrast = cell(1); target_contrast = cell(1); 
contrast_direction = cell(1); contrast_change_tracker = cell(1); pp_ram = cell(1); reference_ram = cell(1); target_ram = cell(1); 
ram_direction = cell(1); ram_change_tracker = cell(1); key_pressed = cell(1); response_made = cell(1); marked_responses = cell(1); 
reaction_time = cell(1); change_count = cell(1); change_type = cell(1); conjunction = cell(1); response_required = cell(1);    
orientation_seq = cell(1); contrast_seq = cell(1); ram_seq = cell(1); task_condition = cell(1); trial_offsets = cell(1); trial_fixations = cell(1);
trial_fixation_pres = cell(1); target_RFPs = cell(1); target_RFP_pres = cell(1); response_fixations = cell(1); response_fixation_pres = cell(1); 
inter_block_intervals = cell(1); inter_block_interval_pres = cell(1); trial_ends = cell(1); cumulative_trial_starts = cell(1); 
cumulative_trial_ends = cell(1); reference_RFPs = cell(1); reference_RFP_pres = cell(1); final_trial_ends = cell(1); cue_times = cell(1); 
cue_presentations = cell(1); total_trial_ends_scan = cell(1); total_trial_starts_scan = cell(1); accurate_timing_ends = cell(1);
accurate_trial_ends = cell(1); accurate_cumulative_ends = cell(1);

orientation_starts = zeros(); contrast_starts = zeros(); shape_starts = zeros(); passive_starts = zeros(); orientation_code = zeros(); 
contrast_code = zeros(); shape_code = zeros(); passive_code = zeros(); orientation_duration = zeros(); contrast_duration = zeros(); shape_duration = zeros(); 
passive_duration = zeros(); inter_block_starts = zeros(); inter_block_duration = zeros(); inter_block_ends_time = zeros(); inter_block_code = zeros(); 
cue_starts = zeros(); cue_duration = zeros(); cue_code = zeros(); trial_end_counting_scan = zeros(); trial_end_counting = zeros(); 
no_change_starts = zeros(); no_change_code = zeros(); no_change_duration = zeros(); shape_change_starts = zeros(); shape_change_code = zeros(); 
shape_change_duration = zeros(); contrast_change_starts = zeros(); contrast_change_code = zeros(); contrast_change_duration = zeros();
orientation_change_starts = zeros(); orientation_change_code = zeros(); orientation_change_duration = zeros(); hit_starts = zeros(); hit_code = zeros(); 
hit_duration = zeros(); miss_starts = zeros(); miss_code = zeros(); miss_duration = zeros(); FA_starts = zeros(); FA_code = zeros(); FA_duration = zeros(); 
CR_starts = zeros(); CR_code = zeros(); CR_duration = zeros(); no_response_starts = zeros(); no_response_code = zeros(); no_response_duration = zeros();

%% ------------------------- GENERAL SETUP ----------------------------- %%

% CHECK THESE PIXELS PER DEGREE MATCH NUMBER OF PIXELS PER DEGREE IN RFP_GENERATOR FUNCTION. 
display_ppd = 46.4; % Correct for propixx projector (in mri scanner). 

% Create and present dialog box that gathers user input.
prompt = {'Enter Participant Number:','Enter Scan Number:'};
dlg_title = 'Participant Information';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

pid = str2double(answer{1}); scan_number = str2double(answer{2});   % Store participant and scan numbers as variables from GUI input. 

% Vital experimental information - check this is correct before running!. 
trials_desired = 10;                 % Number of trials.
task = 'T';                         % Type of task - 'T' = temporal, 'S' = spatial. 
threshold_repeat = 3;               % Number of threshold runs (ensures you load in correct thresholds). 
condition_task = 'temporal';        % Task condition string (for data output)- 'temporal'/'spatial'. 

block_repeats_desired = 4;          % Number of block repeats (one scan). 
baseline_length = 7.5;              % Length of baseline period between blocks (seconds). 
prob_change = 0.8;                  % Probability of change in stimulus feature (0.8 = 20% likelihood change). 
stimulus_presentation_time = 0.2;   % Stimulus presentation time (s). 
response_presentation_time = 0.8;   % Response period duration (s). 
total_trial_length = 1.5;           % Total trial duration (s). 
dummy_volume_length = 9;            % Length of dummy period (TRs to 'throw away at start of scan) (s).
dummy_volume_counter = 1;           % Establish counter variable (ensures only present dummy volumes once per scan). 
trial_cue_length = 1.5;             % Total length of single trial (s).
dummy_trials_desired = dummy_volume_length /total_trial_length;     % Calculate number of trials required to fill dummy period length. 
between_block_interval = baseline_length + trial_cue_length;        % Total length of between-block period (s).

%Calculate total trial duration if including an inter-block interval (s). 
total_inter_block_interval_trial_length = total_trial_length + baseline_length; 
total_trial_cue_length = total_trial_length + between_block_interval;

% Sets up directories for data storage and access to functions. 
%  root_directory = ('/Users/Kirstie/Documents/rfp_code_data/code');    % Top-level directory- mac.
root_directory = ('\Users\Ryan Maloney\Documents\Kirstie\orientation_contrast_shape\');    % Shuttle. 

if pid < 10; % If the particpant number is not two-digits (i.e. less than 10):
    matlab_data = sprintf('fmri_temporal_pp_0%d_scan_0%d_data', pid, scan_number);   % Specify file and directory names with a leading 0 for pp number. 
    pid_directory = (sprintf('participant_0%d', pid)); 
    
    % Shuttle:
    file_path = strcat(root_directory,sprintf('\\data_output\\experimental_thresholds\\0%d_experimental_thresholds',pid)); cd(file_path);  
       
    % Mac:
    file_path = strcat(root_directory,sprintf('/data_output/experimental_thresholds/0%d_experimental_thresholds',pid)); cd(file_path); 
else         % Else
    matlab_data = sprintf('fmri_temporal_pp_%d_scan_0%d_data', pid, scan_number);    % Specify file and directory names without a leading 0 for pp number. 
    pid_directory = (sprintf('participant_%d', pid));
    
    % Shuttle:
    file_path = strcat(root_directory,sprintf('\\data_output\\experimental_thresholds\\%d_experimental_thresholds',pid)); cd(file_path); 
    
    % Mac:
    file_path = strcat(root_directory,sprintf('/data_output/experimental_thresholds/%d_experimental_thresholds',pid)); cd(file_path); 
end

% data_storage = strcat(root_directory,(sprintf('/data_output/fmri_%s_data',condition_task)));    % Set data storage directory - mac. 

                                             % Shuttle.          
data_storage = strcat(root_directory,(sprintf('\\data_output\\fmri_%s_data\\',condition_task)));    % Shuttle.            

full_pid_directory = strcat(data_storage,'\',pid_directory);    % Specify full path to unique participant folder - Shuttle.
% full_pid_directory = strcat(data_storage,'/',pid_directory);    % Specify full path to unique participant folder - Mac.

scan_dir = sprintf('scan_0%d',scan_number);
full_scan_directory = strcat(data_storage,'\',pid_directory,'\',scan_dir);  % Specify full path to unique scan directory - Shuttle. 
%full_scan_directory = strcat(data_storage,'/',pid_directory,'/',scan_dir);  % Specify full path to unique scan directory - Mac. 

orientation_threshold = 0; contrast_threshold = 0; ram_threshold = 0;   % Establish counter variables for average participant threshold values. 

for iii  = 1:threshold_repeat    % for one to the number of threshold repeats each participant completed:
    
    if pid < 10;
        load(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat', condition_task, iii, pid));  % Load in participants threshold values.
    else
        load(sprintf('thresholds_%s_scan_0%d_pp_%d.mat', condition_task, iii, pid));  % Load in participants threshold values.
    end
    
    % Add all thresholds for each feature to single variable (running total). 
    orientation_threshold = orientation_threshold + orientationT;   contrast_threshold = contrast_threshold + contrastT; ram_threshold = ram_threshold + RAMT;   
end

% Calculate the average of these thresholds and double these values. 
orientation_threshold = orientation_threshold / threshold_repeat; orientation_threshold = orientation_threshold * 2;
contrast_threshold = contrast_threshold / threshold_repeat; contrast_threshold = contrast_threshold * 2;
ram_threshold = ram_threshold / threshold_repeat; ram_threshold = ram_threshold * 2;

% cd(strcat(root_directory,'/supporting_files'));                                                 % Return to relevant function directory - mac.
cd(strcat(root_directory,'\supporting_files')); 

% General experiment 'housekeeping'. 
block_label = {'ORIENTATION' 'CONTRAST' 'SHAPE' 'PASSIVE'};         % Specify block labels. 
blocks_desired = [1 2 3 4];                                         % Specify block numbers.
fixation_label = {'O' 'C' 'S' 'P'};                                 % Specify fixation letters. 
fixation_colours = {[1 0 0], [0 1 0], [0 0 1], [1 1 0]};            % Specify fixation colours. 

block_order = [];   % Establish empty matrix for randomisation of block conditions. 

for t = 1:block_repeats_desired;   % For 1 to the number of block repeats desired:
    block = blocks_desired(randperm(length(blocks_desired)));   % Randomise the block order within each block (1 repeat of each condition). 
    block_order = [block_order, block];                         % Add these four randomised conditions to an overall block order matrix. 
end

space = KbName('space'); escape = KbName('escape');                 % Set keynames for potential response keypresses.
trigger = KbName('5%');  change1 = '1!'; change2 = '2@'; no_change1 = '3#'; no_change2 = '4$'; 

PsychDefaultSetup(2); screens = Screen('Screens');  screenNumber = max(screens);    % Call default PTB settings & draw to max screen.
black = BlackIndex(screenNumber); white = WhiteIndex(screenNumber);                 % Define black (0) & white (1).
  
%% -------------------------- SCREEN SETUP ----------------------------- %%
% PsychDebugWindowConfiguration;      % Screen transparency toggle for debug mode. 

[window, windowRect] = PsychImaging('OpenWindow', screenNumber, (128/255));     % Open gray screen window. 
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');      % Alpha-blending for anti-aliased lines.

% deal with the mouse cursor
SetMouse(5000,0); % Move mouse to top right corner
HideCursor(window); % Make sure the mouse is hidden

ListenChar(0);       % Hide mouse cursor & disable keypress output to the command window. 
% PsychHID('KbQueueCreate');                                                      % Establish keyboard queue variable (needed to flush events later). 
%-------------------------------------------------------------------------
% set allowed keys

% allowing any button 1:4 so can use whichever is comfiest
% also manually allowing keys:
%   11:14 - 1,2,3,4 on linux system
%   1!, 2@, 3# and 4$ (49-52 on a win system)
allowedKeys = zeros(1,256);
targetKeys = [KbName({'1','2','3','4'}) 11:14 49:52];
quitKey = KbName('q');
allowedKeys([targetKeys quitKey]) = 1;

% create keyboard queue (but don't start it yet!)
KbQueueCreate([],allowedKeys);

%-------------------------------------------------------------------------


[screenXpixels, screenYpixels] = Screen('WindowSize', window);                  % Window size (pixels). 
[xCenter, yCenter] = RectCenter(windowRect);                                    % Centre of screen (pixels). 

% No longer controls stimulus size (this is hard-coded in Rich's RFP_generator function), but controls position of spatial targets. 
baseRect = [0 0 500 500]; 

xCenter_left = xCenter - (5 * display_ppd) ; xCenter_right = xCenter + (5 * display_ppd) ; % Position RFPs 5° visual angle from centre.

% Set positions of target % reference RFP (reference - left, target - right).  
rect_ref_left = CenterRectOnPointd(baseRect, xCenter_left, yCenter);
rect_target_right = CenterRectOnPointd(baseRect, xCenter_right, yCenter); 

%% ----------------------- INTRODUCTION SCREEN ------------------------- %%

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial');    % Set text parameters. 

% Present relevant instructions for selected condition. 
line1 = 'Waiting for trigger\n\n'; line2 = '1 or 2 = change\n';
line3 = '3 or 4 = no change';

DrawFormattedText(window, [line1 line2 line3], 'center', 'center', white);
Screen('Flip',window); [secs,KeyCode] = KbWait; % Draw text to screen and wait for keypress.    

% % % % % Wait for trigger response to progress. 
% % % % while KeyCode(trigger) == 0;  % While trigger is not pressed.
% % % %     [secs,KeyCode] = KbWait;    % Wait for trigger response. 
% % % %     if KeyCode(escape);         % If trigger pressed:
% % % %         sca; ListenChar(); ShowCursor;  % Close all, re-enable keypresses to command window and make mouse cursor visible. 
% % % %         RestrictKeysForKbCheck([]); % Re-allow keypresses from any key, rather than response keypresses. 
% % % %     end
% % % % end

keyPress = 'not_five';
while keyPress ~= '5'
    [~,key] = KbWait;
    keyPress = KbName(key);
    if iscell(keyPress)
        keyPress = keyPress{1};
    end
end

% start recording keypresses
KbQueueStart;


%% ----------------- RADIAL FREQUENCY PATTERN SETUP -------------------- %%

% Set parameters for RFP. 
C = 0.5;	% Contrast.  
A = 0.5;    % Radial modulation amplitude (A < 1).  
phi = 0;    % Phi - angular rotation. 
        
distances = RFP_generator(phi,C,A); % Create reference RFP. 
scn = Screen('MakeTexture', window, distances); % Convert reference RFP number matrix to texture.
   
%% ------------------------- EXPERIMENT LOOPS -------------------------- %%

topPriorityLevel = MaxPriority(window); Priority(topPriorityLevel);                 % Set script to maximum processing priority. 
% % % RestrictKeysForKbCheck([KbName('1!'),KbName('2@'),KbName('3#'),KbName('4$')]);      % Code will only look for 1,2,3 & 4 keypresses (matches 4-button box). 
% % % PsychHID('KbQueueStart');                                                           % Establish keyboard queue start (needed to flush keyboard events later).

Screen('TextSize', window, 15); % Set smaller text size for fixation cross. 

%% ---------------------------- BLOCK LOOP ----------------------------- %%

% Establish counter variables. 
i = 1; y = 1; oo = 1;
trial_end_in_scan = 0; total_trial_start_scan = 0; cumulative_trial_end = 0;
no_change_counter = 1;  orientation_change_counter = 1; contrast_change_counter = 1; shape_change_counter = 1;
orientation_starts_counter = 1; contrast_counter = 1; shape_counter = 1; passive_counter = 1; 
orientation_counter = 1; contrast_counter = 1; shape_counter = 1; passive_counter = 1;
hit_counter = 1; miss_counter = 1; CR_counter = 1; FA_counter = 1; no_response_counter = 1; between_block_counter = 1; cue_time = 0;

for block_number = 1:(length(block_order));    % For each of the blocks within the randomised 'block_order' array:
    
%     block_order = [];   % Establish empty matrix for randomisation of block conditions. 

%     for t = 1:block_repeats_desired;   % For 1 to the number of block repeats desired:
%         block = blocks_desired(randperm(length(blocks_desired)));   % Randomise the block order within each block (1 repeat of each condition). 
%         block_order = [block_order, block];                         % Add these four randomised conditions to an overall block order matrix. 
%     end
    
    block_condition = block_order(block_number);   % Select block number from randomised array.
    
    for a = 1:trials_desired;           % For number of trials desired:

    o_randomiser = rand(); c_randomiser = rand(); r_randomiser = rand();    % Establish random numbers for loop iterations. 
    orientation_seq{a} = phi; contrast_seq{a} = C; ram_seq{a} = A;          % Fill cell array position with default reference value. 
    
        if o_randomiser > prob_change;                      % If random value is above 0.8 (20% chance).
            orientation_seq{a} = orientation_threshold;         % Add orientation value to orientation trial sequence.  
        end                                                     % Repeat process for contrast & shape. 
        if c_randomiser > prob_change;   
            contrast_seq{a} = contrast_threshold;
        end
        if r_randomiser > prob_change;
            ram_seq{a} = ram_threshold;
        end
    end 
    
    % Shuffle trial order within these sequences.
    orientation_seq = orientation_seq(randperm(length(orientation_seq))); 
    contrast_seq = contrast_seq(randperm(length(contrast_seq)));  
    ram_seq = ram_seq(randperm(length(ram_seq)));  
    
    % Select fixation letter & trial cue matching block. 
    fixation = fixation_label{block_condition};  
    trial_cue = block_label{block_order(block_number)};
    
    if block_number < length(block_order);                     % If j less than the number of specified blocks (in total). 
       trial_cue_2 = block_label{block_order(block_number+1)};       % Set the trial cue to the following block (i.e. set the trial cue to the block coming next).
    end
    
    trial_end_in_block = 0;             % Establish counter variable for cumulative trial ends (refreshed every block).
    
%% ----------- DUMMY VOLUMES & INTER-BLOCK INTERVAL + TRIAL CUE -------- %%

    if dummy_volume_counter == 1;       % If the dummy counter is equal to 1 (if dummy volumes have not yet been presented):  
        
        ind_trial_start_dummy = GetSecs();  % Get time at start of dummy volume presentation. 
        
        for q = 1:dummy_trials_desired;     % For 1 to number of dummy trials specied:
            
            Screen('DrawTexture', window, scn);                                 % Draw reference RFP texture to screen. 
            DrawFormattedText(window, 'P', 'center', 'center', white);          % Draws condition-specific fixation letter. 
            
            t0 = GetSecs; Screen('Flip', window); reference_RFP = GetSecs();    % Flip RFP to screen and record time. 
            time_diff = reference_RFP - t0;                                     % Calculate time taken to flip screen. 

            waitTime = (stimulus_presentation_time - time_diff);                % Specify stimulus presentation time (- time taken to flip screen).        
            endTime = waitTime + reference_RFP;                                 % Specify end time as time just after flip + stimulus presentation time.  

            while GetSecs < endTime;                                            % While current time is less than the specified end time:
                DrawFormattedText(window, 'P', 'center', 'center', white);      % Prepare the temporal between-RFP fixation cross. 
            end

            t0 = GetSecs; Screen('Flip', window); trial_fixation = GetSecs;             % Flip fixation cross to screen and record time. 
            time_diff = trial_fixation - t0;                                            % Calculate time taken to flip screen. 

            waitTime = (stimulus_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen).
            endTime = waitTime + trial_fixation;                                        % Specify end time as time just after flip + stimulus presentation time. 

            while GetSecs < endTime;                        % While current time is less than the specified end time:
                if exist('distances_target','var') == 0;        % If the target RFP texture does not exist: 

                o_target = orientation_seq{q}; c_target = contrast_seq{q}; r_target = ram_seq{q};   % Save target values as variables for function below. 

                % Call 'calculate_fmri_value' function to calculate feature-value for test RFP and categorise this change (e.g. clockwise, lower contrast...).
                [o_final_target, o_tracker, o_direction, c_final_target, c_tracker, c_direction, r_final_target, r_tracker, r_direction]...
                = calculate_fmri_value(o_target, c_target, r_target);

                % Store final target values in relevant sequence arrays. 
                orientation_seq{q} = o_final_target; contrast_seq{q} = c_final_target; ram_seq{q} = r_final_target; 
                
                % Create a target RFP with change(s) in the relevant feature dimension and create target RFP texture. 
                distances_target = RFP_generator(orientation_seq{q},contrast_seq{q},ram_seq{q});
                sn = Screen('MakeTexture', window, distances_target);                    
                end
            end

            Screen('DrawTexture', window, sn);                              % Draw target RFP to the centre of the screen
            DrawFormattedText(window, 'P', 'center', 'center', white);      % Draws condition-specific fixation letter. 

            t0 = GetSecs; Screen('Flip', window); target_RFP = GetSecs;     % Flip fixation cross and RFP to screen and record time. 
            time_diff = target_RFP - t0;                                    % Calculate time taken to flip screen. 

            waitTime = (stimulus_presentation_time - time_diff);            % Specify stimulus presentation time (- time taken to flip screen).
            endTime = waitTime + target_RFP;                                % Specify end time as time just after flip + stimulus presentation time. 

            while GetSecs < endTime;                                        % While current time is less than the specified end time:
                DrawFormattedText(window, 'P', 'center', 'center', white);      % Prepare the response period fixation cross. 
            end

            t0 = GetSecs; Screen('Flip', window); response_fixation = GetSecs;  % Flip response fixation to screen and record time. 
            time_diff = response_fixation - t0;                                 % Calculate time taken to flip screen. 

            waitTime = (response_presentation_time - time_diff);                % Specify response presentation time (- time taken to flip screen). 
            endTime = waitTime + response_fixation;                             % Specify end time as time just after flip + response presentation time.

            while GetSecs < endTime;                                            % While current time is less than the specified end time:
                DrawFormattedText(window, 'X', 'center', 'center', black);          % Draw black fixation cross. 
            end
            
            if q == dummy_trials_desired;       % If counter variable is equal to total number of dummy trials desired:
       
                % Specify end time as length of inter_block interval (fixation) + length of dummy volume trials + start of dummy trials:
                endTime = baseline_length + dummy_volume_length + ind_trial_start_dummy;    

                Screen('Flip', window);  % Flip inter-block fixation to screen and record time. 

                while GetSecs < endTime;                                % While current time is less than the specified end time:
                    instructions = sprintf('%s',trial_cue);       % Prepare cue with relevant upcoming condition name. 
                    DrawFormattedText(window, sprintf(instructions), 'center', 'center', white);    % Draw this text to the back screen buffer.
                end

                % Specify end time as length of cue period + inter_block_interval + length of dummy volume trials + start of dummy trials:
                endTime = trial_cue_length + baseline_length + dummy_volume_length + ind_trial_start_dummy;

                Screen('Flip', window);     % Flip trial_cue to screen and record time. 

                while GetSecs < endTime;    % Do nothing until current time is equal to the specified end time:
                end
                
                cue_final_end = GetSecs();  % Record time at end of of cue presentation period. 
            end
        end  
        
        % Set end of cue time as end of cue presentation time - the time at the start of the dummy trials + the inter-block (fixation) length. 
            % i.e: works out how long the trial cue was presented for. 
        cue_presentation = cue_final_end - (ind_trial_start_dummy + baseline_length);

        dummy_volume_counter = 2;   % Increment dummy volume counter to 2 (means there will be no further repetitions of the dummy volumes). 
    end
    
    total_trial_start = GetSecs();  % Record time at start of trials loop (time is refreshed with each block). 
    
%% ----------------- TRIAL LOOP (REPEATED EACH BLOCK) ------------------ %% 

    for k = 1:trials_desired;            % For each trial to the total number of trials desired. 

        ind_trial_start = GetSecs;       % Record time at start of each trial. 
        tt = 1; pp = 1; ppp = 1;         % Establish counter variables. 

        if k == 1;                       % If this is the first trial:
            cumulative_trial_start = 0;     % Set start time of first trial to 0 (for cumulative trial time counter).   
        else 
            % Set cumulative time of trial to start time of individual trial - start time of block of trials. 
            cumulative_trial_start = ind_trial_start - total_trial_start; 
        end

        clear distances_target; % Clear the target RFP each loop, so a new one is generated each time. 

        % Add participant thresholds to array. 
        pp_orientation{i} = orientation_threshold;  pp_contrast{i} = contrast_threshold; pp_ram{i} = ram_threshold;     
        ref=[phi C A];  inp=[orientation_seq{k} contrast_seq{k} ram_seq{k}];        % Vectors of reference and input conditions.                   
        condition{i} = lower(block_label(block_condition));                         % Store attention condition of block. 
        trial_condition = char(condition{i});                                       % Save attention condition as variable for later function. 
        
        % Store start time of trial cue and presentation period in relevant cell arrays.
        cue_times{i} = cue_time * 1000; cue_presentations{i} = cue_presentation * 1000;      
        
        % Set default response parameters. 
        response = 'no'; RT = 'no_response'; keyLabel = 'N/A'; co_occur = 'no'; num_change = 0; type_changes = 'no_change';
        orientation_tracker = 0; contrast_tracker = 0; ram_tracker = 0;  response_needed= 'N'; inter_block_interval_presentation = 0;
        
        if task == 'T';                                                         % If the temporal task has been specified:
            Screen('DrawTexture', window, scn);                                     % Draw reference RFP texture to screen. 
            DrawFormattedText(window, fixation, 'center', 'center', white);         % Draw condition-specific fixation letter. 

            t0 = GetSecs; Screen('Flip', window); reference_RFP = GetSecs();    % Flip RFP to screen and record time. 
            trial_offset = reference_RFP - ind_trial_start;                     % Calculate time between start of trial loop and RFP flip. 
            time_diff = reference_RFP - t0;                                     % Calculate time taken to flip screen. 

            waitTime = (stimulus_presentation_time - time_diff);                % Specify stimulus presentation time (- time taken to flip screen).        
            endTime = waitTime + reference_RFP;                                 % Specify end time as time just after flip + stimulus presentation time. 

            trial_start_ref_RFP = 0;                                            % Set start time of individual trial to 0.     

            while GetSecs < endTime;                                                  % While current time is less than the specified end time:
                    DrawFormattedText(window, fixation, 'center', 'center', white);     % Prepare the temporal between-RFP fixation letter. 
            end
        end

        t0 = GetSecs; Screen('Flip', window); trial_fixation = GetSecs;             % Flip fixation letter to screen and record time. 
        time_diff = trial_fixation - t0;                                            % Calculate time taken to flip screen. 

        waitTime = (stimulus_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen).
        endTime = waitTime + trial_fixation;                                        % Specify end time as time just after flip + stimulus presentation time. 
        trial_fixation = trial_fixation - (ind_trial_start + trial_offset);         % Record time temporal between-RFP fixation cross was presented. 

        while GetSecs < endTime;                        % While current time is less than the specified end time:
            if exist('distances_target','var') == 0;        % If the target RFP texture does not exist: 

            changes=inp-ref~=0;                             % Logical of nonzero locations (where were input and reference values unequal). 
            num_change=sum(changes);                        % Sum these logicals together. 

            if num_change > 0;                                      % If some change from reference values occured:
                fmt=[repmat('%s & ',1,num_change-1) '%s'];              % Create a dynamic format string
                type_changes = sprintf(fmt,block_label{changes});       % Build 'type_changes' variable specifying the change. 
                response_needed = 'Y';                                  % Mark as response needed. 
            end;

            if num_change > 1;                                      % If more than one change occured:
                co_occur = 'yes';                                       % Mark as feature co_occurence. 
            end

            o_target = orientation_seq{k}; c_target = contrast_seq{k}; r_target = ram_seq{k};   % Save target values as variables for function below. 

            % Call 'calculate_fmri_value' function to calculate feature-value for test RFP and categorise this change (e.g. clockwise, lower contrast...).
            [o_final_target, o_tracker, o_direction, c_final_target, c_tracker, c_direction, r_final_target, r_tracker, r_direction]...
            = calculate_fmri_value(o_target, c_target, r_target);

            % Add variables to relevant cell arrays.
            orientation_direction{i} = o_direction; orientation_seq{k} = o_final_target; contrast_direction{i} = c_direction; ram_change_tracker{i} = r_tracker; 
            contrast_seq{k} = c_final_target; ram_direction{i} = r_direction; ram_seq{k} = r_final_target; target_orientation{i} = orientation_seq{k};
            target_contrast{i} = contrast_seq{k}; target_ram{i} = ram_seq{k};  orientation_change_tracker{i} = o_tracker; contrast_change_tracker{i} = c_tracker;
            
            % Create a target RFP with change(s) in the relevant feature dimension and create target RFP texture. 
            distances_target = RFP_generator(orientation_seq{k},contrast_seq{k},ram_seq{k});
            sn = Screen('MakeTexture', window, distances_target);                    
            end
        end

        if task == 'T';                                                             % If temporal task:
            Screen('DrawTexture', window, sn);                                          % Draw target RFP to the centre of the screen. 
        elseif task == 'S';                                                         % If spatial task:
            Screen('DrawTexture', window, scn, [], rect_ref_left);                      % Draw reference RFP to left of the screen.  
            Screen('DrawTexture', window, sn, [], rect_target_right);                   % Draw target RFP to right of the screen. 
        end
        DrawFormattedText(window, fixation, 'center', 'center', white);             % Draw relevant fixation letter & colour on top of RFPs.  

        t0 = GetSecs; Screen('Flip', window); target_RFP = GetSecs;                 % Flip fixation cross and RFPs to screen and record time. 
        time_diff = target_RFP - t0;                                                % Calculate time taken to flip screen. 

        waitTime = (stimulus_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen).
        endTime = waitTime + target_RFP;                                            % Specify end time as time just after flip + stimulus presentation time. 

        if task == 'S';                                                             % If spatial task:
            trial_offset = target_RFP - ind_trial_start;                                % Calculate time between start of trial loop and RFP flip. 
        end

        target_RFP = target_RFP - (ind_trial_start + trial_offset);                 % Record time target RFP fixation stimulus was presented. 

        while GetSecs < endTime;                                                    % While current time is less than the specified end time:
            DrawFormattedText(window, fixation, 'center', 'center', white);         % Prepare the response period fixation letter. 
        end

        t0 = GetSecs; Screen('Flip', window); response_fixation = GetSecs;          % Flip response fixation to screen and record time. 
        time_diff = response_fixation - t0;                                         % Calculate time taken to flip screen. 

        waitTime = (response_presentation_time - time_diff);                        % Specify stimulus presentation time (- time taken to flip screen). 
        endTime = waitTime + response_fixation;                                     % Specify end time as time just after flip + stimulus presentation time.

        trialStartTime = GetSecs;
        
        while GetSecs < endTime;   %&&&                                                 % While current time is less than the specified end time:
% % %         [keyisDown, secs, KeyCode] = KbCheck;                                       % Wait for keypress.  
% % %             if keyisDown;                                                           % If a key has been pressed. 
% % %                 keyLabel = KbName(KeyCode);                                             % Store keypress name corresponding to keyboard letter.              
% % %                 RT = secs - response_fixation;                                          % Record response reaction time. 
% % %                 if isequal(keyLabel,'ESCAPE');                                      % If participant pressed the escape key. 
% % %                     sca; ListenChar(); ShowCursor;                                      % Close the screen and re-enable keyboard and mouse input.                                                              
% % %                 end  
% % %             end    


            % check to see if keyboard been pressed
            [pressed,firstPress,~,~,~] = KbQueueCheck;
            if pressed
                % first check if they've pressed any of the response keys
                if any(firstPress(targetKeys))
                    keyLabel = KbName(firstPress);
                                    
                    %% TODO - make this more elegant using firstPress
                    RT = GetSecs - trialStartTime; 
                    %%%
                    if iscell(keyLabel)
                        keyLabel = keyLabel{1};
                    end
                end
                % now check if they've hit quit
                if firstPress(quitKey)
                    finished = true;
                    currMessage = 'Finished';
                end

                KbQueueFlush; %clear out so we don't fill up the buffer

            end

        end %&&&

        
        response_fixation = response_fixation - (ind_trial_start + trial_offset);   % Record time response fixation was presented. 
        
        endTime = total_trial_length + ind_trial_start;             % Specify end time as total trial length + trial start time.  

% % %         FlushEvents(['keyisDown']);                                  % Flush events stored in keyboard buffer. 

        while GetSecs < endTime;                                    % While current time is less than the specified end time:

            if isequal(keyLabel,'N/A');                                         % If no response was made:
                key_pressed{i} = keyLabel; marked_responses{i} = 'Miss';            % Add default values to relevant cell arrays.       
            end

            if isequal(keyLabel, change1) == 1 || isequal(keyLabel, change2) == 1;              % If they pressed '1'/'2' (indicated change). 
                response = 'change'; key_pressed{i} = keyLabel;                                     % Add response values to relevant cell arrays. 
            elseif isequal(keyLabel, no_change1) == 1 || isequal(keyLabel, no_change2) == 1;    % If they pressed '3'/'4' (indicated no change).
                response = 'no_change'; key_pressed{i} = keyLabel;                                  % Add response values to relevant cell arrays.     
            end

            % Classify response as hit, miss, correct rejection or false alarm.
            [marked_response] = d_prime_classification(response, type_changes, block_condition);     
            marked_responses{i} = marked_response;  % Add value to relevant cell array.
             
            % Set time of trial end (not genuine time), but specified length of trial. 
                % This is needed to write out log files within loops successfully (which is in turn needed for accurate timing). 
                % Check trial end values in the excel timing log to ensure that this specified time is accurate.     
            trial_end = total_trial_length; 
            
            if k ~= trials_desired;     % If trial number is not equal to the total number of trials desired per block:
                
                if pp  <  2;            % If the counter variable is less than 2 (ensures loop only completed once). 
                    
                    inter_block_interval = 0;       % Set inter-block interval to default value. 
                    final_trial_end = trial_end;    % Set final trial end to the default trial end specified time (from above). 
            
                    % Again, these do not reflect accurate timings, just specified times for each trial/ block (needed for log files). Excel timing log contains accurate timing. 
                        % If the excel file has correct timing, then log files are accurate representations of trial/block timing. 
                    trial_end_in_block = trial_end_in_block + final_trial_end;      % Calculate total duration of trials so far within a block (refreshed each block).
                    trial_end_in_scan = trial_end_in_scan + final_trial_end;        % Calculate total duration of trials so far within a scan.

                    % If this is the first trial of the scan or if the number of trials is equal to the total number of trials specified within each block.
                    if i == 1;                                                        
                        total_trial_start_scan = total_trial_start_scan + baseline_length + trial_cue_length;       % Set variable to time taken for the first inter-block fixation and trial cue. 
                    else                                                                                        % Otherwise, if this is not the first trial of the scan:
                        total_trial_start_scan = total_trial_start_scan + trial_end;                                % Add the length of the trial to the cumulative counter variable.
                    end
               
                    % Add these values to the relevant cell arrays. 
                    total_trial_starts_scan{i} = total_trial_start_scan; total_trial_ends_scan{i} = trial_end_in_scan * 1000;

                    response_fixation_presentation = trial_end - response_fixation;         % Calculate presentation period of reponse fixation from trial end. 
                    ref_RFP_presentation = trial_fixation - trial_start_ref_RFP;            % Calculate presentation period of reference RFP (period it was visible).  

                    % Convert reference RFP start and presentation values to ms and add to relevant cell arrays. 
                    reference_RFPs{i} = trial_start_ref_RFP * 1000; 
                    reference_RFP_pres{i} = ref_RFP_presentation * 1000;

                    % Calculate presentation periods for fixation between reference and target RFPs and target RFP.
                    trial_fixation_presentation = target_RFP - trial_fixation; target_RFP_presentation = response_fixation - target_RFP;

%% ------- CREATION OF FSL AND FREESURFER/MR VISTA TIMING ARRAYS ------- %%

                    if i == 1;  % If this is the first trial of the scan:
                        
                        % Set the start of the first inter-block fixation to 0s & specify code.
                        inter_block_starts(between_block_counter) = 0;  inter_block_code(between_block_counter) = 0; 
                        % Set duration of inter-block fixation to time specified by baseline length variable (s).
                        inter_block_duration(between_block_counter) = baseline_length;  

                        % Set end time of inter-block fixation to start time + length of inter-block fixation.
                        inter_block_ends_time(between_block_counter) = inter_block_starts(between_block_counter) + baseline_length;

                        % Set start time of first trial cue as end of inter-block fixation.
                        cue_starts(between_block_counter) = baseline_length; cue_code(between_block_counter) = 1; 
                        % Set duration of trial cue to time specified by trial cue length variable (s).
                        cue_duration(between_block_counter) = trial_cue_length; 
                    end

                    if isequal(type_changes,'no_change');                                   % If no change between reference and target RFP values occured:

                        % Specify start of trial as cumulative trial counter time (not updated until lower in the loop). 
                        no_change_starts(no_change_counter) = total_trial_starts_scan{i};       

                        no_change_duration(no_change_counter) = trial_end;                      % Specify no change duration as length of singular trial. 

                         % Specify no change trial code and increment counter variable.   
                        no_change_code(no_change_counter) = 5; no_change_counter = no_change_counter + 1;  
                    end

                                                % Repeat this process for trials including orientation, contrast and shape changes.  
                    if strfind(type_changes,'ORIENTATION') > 0;
                        orientation_change_starts(orientation_change_counter) = total_trial_starts_scan{i};
                        orientation_change_duration(orientation_change_counter) = trial_end;
                        orientation_change_code(orientation_change_counter) = 2; orientation_change_counter = orientation_change_counter + 1;
                     end

                    if strfind(type_changes,'CONTRAST') > 0;
                        contrast_change_starts(contrast_change_counter) = total_trial_starts_scan{i};
                        contrast_change_duration(contrast_change_counter) = trial_end;
                        contrast_change_code(contrast_change_counter) = 3; contrast_change_counter = contrast_change_counter + 1;
                    end
         
                    if strfind(type_changes,'SHAPE') > 0;
                        shape_change_starts(shape_change_counter) = total_trial_starts_scan{i};
                        shape_change_duration(shape_change_counter) = trial_end;
                        shape_change_code(shape_change_counter) = 4; shape_change_counter = shape_change_counter + 1;
                    end

                    if isequal(marked_response,'Hit');  % If participants response was marked as a hit:

                        % Specify the time at the start of this trial as cumulative trial counter time (not updated until lower in the loop).
                        hit_starts(hit_counter) = (total_trial_starts_scan{i}); 

                         % Specify duration as length of singular trial, set response code & increment counter variable. 
                        hit_duration(hit_counter) = trial_end; hit_code(hit_counter) = 2; hit_counter = hit_counter + 1;
                    end

                                                      % Repeat this process for trials with miss, correct rejection, false alarm and no response data. 
                    if isequal(marked_response,'Miss');
                        miss_starts(miss_counter) = (total_trial_starts_scan{i});
                        miss_duration(miss_counter) = trial_end; miss_code(miss_counter) = 3; miss_counter = miss_counter + 1;
                    end
                    
                    if isequal(marked_response,'Correct Rejection');
                        CR_starts(CR_counter) = (total_trial_starts_scan{i});
                        CR_duration(CR_counter) = trial_end; CR_code(CR_counter) = 4; CR_counter = CR_counter + 1;
                    end

                    if isequal(marked_response,'False Alarm'); 
                        FA_starts(FA_counter) = (total_trial_starts_scan{i});
                        FA_duration(FA_counter) = trial_end; FA_code(FA_counter) = 5; FA_counter = FA_counter + 1;
                    end
             
                    if isequal(marked_response,'N/A');
                        no_response_starts(no_response_counter) = (total_trial_starts_scan{i});
                        no_response_duration(no_response_counter) = trial_end; no_response_code(no_response_counter) = 6; no_response_counter = no_response_counter + 1;
                    end

%% ---------------------- STORE TRIAL INFORMATION ---------------------- %%         

                    % Convert values to milliseconds and add to relevant cell arrays. 
                    task_condition{i} = trial_condition; trial_offsets{i} = trial_offset * 1000; trial_fixations{i} = trial_fixation * 1000; 
                    trial_fixation_pres{i} = trial_fixation_presentation * 1000; target_RFPs{i} = target_RFP * 1000;
                    target_RFP_pres{i} = target_RFP_presentation * 1000; response_fixations{i} = response_fixation * 1000;
                    response_fixation_pres{i} = response_fixation_presentation * 1000; inter_block_intervals{i} = inter_block_interval * 1000;
                    inter_block_interval_pres{i} = inter_block_interval_presentation * 1000; trial_ends{i} = trial_end * 1000;
                    cumulative_trial_starts{i} = cumulative_trial_start * 1000; cumulative_trial_ends{i} = trial_end_in_scan * 1000;
                    final_trial_ends{i} = final_trial_end * 1000;

                    % Add variables to relevant arrays and increment counter.
                    trial_number{i} = i; reference_orientation{i} = phi; reference_contrast{i} = C; reference_ram{i} = A; response_made{i} = response;
                    reaction_time{i} = RT; change_count{i} = num_change; change_type{i} = lower(type_changes); conjunction{i} = co_occur; 
                    response_required{i} = response_needed; i = i + 1; 

                    cd(data_storage);   % Change directory to data storage directory (specified above). 
                    
                    if exist(pid_directory, 'dir') > 1;    % If the unique participant folder already exists: 
                        cd(full_pid_directory);                 % Index into this directory.
                    else                                    % Else, if this directory does not already exist:
                        mkdir(pid_directory); cd(full_pid_directory);  % First create, then index into this directory. 
                    end 
                      
                    if exist(full_scan_directory,'dir');            % If the unique scan folder already exists:
                        cd(scan_dir);                                   % Index into this directory. 
                    elseif exist(full_scan_directory,'dir') == 0;   % Else, if the unique scan folder does not already exist:
                        mkdir(scan_dir); cd(scan_dir);                  % First create, then index into this directory. 
                    end

                    save(matlab_data); % Save all workspace variables created up to this point in loop.
                    
                    cd(strcat(root_directory, '\supporting_files'));    % Index back into function directory - Shuttle.
%                   cd(strcat(root_directory,'/supporting_files'));   % Index back into function directory - Mac. 
 
                    pp = pp + 1;    % Increment counter variable (ensures loop is only completed once whilst waiting for trial duration to be reached).         
                    
                end
            end

            real_trial_end = GetSecs;                               % Record accurate trial end time. 
            real_trial_end = real_trial_end - ind_trial_start;      % Set accurate trial end time minus time at start of trial (gets trial length).   
        end
      
%% ---------- INTER-BLOCK INTERVAL AND TRIAL CUE PRESENTATION ---------- %%

        if k == trials_desired;     % Else if trial number is equal to total number of trials desired per block:
            
            DrawFormattedText(window, 'X', 'center', 'center', black);    % Draw black fixation cross. 

            % Specify end time as total trial with inter-block interval length + trial start time. 
            endTime = total_inter_block_interval_trial_length + ind_trial_start;   

            Screen('Flip',window); inter_block_interval = GetSecs(); % Flip inter-block interval fixation to screen and record time.

            while GetSecs < endTime;                            % While current time is less than specified end time:
                if block_number == length(block_order);            % If block number is equal to the total number of blocks desired:
                    
                    % Draw black fixation cross (extend inter_block_interval to 9s of fixation cross rather than 7.5s fixation cross and 1.5s trial cue). 
                    DrawFormattedText(window, 'X', 'center', 'center', black);    
                    
                elseif block_number < length(block_order);         % Else if block number is less than total number of blocks desired:
                    instructions = sprintf('%s',trial_cue_2);           % Set trial cue to upcoming block condition.  
                    DrawFormattedText(window, sprintf(instructions), 'center', 'center', white);  % Draw white trial cue text.   
                end
            end
               
            % Specify end time as total trial with inter-block interval length + trial start time + trial cue length. 
            endTime = total_inter_block_interval_trial_length + ind_trial_start + trial_cue_length;
            
            Screen('Flip', window);  trial_cue_interval = GetSecs();   % Flip trial cue to screen. 
            
            while GetSecs < endTime;    % Wait until current time is equal to specified end time. 
                
                if ppp < 2; % Whle counter variable is less than 2 (ensures loop is only completed once). 
                
                % Record time at end of each trial (to record total between block interval time). 
                    % Again, this is not an accurate indication, but is required to create log files. See excel output for accurate timings. 
                    final_trial_end = GetSecs();    
                    final_trial_end = final_trial_end - ind_trial_start;

                    trial_end_in_block = trial_end_in_block + final_trial_end;      % Calculate total duration of trials so far within a block (refreshed each block).
                    trial_end_in_scan = trial_end_in_scan + final_trial_end + total_trial_length;        % Calculate total duration of trials so far within a scan.

                    total_trial_start_scan = total_trial_start_scan + trial_end;       % Set variable to time taken for the first inter-block fixation and trial cue. 

                    % Add these values to the relevant cell arrays. 
                    total_trial_starts_scan{i} = total_trial_start_scan; total_trial_ends_scan{i} = trial_end_in_scan * 1000;

                    % If i is a multiple of trials desired (i.e. if this is the last trial of a block):
                    if i == trials_desired * oo;

                        % Add inter-block interval and trial cue durations to the value in total trial start scan (will be used in next loop iteration). 
                            % i.e. this ensures trial responses do not appear during inter-block or trial cue periods. 
                        total_trial_start_scan = total_trial_start_scan + baseline_length + trial_cue_length;   
                    end

                    inter_block_interval = inter_block_interval - (ind_trial_start + trial_offset); % Record time inter-block-interval fixation was presented.
                    trial_cue_interval = trial_cue_interval - (ind_trial_start + trial_offset);     % Record time trial cue was presented. 
                    inter_block_interval_presentation = trial_cue_interval;                         % Calculate presentation period of inter-block interval. 
                    response_fixation_presentation = inter_block_interval - response_fixation;      % Calculate presentation period of response fixation from inter-block interval.                                              
                    ref_RFP_presentation = trial_fixation - trial_start_ref_RFP;                    % Calculate presentation period of reference RFP (period it was visible).  

                    % Convert reference RFP start and presentation values to ms and add to relevant cell arrays. 
                    reference_RFPs{i} = trial_start_ref_RFP * 1000; reference_RFP_pres{i} = ref_RFP_presentation * 1000;

                    % Calculate presentation periods for fixation between reference and target RFPs and target RFP. 
                    trial_fixation_presentation = target_RFP - trial_fixation; target_RFP_presentation = response_fixation - target_RFP;

    %% ------- CREATION OF FSL AND FREESURFER/MR VISTA TIMING ARRAYS ------- %%

                    % If trial number is equal to the total number of trials desired per condition:
                    if isequal(block_label{block_order(block_number)}, 'ORIENTATION');        % If the condition was an orientation attention condition:

                        if y == 1;                                                                        % If this is the first condition of the scan:

                            % Set time of condition onset to time at end of first between-block interval. 
                            orientation_starts(orientation_counter) = baseline_length + trial_cue_length;       

                        else        % Otherwise:                                                                

                         % Set condition start time to time at end of last condition (from cumulative trial scan end counter) -  length of 1 block.     
                            % total_trial_ends_scan gives time at end of the block in question, so minus block length to find start time (in ms). 
                            orientation_starts(orientation_counter) = (total_trial_ends_scan{trials_desired*y} - (trials_desired*total_trial_length)*1000) /1000; 
                        end

                        orientation_duration(orientation_counter) = round(total_trial_length * trials_desired,1); % Specify orientation duration as total number of trials desired * trial length. 
                        orientation_code(orientation_counter) = 2;  % Specify orientation code. 
                        orientation_counter = orientation_counter + 1;    % Increment orientation condition counter variable. 

                                          % Repeat this process for conditions with contrast, shape and a passive attentional focus.

                       elseif isequal(block_label{block_order(block_number)}, 'CONTRAST');
                           if y == 1;
                            contrast_starts(contrast_counter) = baseline_length + trial_cue_length;  
                          else
                            contrast_starts(contrast_counter) = (total_trial_ends_scan{trials_desired*y} - (trials_desired*total_trial_length)*1000)/1000; 
                          end

                          contrast_duration(contrast_counter) = round(total_trial_length * trials_desired,1);
                          contrast_code(contrast_counter) = 3;
                          contrast_counter = contrast_counter + 1;

                      elseif isequal(block_label{block_order(block_number)}, 'SHAPE');
                          if y == 1;
                            shape_starts(shape_counter) = baseline_length + trial_cue_length;  
                          else
                            shape_starts(shape_counter) = (total_trial_ends_scan{trials_desired*y} - (trials_desired*total_trial_length)*1000)/1000;
                          end

                          shape_duration(shape_counter) = round(total_trial_length * trials_desired,1);
                          shape_code(shape_counter) = 4;
                          shape_counter = shape_counter + 1;

                      elseif isequal(block_label{block_order(block_number)}, 'PASSIVE');
                          if y == 1;
                            passive_starts(passive_counter) = baseline_length + trial_cue_length;  
                          else
                            passive_starts(passive_counter) = (total_trial_ends_scan{trials_desired*y} - (trials_desired*total_trial_length)*1000)/1000; 
                          end

                          passive_duration(passive_counter) = round(total_trial_length * trials_desired,1);
                          passive_code(passive_counter) = 5;
                          passive_counter = passive_counter + 1;
                    end

                    y = y + 1;   % Increment condition counter variable. 


                    % If trial number is less than total number of trials specified for the entire scan and the trial number is greater than 1:
                    if i < trials_desired * (length(blocks_desired) * block_repeats_desired) && i > 1;

                        % Inter-block start set to cumulative end time at end of inter-block interval, so we subtract the length of the inter-block period to gain inter-block start time.
                        inter_block_starts(between_block_counter+1) = (total_trial_ends_scan{trials_desired*between_block_counter} - baseline_length)/1000;

                        % Set end time of inter-block fixation to start time + length of inter-block fixation. 
                        inter_block_ends_time(between_block_counter+1) = inter_block_starts(between_block_counter+1) + baseline_length;

                        % Set inter-block duration to length specifed by baseline_length variable (s) set inter-block code.
                        inter_block_duration(between_block_counter+1) = baseline_length; inter_block_code(between_block_counter+1) = 0;

                        % Set trial cue start time to time at end of previous inter-block fixation.
                        cue_starts(between_block_counter+1) = inter_block_ends_time(between_block_counter+1);

                        % Specify cue duration as time in trial_cue_length variable (s) and set trial cue code;
                        cue_duration(between_block_counter+1) = trial_cue_length; cue_code(between_block_counter+1) = 1;

                        between_block_counter = between_block_counter + 1; % Increment between-block counter variable.

                    % Else, if trial number is equal to the total number of trials specified for the entire scan:
                    elseif i == trials_desired * (length(blocks_desired) * block_repeats_desired);

                        % Inter-block start set to cumulative end time at end of inter-block interval, so we subtract the length of the inter-block period to gain inter-block start time.
                        inter_block_starts(between_block_counter+1) = (total_trial_ends_scan{trials_desired*between_block_counter} - baseline_length)/1000;

                        % Set inter-block duration to total length of between block interval (no trial cue is presented).
                        inter_block_duration(between_block_counter+1) = baseline_length + trial_cue_length; 

                        % Specify inter-block ends time to start time + length of between block interval.
                        inter_block_ends_time(between_block_counter+1) = inter_block_starts(between_block_counter+1) + baseline_length + trial_cue_length; 
                        inter_block_code(between_block_counter+1) = 0;

                        between_block_counter = between_block_counter + 1; % Increment between-block counter variable. 
                    end  

                    if isequal(type_changes,'no_change');                                   % If no change between reference and target RFP values occured:

                        % Specify start of trial as cumulative trial counter time (not updated until lower in the loop). 
                        no_change_starts(no_change_counter) = total_trial_starts_scan{i};       
                        no_change_duration(no_change_counter) = trial_end;                      % Specify no change duration as length of singular trial. 

                        % Specify no change trial code and increment counter variable.   
                        no_change_code(no_change_counter) = 5; no_change_counter = no_change_counter + 1;   
                    end

                                        % Repeat this process for trials including orientation, contrast and shape changes.  
                    if strfind(type_changes,'ORIENTATION') > 0;
                        orientation_change_starts(orientation_change_counter) = total_trial_starts_scan{i};
                        orientation_change_duration(orientation_change_counter) = trial_end;
                        orientation_change_code(orientation_change_counter) = 2; orientation_change_counter = orientation_change_counter + 1;
                     end

                    if strfind(type_changes,'CONTRAST') > 0;
                        contrast_change_starts(contrast_change_counter) = total_trial_starts_scan{i};
                        contrast_change_duration(contrast_change_counter) = trial_end;
                        contrast_change_code(contrast_change_counter) = 3; contrast_change_counter = contrast_change_counter + 1;
                    end

                    if strfind(type_changes,'SHAPE') > 0;
                        shape_change_starts(shape_change_counter) = total_trial_starts_scan{i};
                        shape_change_duration(shape_change_counter) = trial_end;
                        shape_change_code(shape_change_counter) = 4; shape_change_counter = shape_change_counter + 1;
                    end

                    if isequal(marked_response,'Hit');  % If participants response was marked as a hit:

                        % Specify the time at the start of this trial as cumulative trial counter time (not updated until lower in the loop).
                        hit_starts(hit_counter) = (total_trial_starts_scan{i}); 
                         % Specify duration as length of singular trial, set response code & increment counter variable. 
                        hit_duration(hit_counter) = trial_end; hit_code(hit_counter) = 2; hit_counter = hit_counter + 1;
                    end

                                              % Repeat this process for trials with miss, correct rejection, false alarm and no response data. 
                    if isequal(marked_response,'Miss');
                        miss_starts(miss_counter) = (total_trial_starts_scan{i});
                        miss_duration(miss_counter) = trial_end; miss_code(miss_counter) = 3; miss_counter = miss_counter + 1;
                    end
                    
                    if isequal(marked_response,'Correct Rejection');
                        CR_starts(CR_counter) = (total_trial_starts_scan{i});
                        CR_duration(CR_counter) = trial_end; CR_code(CR_counter) = 4; CR_counter = CR_counter + 1;
                    end

                    if isequal(marked_response,'False Alarm'); 
                        FA_starts(FA_counter) = (total_trial_starts_scan{i});
                        FA_duration(FA_counter) = trial_end; FA_code(FA_counter) = 5; FA_counter = FA_counter + 1;
                    end

                    if isequal(marked_response,'N/A');
                        no_response_starts(no_response_counter) = (total_trial_starts_scan{i});
                        no_response_duration(no_response_counter) = trial_end; no_response_code(no_response_counter) = 6; no_response_counter = no_response_counter + 1;
                    end

    %% ---------------------- STORE TRIAL INFORMATION ---------------------- %%         

                    % Convert values to milliseconds and add to relevant cell arrays. 
                    task_condition{i} = trial_condition; trial_offsets{i} = trial_offset * 1000; trial_fixations{i} = trial_fixation * 1000; 
                    trial_fixation_pres{i} = trial_fixation_presentation * 1000; target_RFPs{i} = target_RFP * 1000;
                    target_RFP_pres{i} = target_RFP_presentation * 1000; response_fixations{i} = response_fixation * 1000;
                    response_fixation_pres{i} = response_fixation_presentation * 1000; inter_block_intervals{i} = inter_block_interval * 1000;
                    inter_block_interval_pres{i} = inter_block_interval_presentation * 1000; trial_ends{i} = trial_end * 1000;
                    cumulative_trial_starts{i} = cumulative_trial_start * 1000; cumulative_trial_ends{i} = trial_end_in_scan * 1000;
                    final_trial_ends{i} = (final_trial_end + total_trial_length) * 1000 ;

                    % Add variables to relevant arrays and increment counter.
                    trial_number{i} = i; reference_orientation{i} = phi; reference_contrast{i} = C; reference_ram{i} = A; response_made{i} = response;
                    reaction_time{i} = RT; change_count{i} = num_change; change_type{i} = lower(type_changes); conjunction{i} = co_occur; 
                    response_required{i} = response_needed; i = i + 1; 

                    cd(data_storage);   % Change directory to data storage directory (specified above).  

                    if exist(full_pid_directory, 'dir');    % If the unique participant folder already exists: 
                        cd(full_pid_directory);                 % Index into this directory.

                    else                                    % Else, if this directory does not already exist:
                        mkdir(pid_directory); cd(full_pid_directory);  % First create, then index into this directory. 
                    end 

                    if exist(full_scan_directory,'dir');            % If the unique scan folder already exists:
                        cd(scan_dir);                                   % Index into this directory. 
                    elseif exist(full_scan_directory,'dir') == 0;   % Else, if the unique scan folder does not already exist:
                        mkdir(scan_dir); cd(scan_dir);                  % First create, then index into this directory. 
                    end

                    save(matlab_data); % Save all workspace variables created up to this point in loop.

                    cd(strcat(root_directory, '\supporting_files'));    % Index back into function directory - Shuttle.
%                   cd(strcat(root_directory,'/supporting_files'));     % Index back into function directory - Mac. 

                    ppp = ppp + 1; % Increment counter variable (ensures loop is only completed once whilst waiting to reach specified duration). 
                end
            end
            
            oo = oo + 1; % Increment counter variable (Increases by one every block completion). 
        end
        
        actual_trial_end = GetSecs();   % Record accurate time for end of trial with inter-block period and trial cue. 
        actual_trial_end = actual_trial_end - ind_trial_start;  % Subtract time at start of this trial to get trial duration. 
            
        if k ~= trials_desired;         % If the trial is not equal to the number of trials desired:
            % Add the accurate trial duration (without inter-block) to the relevant array (in ms).
            accurate_timing_ends{i-1} = real_trial_end * 1000;       
            
        else                            % Else, if the trial number is equal to the number of trials desired:
            % Add the accurate trial duration (with inter-block) to the relevant array (in ms).
            accurate_timing_ends{i-1} = actual_trial_end * 1000;        
        end
        
        accurate_trial_ends{i-1} = real_trial_end * 1000;   % Add accurate trial end (without inter-block) to accurate trial end counter. 
        
        % Calculate cumulative trial duration with each trial repetition (includes trials with inter-block interval).
        cumulative_trial_end = cumulative_trial_end + accurate_timing_ends{i-1} ;     
        accurate_cumulative_ends{i-1} = cumulative_trial_end;   % Store this cumulative trial duration value in relevant cell array. 
        
%         jheapcl;
    end
end

% release KbQueue
KbQueueRelease;

% % % RestrictKeysForKbCheck([]); % Re-enable input from all keypresses (not just keys corresponding to button-box keypresses. 
sca; ShowCursor; % Close screen and re-enable keypresses and mouse input. 
%ListenChar(); 
cd(data_storage); cd(full_pid_directory);  % Change directory to data storage and unique participant directories (specified above). 

%% --------------------- BLOCK TIMING INFORMATION ---------------------- %%

              %% ------------- FREESURFER ------------- %%

cd(data_storage); cd(full_pid_directory); cd(scan_dir);    % Move into unique participant scan directory.

% Specify the name of the freesurfer directory. 
freesurfer_directory = (sprintf('freesurfer_log_files')); 

% full_freesurfer_directory = strcat(data_storage,'/',pid_directory,'/',scan_dir,'/',freesurfer_directory); % Specify full directory name of freesurfer directory - Mac. 
full_freesurfer_directory = strcat(data_storage,'\',pid_directory,'\',scan_dir,'\',freesurfer_directory);   % Specify full directory name of freesurfer directory - Shuttle. 

if exist(full_freesurfer_directory, 'dir');                 % If the participant's freesurfer directory already exists:
    cd(freesurfer_directory);                                   % Index into this directory. 
else                                                        % Else, if this directory does not already exist:
    mkdir(freesurfer_directory); cd(freesurfer_directory);      % Create, and then index into this directory. 
end 

% Create matrices of block timing information, and create a final matrix containing all block timing information across scan. 
total_inter_block = [floor(inter_block_starts * 2)/2; inter_block_code; inter_block_duration; ones(1,length(inter_block_duration))]';
total_cue = [floor(cue_starts * 2)/2; cue_code; cue_duration; ones(1,length(cue_duration))]';
total_orientation_block = [floor(orientation_starts * 2)/2; orientation_code; orientation_duration; ones(1,length(orientation_duration))]';
total_contrast_block = [floor(contrast_starts * 2)/2; contrast_code; contrast_duration; ones(1,length(contrast_duration))]';
total_shape_block = [floor(shape_starts * 2)/2; shape_code; shape_duration; ones(1,length(shape_duration))]';
total_passive_block = [floor(passive_starts * 2)/2; passive_code; passive_duration; ones(1,length(passive_duration))]';
all_block_data_freesurfer = [total_inter_block; total_cue; total_orientation_block; total_contrast_block; total_shape_block; total_passive_block];

all_block_data_freesurfer = sortrows(all_block_data_freesurfer);  % Sort timings within this matrix into time order (starting at 0 to end of scan). 

% Save this data to a tab-delimited text file with 2 sig. figures precision. 
if pid < 10;
    dlmwrite(sprintf('freesurfer_block_log_pp_0%d_scan_%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_block_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
else
    dlmwrite(sprintf('freesurfer_block_log_pp_%d_scan_%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_block_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
end

              %% ------------- MR VISTA ------------- %%
              
cd(data_storage); cd(full_pid_directory); cd(scan_dir);    % Move into unique participant scan directory.

% Specify the name of the mr_vista directory. 
mrvista_directory = sprintf('mrvista_log_files'); 

% full_mrvista_directory = strcat(data_storage,'/',pid_directory,'/',scan_dir,'/',mrvista_directory); % Specify full directory name of mr vista directory - Mac. 
full_mrvista_directory = strcat(data_storage,'\',pid_directory,'\',scan_dir,'\',mrvista_directory);     % Specify full directory name of mr vista directory - Shuttle. 

if exist(full_mrvista_directory, 'dir');                        % If the participant's mr vista directory already exists:
    cd(full_mrvista_directory);                                     % Index into this directory. 
else                                                            % Else, if this directory does not already exist:
    mkdir(full_mrvista_directory); cd(full_mrvista_directory);      % Create, and then index into this directory. 
end 
     
% Repeat the same process for mr vista block timings. 
total_inter_block_mrvista = [floor(inter_block_starts * 2)/2; inter_block_code]';
total_cue_mrvista = [floor(cue_starts * 2)/2; cue_code]';
total_orientation_mrvista = [floor(orientation_starts * 2)/2; orientation_code]';
total_contrast_mrvista = [floor(contrast_starts * 2)/2; contrast_code]';
total_shape_mrvista = [floor(shape_starts * 2)/2; shape_code]';
total_passive_mrvista = [floor(passive_starts * 2)/2; passive_code]';
all_block_data_mrvista = [total_inter_block_mrvista; total_cue_mrvista; total_orientation_mrvista; total_contrast_mrvista; total_shape_mrvista; total_passive_mrvista];

all_block_data_mrvista = sortrows(all_block_data_mrvista);  % Sort timings within this matrix into time order (starting at 0 to end of scan). 

% Save this data to a tab-delimited text file with 2 sig. figures precision. 
if pid < 10;
    dlmwrite(sprintf('mrvista_block_log_pp_0%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_block_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');    
else
    dlmwrite(sprintf('mrvista_block_log_pp_%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_block_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');    
end
%% -------------------- EVENT TIMING INFORMATION ----------------------- %%

% Establish base matrices with inter-block and cue timing information.
all_event_data_freesurfer = [total_inter_block; total_cue];   
all_event_data_mrvista = [total_inter_block_mrvista; total_cue_mrvista];
              
if no_change_starts ~= 0;   % If at least one trial in scan had no change between reference and target RFP values:
    
    % Collate timing and coding information for these trials into one overall matrix. 
    total_no_change_freesurfer = [no_change_starts; no_change_code; no_change_duration; ones(1,length(no_change_duration))]';  
    total_no_change_mrvista = [no_change_starts; no_change_code]';
    
    % Add this matrix to the overall matrix storing all event timing data. 
    all_event_data_freesurfer = [all_event_data_freesurfer; total_no_change_freesurfer];   
    all_event_data_mrvista = [all_event_data_mrvista; total_no_change_mrvista];
end
                             % Repeat this process for trials including orientation, contrast and shape changes.  
                             
if orientation_change_starts ~= 0;
   total_orientation_change_freesurfer = [orientation_change_starts; orientation_change_code; orientation_change_duration; ones(1,length(orientation_change_duration))]';
   total_orientation_change_mrvista = [orientation_change_starts; orientation_change_code]';
   all_event_data_freesurfer = [all_event_data_freesurfer; total_orientation_change_freesurfer];
   all_event_data_mrvista = [all_event_data_mrvista; total_orientation_change_mrvista];
end

if contrast_change_starts ~= 0;
    total_contrast_change_freesurfer = [contrast_change_starts; contrast_change_code; contrast_change_duration; ones(1,length(contrast_change_duration))]';
    total_contrast_change_mrvista = [contrast_change_starts; contrast_change_code]';
    all_event_data_freesurfer = [all_event_data_freesurfer; total_contrast_change_freesurfer];
    all_event_data_mrvista = [all_event_data_mrvista; total_contrast_change_mrvista];
end

if shape_change_starts ~= 0;
    total_shape_change_freesurfer = [shape_change_starts; shape_change_code; shape_change_duration; ones(1,length(shape_change_duration))]';
    total_shape_change_mrvista = [shape_change_starts; shape_change_code]';
    all_event_data_freesurfer = [all_event_data_freesurfer; total_shape_change_freesurfer];
    all_event_data_mrvista = [all_event_data_mrvista; total_shape_change_mrvista];
end

cd(data_storage); cd(full_pid_directory);cd(scan_dir); cd(freesurfer_directory); % Move into particpant's freesurfer directory. 

% Sort and write out data to tab-delimited .txt file. 
all_event_data_freesurfer = sortrows(all_event_data_freesurfer); 

if pid < 10;
    dlmwrite(sprintf('freesurfer_event_log_pp_0%d_scan_0%d_%s.txt',pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_event_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
else
    dlmwrite(sprintf('freesurfer_event_log_pp_%d_scan_0%d_%s.txt',pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_event_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
end
 
cd(data_storage); cd(full_pid_directory);cd(scan_dir); cd(mrvista_directory);   % Move into particpant's mr vista directory. 

% Sort and write out data to tab-delimited .txt file. 
all_event_data_mrvista = sortrows(all_event_data_mrvista);   

if pid < 10;
    dlmwrite(sprintf('mrvista_event_log_pp_0%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')) ,all_event_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');   
else
     dlmwrite(sprintf('mrvista_event_log_pp_%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')) ,all_event_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');  
end

%% ------------------- RESPONSE TIMING INFORMATION --------------------- %%

% Establish base matrices with inter-block and cue timing information.
all_response_data_freesurfer = [total_inter_block; total_cue];  
all_response_data_mrvista = [total_inter_block_mrvista; total_cue_mrvista];

                     % If participant made a response marked as a hit for at least one trial:
if hit_starts ~= 0;
    
    % Collate timing and coding information for these trials into one overall matrix. 
    total_hit_freesurfer = [hit_starts; hit_code; hit_duration; ones(1,length(hit_duration))]';
    total_hit_mrvista = [hit_starts; hit_code]';
    
    % Add this matrix to the overall matrix storing all event timing data. 
    all_response_data_freesurfer = [all_response_data_freesurfer; total_hit_freesurfer];     
    all_response_data_mrvista = [all_response_data_mrvista; total_hit_mrvista];
end
                        % Repeat this process for trials with miss, correct rejection, false alarm and no response data.
                        
if miss_starts ~= 0;
    total_miss_freesurfer = [miss_starts; miss_code; miss_duration; ones(1,length(miss_duration))]';
    total_miss_mrvista = [miss_starts; miss_code]';
    all_response_data_freesurfer = [all_response_data_freesurfer; total_miss_freesurfer];
    all_response_data_mrvista = [all_response_data_mrvista; total_miss_mrvista];
end

if CR_starts ~= 0;
    total_CR_freesurfer = [CR_starts; CR_code; CR_duration; ones(1,length(CR_duration))]';
    total_CR_mrvista = [CR_starts; CR_code]';
    all_response_data_freesurfer = [all_response_data_freesurfer; total_CR_freesurfer];
    all_response_data_mrvista = [all_response_data_mrvista; total_CR_mrvista];
end

if FA_starts ~= 0;
    total_FA_freesurfer = [FA_starts; FA_code; FA_duration; ones(1,length(FA_duration))]';
    total_FA_mrvista = [FA_starts; FA_code]';
    all_response_data_freesurfer = [all_response_data_freesurfer; total_FA_freesurfer];
    all_response_data_mrvista = [all_response_data_mrvista; total_FA_mrvista];
end
    
if no_response_starts ~= 0;
    
    total_no_response_freesurfer = [no_response_starts; no_response_code; no_response_duration; ones(1,length(no_response_duration))]';
    total_no_response_mrvista = [no_response_starts; no_response_code]';
    all_response_data_freesurfer = [all_response_data_freesurfer; total_no_response_freesurfer];
    all_response_data_mrvista = [all_response_data_mrvista; total_no_response_mrvista];
end

cd(data_storage); cd(full_pid_directory); cd(scan_dir); cd(freesurfer_directory);   % Move into particpant's freesurfer directory.

% Sort and write out data to tab-delimited .txt file.
all_response_data_freesurfer = sortrows(all_response_data_freesurfer);

if pid < 10;
    dlmwrite(sprintf('freesurfer_response_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), all_response_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
else
    dlmwrite(sprintf('freesurfer_response_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), all_response_data_freesurfer,'delimiter','\t','precision','%.2f', 'newline', 'pc');
end
 
cd(data_storage); cd(full_pid_directory);cd(scan_dir); cd(mrvista_directory);       % Move into particpant's mr vista directory.

% Sort and write out data to tab-delimited .txt file.
all_response_data_mrvista = sortrows(all_response_data_mrvista);

if pid < 10;
    dlmwrite(sprintf('mrvista_response_log_pp_0%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_response_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');  
else
    dlmwrite(sprintf('mrvista_response_log_pp_%d_scan_0%d_%s.par', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')),all_response_data_mrvista,'delimiter','\t','precision','%.2f', 'newline', 'pc');  
end

%% ---------------------- FSL LOG FILE OUTPUT -------------------------- %%

cd(data_storage); cd(full_pid_directory);  cd(scan_dir);  % Move into unique participant scan directory. 

% Specify name of fsl directory. 
fsl_directory = (sprintf('fsl_log_files')); 

% full_fsl_directory = strcat(data_storage,'/', pid_directory, '/', scan_dir,'/',fsl_directory);  % Specify full directory name of fsl directory - Mac. 
full_fsl_directory = strcat(data_storage,'\', pid_directory, '\', scan_dir,'\',fsl_directory);    % Specify full directory name of fsl directory - Shuttle. 

if exist(full_fsl_directory, 'dir');                % If the participant's fsl directory already exists:
    cd(full_fsl_directory);                             % Index into this directory. 
else                                                % Else, if this directory does not already exist:
    mkdir(full_fsl_directory); cd(full_fsl_directory);  % Create, and then index into this directory. 
end 

%% --------------------- BLOCK TIMING INFORMATION ---------------------- %%

% Create matrix of block timing information.
orientation_blocks = [(floor(orientation_starts * 2)/2)', orientation_duration', ones(1,length(orientation_duration))'];
contrast_blocks = [(floor(contrast_starts * 2)/2)', contrast_duration', ones(1,length(contrast_duration))'];
shape_blocks = [(floor(shape_starts * 2)/2)', shape_duration', ones(1,length(shape_duration))'];
passive_blocks = [(floor(passive_starts * 2)/2)', passive_duration', ones(1,length(passive_duration))'];
inter_block_blocks = [(floor(inter_block_starts * 2)/2)', inter_block_duration', ones(1,length(inter_block_duration))'];
cue_blocks = [(floor(cue_starts * 2)/2)', cue_duration', ones(1,length(cue_duration))'];

% Save this data to a tab-delimited text file with 2 sig. figures precision. 
if pid < 10;
    dlmwrite(sprintf('FSL_orientation_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), orientation_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 

                                % Repeat this process for blocks with contrast, shape and passive attentional focus.
    dlmwrite(sprintf('FSL_contrast_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), contrast_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_shape_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), shape_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_passive_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), passive_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_inter_block_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), inter_block_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_cue_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), cue_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');

else
    dlmwrite(sprintf('FSL_orientation_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), orientation_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    dlmwrite(sprintf('FSL_contrast_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), contrast_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_shape_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), shape_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_passive_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), passive_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_inter_block_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), inter_block_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    dlmwrite(sprintf('FSL_cue_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), cue_blocks,'delimiter','\t','precision','%.2f', 'newline', 'pc');
end

%% -------------------- EVENT TIMING INFORMATION ----------------------- %%

if no_change_starts ~= 0;   % If at least one trial in scan had no change between reference and target RFP values:
    
    % Create matrix of orientation block timing information.
    no_change_trials = [round(no_change_starts,1)', no_change_duration', ones(1,length(no_change_duration))'];    
    
    % Save this data to a tab-delimited text file with 2 sig. figures precision. 
    if pid < 10;
        dlmwrite(sprintf('FSL_no_change_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), no_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');  
    else
        dlmwrite(sprintf('FSL_no_change_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), no_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');    
    end
end
                            % Repeat this process for trials including orientation, contrast and shape changes.  
                            
if orientation_change_starts ~= 0;
    orientation_change_trials = [round(orientation_change_starts,1)', orientation_change_duration', ones(1,length(orientation_change_duration))'];
    
    if pid < 10;
        dlmwrite(sprintf('FSL_orientation_change_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), orientation_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    else
        dlmwrite(sprintf('FSL_orientation_change_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), orientation_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    end
end

if contrast_change_starts ~= 0;
    contrast_change_trials = [round(contrast_change_starts,1)', contrast_change_duration', ones(1,length(contrast_change_duration))'];
    
    if pid < 10;
        dlmwrite(sprintf('FSL_contrast_change_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), contrast_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    else
         dlmwrite(sprintf('FSL_contrast_change_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), contrast_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    end
end

if shape_change_starts ~= 0;
    shape_change_trials = [round(shape_change_starts,1)', shape_change_duration', ones(1,length(shape_change_duration))'];
    
    if pid < 10;
        dlmwrite(sprintf('FSL_shape_change_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), shape_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    else
        dlmwrite(sprintf('FSL_shape_change_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), shape_change_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    end
end

%% ------------------- RESPONSE TIMING INFORMATION --------------------- %%

if hit_starts ~= 0; % If participant made a response marked as a hit for at least one trial:
    
    % Create matrix of orientation block timing information.
    hit_trials = [round(hit_starts,1)', hit_duration', ones(1,length(hit_duration))'];   
    
    % Save this data to a tab-delimited text file with 2 sig. figures precision. 
    if pid < 10;
        dlmwrite(sprintf('FSL_hit_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number,datestr(now,'ddmmyyyy_HHMM')), hit_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');    
    else
        dlmwrite(sprintf('FSL_hit_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), hit_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc'); 
    end
end
              % Repeat this process for trials with miss, correct rejection, false alarm and no responses.  
  
if miss_starts ~= 0;
    miss_trials = [round(miss_starts,1)', miss_duration', ones(1,length(miss_duration))'];    
    if pid < 10;
        dlmwrite(sprintf('FSL_miss_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), miss_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');
    else
        dlmwrite(sprintf('FSL_miss_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), miss_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');     
    end
end    
     
if CR_starts ~= 0;
    CR_trials = [round(CR_starts,1)', CR_duration', ones(1,length(CR_duration))']; 
    
    if pid < 10;
        dlmwrite(sprintf('FSL_CR_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), CR_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');     
    else
        dlmwrite(sprintf('FSL_CR_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), CR_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');     
    end
end    
       
if FA_starts ~= 0;
    FA_trials = [round(FA_starts,1)', FA_duration', ones(1,length(FA_duration))'];          
    
    if pid < 10;
        dlmwrite(sprintf('FSL_FA_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), FA_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');  
    else
        dlmwrite(sprintf('FSL_FA_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), FA_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');     
    end
end    
  
if no_response_starts ~= 0;
    no_response_trials = [round(no_response_starts,1)', no_response_duration', ones(1,length(no_response_duration))']; 
    
    if pid < 10;
        dlmwrite(sprintf('FSL_no_response_log_pp_0%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), no_response_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');    
    else
        dlmwrite(sprintf('FSL_no_response_log_pp_%d_scan_0%d_%s.txt', pid, scan_number, datestr(now,'ddmmyyyy_HHMM')), no_response_trials,'delimiter','\t','precision','%.2f', 'newline', 'pc');     
    end
end    

%% ------------------------ EXCEL DATA OUTPUT -------------------------- %%
cd(data_storage); cd(full_pid_directory); cd(scan_dir);      % Move into unique participant scan directory. 

% Store all response/feature data in one transposed matrix. 
all_data = [trial_number; condition; pp_orientation; reference_orientation; target_orientation; orientation_direction; ...
    orientation_change_tracker; pp_contrast; reference_contrast; target_contrast; contrast_direction; contrast_change_tracker; ...
    pp_ram; reference_ram; target_ram; ram_direction; ram_change_tracker; key_pressed; response_made; marked_responses;... 
    reaction_time; change_count; change_type; conjunction; response_required; condition]';          
        
% Convert response/feature data to table with headings.
data_table = cell2table(all_data, 'VariableNames', {'Trial' 'Condition' 'O_Threshold' 'O_Reference' 'O_Target' 'O_Direction' 'O_Tracker'...
    'C_Threshold' 'C_Reference' 'C_Target' 'C_Direction' 'C_Tracker' 'S_Threshold' 'S_Reference' 'S_Target' 'S_Direction' 'S_Tracker'...
    'Key_Press' 'Response' 'Detection' 'RT' 'Change_Count' 'Change_Type' 'Conjunction' 'Change_Response_Required' 'Condition_2'});  

if task == 'T';         %If temporal task:
    
    % Store all timing data in one transposed matrix (temporal stores reference RFP and inter-trial fixation times). 
    all_times = [trial_number; task_condition; cumulative_trial_starts; trial_offsets; reference_RFPs; reference_RFP_pres; trial_fixations; trial_fixation_pres;...
        target_RFPs; target_RFP_pres; response_fixations; response_fixation_pres; accurate_trial_ends; inter_block_intervals; inter_block_interval_pres; accurate_timing_ends;...
        accurate_cumulative_ends]';    
    
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames', {'Trial_Number' 'Condition' 'C_Trial_Start' 'Trial_Offset' 'Trial_Start'...
                                                        'Reference_RFP_Presentation' 'Trial_Fixation' 'Trial_Fixation_Presentation' 'Target_RFP'...
                                                        'Target_RFP_Presentation' 'Response_Fixation' 'Response_Fixation_Presentation' 'Trial_End'...
                                                        'Inter_Block_Interval' 'Inter_Block_Interval_Presentation' 'Final_Trial_End' 'C_Trial_End'});  
elseif task == 'S';     % If spatial task:

    % Store all timing data in one transposed matrix (spatial does not store reference RFP or inter-trial fixation times). 
    all_times = [trial_number; task_condition; cumulative_trial_starts; trial_offsets; target_RFPs; target_RFP_pres;...
        response_fixations; response_fixation_pres; trial_ends; inter_block_intervals; inter_block_interval_pres; ...
        cumulative_trial_ends]';
    
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames',{'Trial_Number' 'Condition' 'C_Trial_Start' 'Trial_Offset' 'Trial_Start'...
                                                        'RFP_Presentation' 'Response_Fixation' 'Response_Fixation_Presentation' 'Trial_End'...
                                                        'Inter_Block_Interval' 'Inter_Block_Interval_Presentation' 'C_Trial_End_In_Block'});
end

if pid < 10;
    writetable(data_table, sprintf('fmri_responses_%s_pp_0%d_scan_0%d_%s.csv', condition_task, pid, scan_number, datestr(now,'ddmmyyyy_HHMM')));      % Save response data as .csv file.
    writetable(time_table, sprintf('fmri_timing_log_%s_pp_0%d_scan_0%d_%s.csv', condition_task, pid, scan_number, datestr(now,'ddmmyyyy_HHMM')));     % Save timing data as .csv.
else
    writetable(data_table, sprintf('fmri_responses_%s_pp_%d_scan_%d_%s.csv', condition_task, pid, scan_number, datestr(now,'ddmmyyyy_HHMM')));      % Save response data as .csv file.
    writetable(time_table, sprintf('fmri_timing_log_%s_pp_%d_scan_%d_%s.csv', condition_task, pid, scan_number, datestr(now,'ddmmyyyy_HHMM')));     % Save timing data as .csv.
end
   
if pid < 10;
    matlab_data = sprintf('fmri_temporal_pp_0%d_scan_0%d_data_%s', pid, scan_number,datestr(now,'ddmmyyyy_HHMM'));
else
    matlab_data = sprintf('fmri_temporal_pp_%d_scan_0%d_data_%s', pid, scan_number,datestr(now,'ddmmyyyy_HHMM'));
end
save(matlab_data);

cd(start_dir);  % Move back into start directory (directory code is stored in). 
% % % jheapcl;
%% --------------------------------------------------------------------- %%