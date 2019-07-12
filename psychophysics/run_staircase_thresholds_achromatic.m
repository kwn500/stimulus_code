%% ------------------- RFP 2AFC DISCRIMINATION TASK -------------------- %%

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

% Participant must judge whether target RFP is different in a cued feature 
% with respect to a reference RFP. Participant must respond with 'A' or 
% 'L' button presses indicating the direction of change. 

% Run orientation 'feature_task' condition first. 

%% ---------------- PRE-ALLOCATE CELL ARRAYS FOR SPEED ----------------- %%  

sca; close all ; clearvars; jheapcl;
KbName('UnifyKeyNames');   
start_dir = pwd;

trial_number = cell(1); test_values = cell(1); reference_value = cell(1);  
movement_direction = cell(1); responses = cell(1); reaction_time = cell(1); 
answers = cell(1); key_pressed = cell(1); movement_value = cell(1);
reference_RFPs = cell(1); reference_RFP_pres = cell(1);
task_condition = cell(1); trial_offsets = cell(1); trial_fixations = cell(1);
trial_fixation_pres = cell(1); target_RFPs = cell(1); target_RFP_pres = cell(1);
response_fixations = cell(1); response_fixation_pres = cell(1); 
trial_ends = cell(1); cumulative_trial_starts = cell(1); cumulative_trial_ends = cell(1);

%% ------------------------ GENERAL SETUP ------------------------------ %%

% CHECK THESE SETTINGS MATCH PIXELS PER DEGREE SETTINGS IN RFP_GENERATOR & SELECTIVE VS NON-SELECTIVE. 
% Calculate number of pixels per degree to control later stimulus sizes. 
display_width = 52.4;       % Monitor width (cm).
display_horRes = 1920;      % Horizontal resolution of monitor (pixels).
display_viewDist = 57;      % Participant viewing distance (cm).
display_ppd = display_horRes/(2*atand((display_width/2)/display_viewDist)); % Calculates pixels/°. 

% Create and present dialog box that gathers user input.
prompt = {'Enter Participant Number:','Enter Scan Number:', 'Enter Feature Task:' ...
    'Enter Threshold Estimate'};
dlg_title = 'Participant Information';
num_lines = 1;
answer = inputdlg(prompt,dlg_title,num_lines);

% Set participant id, scan number and feature task variables from GUI input. 
pid = str2double(answer{1}); scan_number = str2double(answer{2}); feature_task = answer{3};

% Sets up directories for data storage and access to functions. 
% root_directory = ('/Users/Kirstie/Documents/rfp_code_data/code');    % Top-level directory.
% base_directory = strcat(root_directory,'/supporting_files');            % Path to relevant functions.

% Sets up directories for data storage and access to functions. 
root_directory = ('\Users\PS-RMLab\Documents\Kirstie\updated');    % Top-level directory- Shuttle.
base_directory = strcat(root_directory,'\supporting_files');             % Path to relevant functions- Shuttle.
cd(base_directory);                                                     % Move into this function directory.        

% Vital experimental information - check this is correct before running!. 
trials_desired = 50;                 % Number of trials.
task = 'T';                         % Type of task - 'T' = temporal, 'S' = spatial. 
condition_task = 'temporal';        % Task condition string (for data output)- 'temporal'/'spatial'.  

practice_trials = 10;
stimulus_presentation_time = 0.2;   % Stimulus presentation time (s). 
response_presentation_time = 0.8;   % Response period duration (s). 
total_trial_length = 1.5;           % Total trial duration (s).  

data_storage = strcat(root_directory,'/data_output/staircase_thresholds/');

% Set fixation characteristics for each condition and specific data storage directory. 
if isequal(feature_task,'Orientation')         % If orientation threshold task. 
    fixation = 'O'; % Set relevant fixation letter and colour. 
%     data_storage = strcat(root_directory,'\data_output\staircase_thresholds\orientation_discrimination/'); % Create directory for storage of thresholds.
%     data_storage = strcat(root_directory,'/data_output/staircase_thresholds/orientation_discrimination/');  %mac
%     condition_dir = 'orientation_discrimination/';
    condition_dir = 'orientation_discrimination\';  % Shuttle. 
elseif isequal(feature_task,'Contrast')
    fixation = 'C';  
%     data_storage = strcat(root_directory,'\data_output\staircase_thresholds\contrast_discrimination/');
%     data_storage = strcat(root_directory,'/data_output/staircase_thresholds/contrast_discrimination/');  %mac
%     condition_dir = 'contrast_discrimination/';
    condition_dir = 'contrast_discrimination\';   % Shuttle. 
elseif isequal(feature_task, 'RAM')
    fixation = 'S'; 
%     data_storage = strcat(root_directory,'\data_output\staircase_thresholds\RAM_discrimination/');
%     data_storage = strcat(root_directory,'/data_output/staircase_thresholds/RAM_discrimination/');  %mac
%     condition_dir = 'RAM_discrimination/';
    condition_dir = 'RAM_discrimination\';    % Shuttle. 
end

% General experiment 'housekeeping'. 
threshold = zeros(1,trials_desired);                        % Create empty column array for storage of threshold values. 
space = KbName('space'); escape = KbName('escape');         % Set keynames for potential response keypresses.                                 

PsychDefaultSetup(2); screens = Screen('Screens');  screenNumber = max(screens);    % Call default PTB settings & draw to max screen.
black = BlackIndex(screenNumber); white = WhiteIndex(screenNumber);                 % Define black (0) & white (1).

%% -------------------------INITALISE QUEST ---------------------------- %%

% Collects threshold and SD estimates & establishes structure for later estimates. 
tGuess=str2double(answer{4}); 
% tGuessSd=str2double(answer{5});
tGuessSd = 0.5;

pThreshold=0.75;                                                                    % Set threshold for % correct responses. 
beta=3.5;delta=0.01;gamma=0.5;
% range=1;  
q = QuestCreate(log10(tGuess),log10(tGuessSd),pThreshold,beta,delta,gamma); 
% range);   % Create struct for measuring threshold.  
q.normalizePdf=1;                                                                        

%% -------------------------- SCREEN SETUP------------------------------ %%

% PsychDebugWindowConfiguration

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

%% --------------------------- SOUND SETUP ----------------------------- %%

InitializePsychSound(1);  % Initialize sounddriver
nrchannels = 2;  freq = 48000;   repetitions = 1; beepLengthSecs = 1; beepPauseTime = 1;  % Set key variables. 

% Wait for device to start & start sound immediately when called. 
startCue = 0; waitForDeviceStart = 1; 

pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels); % Initialise PsychPort audio. 

%% ----------------------- INTRODUCTION SCREEN ------------------------- %%

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial');    % Set text parameters. 

% Present relevant instructions for selected condition. 
line1 = sprintf('%s Discrimination Task', feature_task);

if task == 'T'
    target_shape = 'second shape'; reference_shape = 'first';
    line2 = '\n\n Compare the first and second shapes.';
elseif task == 'S'
    target_shape = 'right shape'; reference_shape = 'left';
end

if isequal(feature_task,'Orientation')    
    line4 = sprintf('\n\n If the %s is rotated ANTICLOCKWISE with respect the %s, press "N"', target_shape, reference_shape);
    line3 = sprintf('\n\n If the %s is rotated CLOCKWISE with respect to the %s, press "U"', target_shape, reference_shape);
elseif isequal(feature_task,'Contrast')
    line3 = sprintf('\n\n If the %s is HIGHER contrast with respect the %s, press "U"', target_shape, reference_shape);
    line4 = sprintf('\n\n If the %s is LOWER contrast with respect the %s, press "N"', target_shape, reference_shape);
elseif isequal(feature_task, 'RAM')
    line3 = sprintf('\n\n If the %s is SPIKIER (higher RAM) than the %s, press "U"', target_shape, reference_shape);
    line4 = sprintf('\n\n If the %s is SMOOTHER (lower RAM) than the %s, press "N"', target_shape, reference_shape);
end

line5 = '\n\n Press the spacebar to continue';  

DrawFormattedText(window, [line1 line2 line3 line4 line5], 'center', 'center', white);
Screen('Flip',window); [secs,KeyCode] = KbWait; % Draw text to screen and wait for keypress.  

% Wait for spacebar response to progress. 
while KeyCode(space) == 0  % While spacebar is not pressed.
    [secs,KeyCode] = KbWait;    % Wait for spacebar response. 
    if KeyCode(escape)         % If escape key pressed:
        sca; ListenChar(); ShowCursor;  % Close all, re-enable keypresses to command window and make mouse cursor visible. 
    end
end

%% ----------------------- RADIAL-FREQ PATTERN ------------------------- %%

% Set parameters for RFP. 
C = 0.5;	% Contrast.  
A = 0.5;    % Radial modulation amplitude (A < 1).  
phi = 0;    % Phi - angular rotation. 

% Create specific reference RFP for selected condition. 
if isequal(feature_task,'Orientation')
    reference = 0; 
    distances = RFP_generator(reference,C,A); 
elseif isequal(feature_task,'Contrast')
    reference = 0.5;
    distances = RFP_generator(phi,reference,A);
elseif isequal(feature_task, 'RAM')
    reference = 0.5; 
    distances = RFP_generator(phi,C,reference); 
end

% Converts specific reference RFP number matrix to texture.
scn = Screen('MakeTexture', window, distances);   

%% -------------------------- QUEST STAIRCASE -------------------------- %% 

% Set script to maximum processing priority. 
topPriorityLevel = MaxPriority(window); Priority(topPriorityLevel); 
    
goodBeep = MakeBeep(500, 0.2, 10000); badBeep= MakeBeep(500, 0.1, 100000);    % Set correct (high) and incorrect (low-pitched) beeps. 
Screen('TextSize', window, 15);                                               % Set smaller text size for fixation cross. 

% tic; total_trial_start = GetSecs();  % Record time at start of trials loop. 
trial_onset = 0;                     % Set individual trial start time to 0.   
cumulative_trial_end = 0;

for k = 1:trials_desired            % For each trial to the total number of trials desired. 
    
    ind_trial_start = GetSecs;       % Record time at start of each trial. 
    
    if k ==1 
        total_trial_start = ind_trial_start;
    end
    
    if k == 1                       % If this is the first trial:
        cumulative_trial_start = 0;     % Set start time of first trial to 0 (for cumulative trial time counter).   
    else 
        % Set cumulative time of trial to start time of individual trial - start time of block of trials. 
        cumulative_trial_start = ind_trial_start - total_trial_start; 
    end
     
    clear distances_target; % Clear the target RFP each loop, so a new one is generated each time. 

    response_made = 0; response = 'no response'; answer = 'no answer'; keyLabel = 'N/A'; RT = 0; key_press_counter = 0; % Set default response parameters. 

    if task == 'T'                                                            % If temporal task:
            Screen('DrawTexture', window, scn);                                       % Draw condition-specific reference RFP texture to screen. 
            DrawFormattedText(window, fixation, 'center', 'center', white); % Draws condition-specific fixation letter in condition colour. 
            
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

    while GetSecs < endTime                        % While current time is less than the specified end time:
        if exist('distances_target','var') == 0    % If the target RFP texture does not exist: 
           
        if k > 1 && k <= practice_trials
            tTest = tTest;
        else
            tTest=QuestQuantile(q);                         % Get the recommended threshold value from Quest (log10 units). 
        end        
        
        % Call 'calculate_value' function to calculate feature-value for test RFP and categorise this change (e.g. clockwise, lower contrast...).
        [change_value, test_value, ii, movement, max] = calculate_value(tTest, reference, feature_task);
        
        % Store these target RFP values in the relevant cell arrays. 
        test_values{k} = test_value; movement_value{k} = movement; threshold(k) = abs(change_value);
        
        % Create a target RFP with the change in the relevant feature dimension and create target RFP texture. 
        if isequal(feature_task,'Orientation')                 % e.g. If orientation task:                    
            distances_target = RFP_generator(test_value,C,A);       % Set phi equivalent to the change value.   
        elseif isequal(feature_task, 'Contrast')
            distances_target = RFP_generator(phi,test_value,A);    
        elseif isequal(feature_task, 'RAM')
            distances_target = RFP_generator(phi,C,test_value);  
        end
        sn = Screen('MakeTexture', window, distances_target);   % Create target RFP texture. 
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
            DrawFormattedText(window, fixation, 'center', 'center', white);       % Prepare the response period fixation cross. 
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
                    if isequal(keyLabel,'u') || isequal(keyLabel,'n')                  % If keypress was valid ('A' or 'L'):
                        response_made = 1;                                                  % Mark as valid response made. 
                    elseif isequal(keyLabel,'ESCAPE')                                  % Else if participant pressed the escape key. 
                        sca; ListenChar(); ShowCursor;                                      % Close the screen and re-enable keyboard and mouse input. 
                    end  
                end
            end    
   

        % If the keypress was not 'A' or 'L', or if no response was made, record default values in the relevant cell arrays. 
        if (isequal(keyLabel,'u') == 0 && isequal(keyLabel,'n') == 0) || response_made == 0    
                trial_number{k} = k; reference_value{k} = reference; responses{k} = response; ...
                    answers{k} = answer; key_pressed{k} = keyLabel; reaction_time{k} = RT * 1000; 
                continue                                                        % Move to next loop iteration- (Don't pass information to Quest). 
        end
    end
    
    response_fixation = response_fixation - (ind_trial_start + trial_offset);   % Record time response fixation was presented. 
    endTime = total_trial_length + ind_trial_start;           % Specify end time as total trial length + trial start time.        


    while GetSecs < endTime                                                    % While current time is less than the specified end time:
        if response_made == 1                                                     % If a valid response was made:
            while isequal(response,'no response')                                      % While the response hasn't been classified (is at default value):
                
            % Call 'response_classify' function to categorise response as correct/incorrect and set relevant feedback beep. 
             
            for p = 1:1             
                [response_value, response, myBeep, answer] = response_classify(keyLabel, ii, feature_task); % Ensures that each beep is only played once within while loop. 
                PsychPortAudio('FillBuffer', pahandle, [myBeep; myBeep]);                       % Fill sound buffer with relevant beep frequency. 
                PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);   % Play beep.
                PsychPortAudio('DeleteBuffer');
            end
                
            if k > practice_trials
                q=QuestUpdate(q,log10(abs(change_value)),response_value);                       % Update Quest with presented threshold and observer response. 
            end
            
            trial_number{k} = k; reference_value{k} = reference; responses{k} = response; ...   % Store trial variables in relevant arrays. 
                answers{k} = answer; key_pressed{k} = keyLabel; reaction_time{k} = RT * 1000; 
            end
        end

   
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
        
        scan_dir = sprintf('scan_%d', scan_number);

        if exist(scan_dir,'dir') == 0
            mkdir(scan_dir);
        end
        
        cd(scan_dir);
        
        if exist(condition_dir,'dir') == 0
            mkdir(condition_dir);
        end
        
        cd(condition_dir);
        
        if task == 'T'                                                         % If temporal task:
            ref_RFP_presentation = trial_fixation - trial_start_ref_RFP;            % Calculate presentation period of reference RFP (period it was visible).  
            reference_RFPs{k} = trial_start_ref_RFP * 1000;                           % Convert values to milliseconds and add to relevant cell arrays. 
            reference_RFP_pres{k} = (trial_fixation - reference_RFP) * 1000;
        end
       
        
        % Calculate presentation periods for inter-trial fixation, target RFP and response fixation.
        trial_fixation_presentation = target_RFP - trial_fixation;             
        target_RFP_presentation = response_fixation - target_RFP;
        
        trial_end = GetSecs();                                  % Record time at end of trial. 
%         cumulative_trial_end = trial_end - total_trial_start;   % Calculate cumulative trial durations at end of each trial.
        
        trial_end = trial_end - ind_trial_start;                % Calculate total length of trial - trial start time. 
        cumulative_trial_end = cumulative_trial_end + trial_end;
        
        response_fixation_presentation = trial_end - response_fixation;

        % Convert values to milliseconds and add to relevant cell arrays. 
        task_condition{k} = feature_task; trial_offsets{k} = trial_offset * 1000; trial_fixations{k} = trial_fixation * 1000; 
        trial_fixation_pres{k} = trial_fixation_presentation * 1000; target_RFPs{k} = target_RFP * 1000;
        target_RFP_pres{k} = target_RFP_presentation * 1000; response_fixations{k} = response_fixation * 1000;
        response_fixation_pres{k} = response_fixation_presentation * 1000; trial_ends{k} = trial_end * 1000;
        cumulative_trial_starts{k} = cumulative_trial_start * 1000; cumulative_trial_ends{k} = cumulative_trial_end * 1000;

        if pid < 10
            save(sprintf('%s_threshold_workspace_pp_0%d_scan_0%d',fixation, pid,scan_number));
        else
            save(sprintf('%s_threshold_workspace_pp_%d_scan_0%d',fixation, pid,scan_number));
        end
        cd(base_directory); % Move into function directory.  
    end
%     jheapcl;
end
% toc     % Calculate elapsed time for all trials. 

t=QuestMean(q); sd=QuestSd(q); % Calculate final threshold mean and standard deviation. 
t = 10^t;  % Convert mean threshold to non-log 10 units. 
fprintf('Final threshold estimate (mean+-sd) is %.2f +- %.2f\n',t,sd) % Print final threshold values to command window.

sca; PsychPortAudio('Close', pahandle);     % Close screen and psych-port audio.        
ListenChar(); ShowCursor;                   % Re-enable keypresses and mouse input. 

%% -------------------------- DATA OUTPUT ------------------------------ %%

cd(data_storage);   % Change directory to data storage directory (specified above). 

% general_threshold_dir = strcat(root_directory,'/data_output/experimental_thresholds');     % Index into a general 'experimental_thresholds' directory.
general_threshold_dir = strcat(root_directory,'\data_output\experimental_thresholds'); % Shuttle. 

cd(general_threshold_dir); 

if pid < 10
    threshold_dir = sprintf('0%d_experimental_thresholds', pid);  % Specify unique participant threshold directory.
else
    threshold_dir = sprintf('%d_experimental_thresholds', pid);  % Specify unique participant threshold directory.
end

if exist(threshold_dir, 'dir') == 0 % If the unique participant threshold directory does not exist;
    mkdir(threshold_dir); % Create participant threshold directory.
end

% Store all response/feature data in one transposed matrix. 
all_data = [trial_number; test_values; reference_value; movement_value; answers; key_pressed; reaction_time; num2cell(threshold)]';  

% Convert response/feature data to table with headings.
data_table = cell2table(all_data, 'VariableNames', {'Trial_Number' 'Target_Value' 'Reference_Value' 'Movement' 'Answer' 'Key_Pressed'...
                                                        'Reaction_Time' 'Threshold'});

plot(abs(threshold(practice_trials+1:end)));   % Plot participants threshold estimations across the experiment. 
ylim([0 max]); ylabel(sprintf('%s',feature_task)); xlim([1 (trials_desired-practice_trials)]); xlabel('Trials');  % Specify x and y labels and axis scales. 

if task == 'T'      %If temporal task:

    % Store all timing data in one transposed matrix (temporal stores reference RFP and inter-trial fixation times). 
    all_times = [trial_number; task_condition; cumulative_trial_starts; trial_offsets; reference_RFPs; reference_RFP_pres; trial_fixations; trial_fixation_pres; target_RFPs; target_RFP_pres; ...
        response_fixations; response_fixation_pres; trial_ends; cumulative_trial_ends]';
    
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames', {'Trial_Number' 'Condition' 'C_Trial_Start' 'Trial_Offset' 'Trial_Start' 'Reference_RFP_Presentation'...
                                                        'Trial_Fixation' 'Trial_Fixation_Presentation' 'Target_RFP' 'Target_RFP_Presentation' ...
                                                        'Response_Fixation' 'Response_Fixation_Presentation' 'Trial_End' 'C_Trial_End'});   
                                                    
elseif task == 'S' % If spatial task:
     
    % Store all timing data in one transposed matrix (spatial does not store reference RFP or inter-trial fixation times). 
    all_times = [trial_number; task_condition; cumulative_trial_starts; trial_offsets; target_RFPs; target_RFP_pres; response_fixations; response_fixation_pres; trial_ends; cumulative_trial_ends]';
    
    % Convert timing data to table with headings. 
    time_table = cell2table(all_times, 'VariableNames',{'Trial_Number' 'Condition' 'C_Trial_Start' 'Trial_Offset' 'Trial_Start' 'RFP_Presentation' 'Response_Fixation'...
        'Response_Fixation_Presentation' 'Trial_End' 'C_Trial_End'});
end

cd(data_storage); cd(pid_directory); cd(scan_dir); cd(condition_dir);                                                                        % Index into individual participant directory. 

if pid < 10
    writetable(data_table, sprintf('%s_response_data_%s_scan_0%d_pp_0%d_%s.csv', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM'))); % Save response data as .csv file. 
    writetable(time_table, sprintf('%s_timings_%s_scan_0%d_pp_0%d_%s.csv', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM'))); % Save timing data as .csv. 
    saveas(gcf,sprintf('%s_thresholds_%s_scan_0%d_pp_0%d_%s.jpg', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM')));                        % Save figure.
else
    writetable(data_table, sprintf('%s_response_data_%s_scan_0%d_pp_%d_%s.csv', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM'))); % Save response data as .csv file. 
    writetable(time_table, sprintf('%s_timings_%s_scan_0%d_pp_%d_%s.csv', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM'))); % Save timing data as .csv. 
    saveas(gcf,sprintf('%s_thresholds_%s_scan_0%d_pp_%d_%s.jpg', lower(feature_task), condition_task, scan_number, pid, datestr(now, 'ddmmyyyy_HHMM')));                        % Save figure.  
end

cd(general_threshold_dir); cd(threshold_dir);   % Index into participant threshold directory.    

if isequal(feature_task, 'Orientation')     % If orientation task:
    orientationT = t; orientationSD = sd;       % Rename final mean and standard deviation for data storage. 
    
    if pid < 10
        save(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat', condition_task, scan_number, pid), 'orientationT', 'orientationSD');   % Save threshold and SD to .mat file. 
    else
        save(sprintf('thresholds_%s_scan_0%d_pp_%d.mat', condition_task, scan_number, pid), 'orientationT', 'orientationSD');   % Save threshold and SD to .mat file. 
    end
    
elseif isequal(feature_task, 'Contrast')   
    contrastT = t; contrastSD = sd;
    
    if pid < 10
        save(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat', condition_task, scan_number, pid), 'contrastT', 'contrastSD', '-append');  % Append threshold and SD to file.
    else
        save(sprintf('thresholds_%s_scan_0%d_pp_%d.mat', condition_task, scan_number, pid), 'contrastT', 'contrastSD', '-append');  % Append threshold and SD to file.
    end
    
elseif isequal(feature_task, 'RAM')  
    RAMT = t; RAMSD = sd;
    
    if pid < 10
        save(sprintf('thresholds_%s_scan_0%d_pp_0%d.mat', condition_task, scan_number, pid), 'RAMT', 'RAMSD', '-append');
    else
        save(sprintf('thresholds_%s_scan_0%d_pp_%d.mat', condition_task, scan_number, pid), 'RAMT', 'RAMSD', '-append');
    end
end

cd(data_storage); cd(pid_directory); cd(scan_dir); cd(condition_dir); 

if pid < 10
    save(sprintf('%s_threshold_workspace_pp_0%d_scan_0%d',fixation, pid,scan_number));
    save(sprintf('%s_threshold_workspace_pp_%d_scan_0%d_%s',fixation, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));
else
    save(sprintf('%s_threshold_workspace_pp_%d_scan_0%d',fixation, pid,scan_number));
    save(sprintf('%s_threshold_workspace_pp_%d_scan_0%d_%s',fixation, pid, scan_number, datestr(now, 'ddmmyyyy_HHMM')));
end


cd(start_dir);
%% --------------------------------------------------------------------- %%