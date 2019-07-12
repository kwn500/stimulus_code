function run_colour_fmri(pid,r_number, run_number)

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

% --------------------- RUN COLOUR FMRI EXPERIMENT ---------------------- %

% fMRI block & event experiment.
% blocks- attend to orientation, contrast, shape, passive.
% blocks- attend to red/green, blue/yellow, luminance.

% events- orientation, contrast, shape x rg, by, luminance.
% also classifies event by detection (hit, miss, false alarm & correct
% rejection)- participant dependent.

% total run time- 5:06 (102 TRs).
% 9s dummy volumes- (3TRs).
% 4 feature blocks, repeated once per colour = 12 blocks (15s/5TRs each).
% 9s between each block (fixation and upcoming block cue).

% KWN 06/04/2018

% -------------------- VITAL EXPERIMENT INFORMATION --------------------- %

exp.trials_desired = 10; % number of trials.
exp.block_repeats = 3;   % number of block repeats (one run).
exp.colour_repeats = 1;  % number of colour repeats (one run).
exp.blocks_desired = [1 2 3 4]; % number of feature blocks (one set).
exp.prob_nochange = 0.8; % probability of no change in stim feature (80%).

exp.stim_duration = 0.2; % experimental timings (seconds).
exp.response_duration = 0.8;
exp.interblock_fix_duration = 7.5;
exp.interblock_cue_duration = 1.5;
exp.trial_length = 1.5;
exp.dummy_duration = 9;

% calculate number of trials required to filly dummy period.
exp.dummy_trials_desired = exp.dummy_duration / exp.trial_length;

exp.feature_cues = {'ORIENTATION', 'CONTRAST', 'SHAPE', 'PASSIVE'}; % specify block conditions.
exp.colour_cues = {'rg', 'by', 'luminance'}; % specify colour conditions.

phi = 0; C= 1; A= 0.5; ref = [phi C A]; % default values for RFP (orientation, contrast, shape).

% ---------------------------- GENERAL SETUP ---------------------------- %

storage.start = pwd; % store location of original directory.
storage.home = 'C:\Users\Ryan Maloney\Documents\Kirstie\colour'; % specify 'base' directory.
storage.psych_data = strcat(storage.home, '\psychophysics_data'); % specify psychophysics data storage directory.
storage.fmri_data = strcat(storage.home, '\fmri_data'); % specify fmri data storage directory.

addpath(genpath(strcat(storage.home, '\code\psychophysics_code\functions\'))); % add function directory.

% specify screen-specific calibration file.
dpy.calib_file = strcat(storage.home, '\calibration_files\','MRI_PROPixx_GE_spectra_gamma.mat');
dpy.backRGB.dir = [1 1 1]'; dpy.backRGB.scale = 0.5; % specify RGB background (mid-gray).

% ---------- LOAD STOCKMAN SHARPE 10-DEGREE CONE FUNDAMENTALS ----------- %

dpy.sensors_raw=load('T_cones_ss10'); % load stockman sharpe 10deg cone fundamentals

% resample cone fundamentals to dimensions needed by getLMS2RGB_gamma and assciated functions.
dpy.sensors=interp1(390:1:830,dpy.sensors_raw.T_cones_ss10',380:1:740);
dpy.sensors(isnan(dpy.sensors))=0; % replace NaNs (from conversion from 390-380) with zeros.

% ------------------- COLLECT PARTICIPANT INFORMATION ------------------- %
wpid = sprintf('%d', pid);
% add leading zeros to participant numbers < 10 (for uniform file naming).
if pid < 10, wpid = sprintf('0%d', pid); end

cd(storage.fmri_data); % index into data storage directory.

% specify, create & index into a participant-specific directory.
storage.pp = sprintf('pp_%s_%s', wpid, r_number); [~, ~] = mkdir(storage.pp); cd(storage.pp);

storage.run = sprintf('run_0%d', run_number); % specify a run-specific directory.

% if exist(storage.run, 'dir') ~= 0 % if this run already exists, exit the experiment.
%     fprintf('this run already exists!\n'); cd(storage.home); return
% else
    [~,~] = mkdir(storage.run); % otherwise, make the run-specific directory.
% end

cd(storage.psych_data); cd(storage.pp);

% -------------------- LOAD PARTICIPANT THRESHOLDS ---------------------- %

% load in participant thetas from isoluminance task.
load(sprintf('rg_isoluminance_theta_average_pp_%s.mat', wpid), 'mean_rg_theta');
load(sprintf('by_isoluminance_theta_average_pp_%s.mat', wpid), 'mean_by_theta');

% % load in 75% correct minimum contrast detection thresholds (rg, by & luminance).
% load(sprintf('lum_averaged_contrast_detection_threshold_pp_%s.mat', wpid), 'lum_min_contrast_mean');
% load(sprintf('rg_averaged_contrast_detection_threshold_pp_%s.mat', wpid), 'rg_min_contrast_mean');
% load(sprintf('by_averaged_contrast_detection_threshold_pp_%s.mat', wpid), 'by_min_contrast_mean');

% load in the participant colour and feature-specific averaged thresholds.
load(sprintf('by_averaged_feature_thresholds_pp_%s.mat', wpid), 'byo_tmean', 'byc_tmean', 'bys_tmean');
load(sprintf('rg_averaged_feature_thresholds_pp_%s.mat', wpid), 'rgo_tmean', 'rgc_tmean', 'rgs_tmean');
load(sprintf('lum_averaged_feature_thresholds_pp_%s.mat', wpid), 'lumo_tmean', 'lumc_tmean', 'lums_tmean');

% double these particpant feature-staircase thresholds.
thresholds.byo = byo_tmean*2; thresholds.byc = byc_tmean*2; thresholds.bys = bys_tmean*2;
thresholds.rgo = rgo_tmean*2; thresholds.rgc = rgc_tmean*2; thresholds.rgs = rgs_tmean*2;
thresholds.lumo = lumo_tmean*2; thresholds.lumc = lumc_tmean*2; thresholds.lums = lums_tmean*2;

% combine these participant thresholds into a singular matrix.
feature_thresholds = [thresholds.rgo, thresholds.rgc, thresholds.rgs; thresholds.byo, thresholds.byc,...
    thresholds.bys; thresholds.lumo, thresholds.lumc, thresholds.lums];

% --------------- CREATE RANDOMISED BLOCK & TRIAL ORDERS ---------------- %

% create a pseudo-randomised order of feature and colour block conditions.
[exp.block_order] = randomise_block_order(exp.feature_cues, exp.colour_cues, exp.block_repeats);

% create a randomised trial order within these blocks (in line with colour-specific orientation, contrast & shape thresholds).
[exp.trial_sequence, exp.dummy_trial_sequence] = randomise_trial_order(exp.trials_desired, exp.dummy_trials_desired,...
    exp.block_order, exp.feature_cues, exp.colour_cues, exp.prob_nochange, feature_thresholds, phi, C, A);

% ------------------------ SCREEN/KEYBOARD SETUP ------------------------ %

PsychDefaultSetup(2);  KbName('UnifyKeyNames'); % match the keynames across operating systems and setup psychtoolbox.
change1 = '1!'; change2 = '2@'; no_change1 = '3#'; no_change2 = '4$'; % set keynames for participant responses.
% space = KbName('space'); escape = KbName('escape');  trigger = KbName('5%');

screens = Screen('Screens');  screenNumber = max(screens);  % draw to max number screen.
white = WhiteIndex(screenNumber); % define black.

[window] = PsychImaging('OpenWindow', screenNumber, (128/255)); % open gray-screen window.
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); % alpha-blending for anti-aliased lines.

% pre-generate all possible RFP stimuli using participant's thresholds.
[dpy, rg, by, lum] = pregenerate_stimuli(exp.colour_cues, dpy, (0.05*2), (0.04*2),...
    (0.15*2), mean_rg_theta, mean_by_theta, feature_thresholds, phi, A, window);

[stim.fixations] = pregenerate_fixations(window, storage.home); % pre-generate all required stimulus fixations and cues.

topPriorityLevel = MaxPriority(window); Priority(topPriorityLevel); % set script to maximum processing priority.
Screen('LoadNormalizedGammaTable', window, dpy.gamma.inverse); % run gamma correction.

SetMouse(5000,0); HideCursor(window); % make sure mouse is hidden.

% from andre: establish keyboard queue variable (needed to flush events later).
% manually allow keys (11:14 - 1,2,3,4 on linux system, 1!, 2@, 3# and 4$ (49-52 on a win system));
allowedKeys = zeros(1,256);
targetKeys = [KbName({'1!','2@','3#','4$'}) 11:14 49:52];
quitKey = KbName('q');
allowedKeys([targetKeys quitKey]) = 1;

KbQueueCreate([],allowedKeys); % create keyboard queue (but don't start it yet!)

ListenChar(-1); % disable keypresses to command window. 

% ----------------------- TIMING INFORMATION ---------------------------- %

exp.ifi = Screen('GetFlipInterval', window); % get screen refresh rate.

% get number of frames corresponding to desired stimulus durations.
stimDurationTimeFrames = round(exp.stim_duration / exp.ifi);
responseDurationTimeFrames = round(exp.response_duration / exp.ifi);
interblockDurationTimeFrames = round(exp.interblock_fix_duration / exp.ifi);
cueDurationTimeFrames = round(exp.interblock_cue_duration / exp.ifi);

waitframes = 1; % number of frames to wait before re-drawing

% ------------------------ INTRODUCTION SCREEN -------------------------- %

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial'); % set text parameters.

% specify instructions for selected condition.
line1 = 'Waiting for trigger\n'; line2 = '1 or 2 = change\n'; line3 = '3 or 4 = no change';

DrawFormattedText(window, [line1 line2 line3], 'center', 'center', white);
Screen('Flip',window); KbWait; % draw text to screen and wait for keypress.

keyPress = 'not_five';
while keyPress ~= '5'
    [~,key] = KbWait; keyPress = KbName(key);
    if iscell(keyPress)
        keyPress = keyPress{1};
    end
end

KbQueueStart; % start recording keypresses

% -------------------------- RUN fMRI EXPERIMENT------------------------- %

trial = 0; counters.dummy = 1; % establish counter variables.

timing.total_starts{1} = 9; % set the start of the timing data to the end of the dummy volumes.

tic
for block_number = 1:(length(exp.block_order)) % for each block in turn:
    
    % --------------------------- DUMMY VOLUMES ----------------------------- %
    
    if counters.dummy == 1 % if the dummy volumes have not yet been presented:
        
        for dummy_trial = 1:exp.dummy_trials_desired % for each of the dummy trials specified:
            
            timing.trial_starts.dummy{dummy_trial} = GetSecs; % record time at the start of the trial.
            
            current_stim.all = lum.stim; % specify luminance stimuli.
            
            Screen('DrawTexture', window, current_stim.all.no_change_colour_tex); % draw reference RFP texture.
            Screen('DrawTexture', window, stim.fixations.passive_fix); % draw condition-specific fixation letter.
            
            % flip these stimuli to the screen and record the time.
            [timing.ref_starts.dummy{dummy_trial}] = deal(Screen('Flip', window));
            for frame = 1:stimDurationTimeFrames % for the number of frames required to match the stimulus duration:
                
                Screen('DrawTexture', window, current_stim.all.no_change_colour_tex); % draw reference RFP texture.
                Screen('DrawTexture', window, stim.fixations.passive_fix); % draw condition-specific fixation letter.
                
                % record timing at end of stimulus presentation period.
                timing.ref_ends.dummy{dummy_trial} = Screen('Flip', window, timing.ref_starts.dummy{dummy_trial} + (waitframes-0.5)*exp.ifi);
            end
            
            Screen('DrawTexture', window, stim.fixations.passive_fix); % draw the condition-specific fixation cross.
            counter = 0; % establish counter variable.
            
            % repeat the same process of drawing stimuli and flipping to the screen each frame as above.
            timing.interstimulus_starts.dummy{dummy_trial} = Screen('Flip', window);
            for frame = 1:stimDurationTimeFrames
                Screen('DrawTexture', window, stim.fixations.passive_fix);
                
                while counter < 1 % to ensure this loop is only completed once:
                    
                    % combine the trial-specific orientation, contrast and shape change variables.
                    feature_changes = [exp.dummy_trial_sequence(3,dummy_trial), exp.dummy_trial_sequence(4,dummy_trial),...
                        exp.dummy_trial_sequence(5,dummy_trial)];
                    
                    % create a matrix of logicals with 1- no change, 0-change in orientation, contrast & shape respectively.
                    changes = feature_changes == ref;
                    
                    % select and draw the relevant pre-generated stimulus based on the trial-specific feature changes.
                    [current_stim.target_draw] = get_target_RFP(current_stim.all, changes);
                    [current_stim.target_tex] = Screen('MakeTexture', window, current_stim.target_draw);
                    
                    counter = counter + 1; % increment the counter variable to prevent repetitions of the loop.
                end
                timing.interstimulus_ends.dummy{dummy_trial} = Screen('Flip', window, timing.interstimulus_starts.dummy{dummy_trial} + (waitframes-0.5)*exp.ifi);
            end
            
            Screen('DrawTexture', window, current_stim.target_tex); % draw target RFP.
            Screen('DrawTexture', window, stim.fixations.passive_fix); % draw condition-specific fixation letter.
            
            % repeat the same timing process as before for the target stimulus.
            timing.target_starts.dummy{dummy_trial} = Screen('Flip',window);
            for frame = 1:stimDurationTimeFrames
                Screen('DrawTexture', window, current_stim.target_tex);
                Screen('DrawTexture', window, stim.fixations.passive_fix);
                
                timing.target_ends.dummy{dummy_trial} = Screen('Flip', window, timing.target_starts.dummy{dummy_trial} + (waitframes-0.5)*exp.ifi);
            end
            
            Screen('DrawTexture', window, stim.fixations.passive_fix); % draw the response fixation cross.
            counter = 0; % establish a counter variable.
            
            % repeat the same process for the response duration period.
            timing.response_starts.dummy{dummy_trial} = Screen('Flip',window);
            for frame = 1:responseDurationTimeFrames
                Screen('DrawTexture', window, stim.fixations.passive_fix);
                
                while counter < 1
                    % calculate the relative timings for each of the stimuli by subtracting the start time.
                    timing.ref_durations.dummy{dummy_trial} = timing.ref_ends.dummy{dummy_trial} - timing.ref_starts.dummy{dummy_trial};
                    timing.ref_starts.dummy{dummy_trial} = timing.ref_starts.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    timing.ref_ends.dummy{dummy_trial} = timing.ref_ends.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    
                    timing.interstimulus_durations.dummy{dummy_trial} = timing.interstimulus_ends.dummy{dummy_trial} - timing.interstimulus_starts.dummy{dummy_trial};
                    timing.interstimulus_starts.dummy{dummy_trial} = timing.interstimulus_starts.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    timing.interstimulus_ends.dummy{dummy_trial} = timing.interstimulus_ends.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    
                    timing.target_durations.dummy{dummy_trial} = timing.target_ends.dummy{dummy_trial} - timing.target_starts.dummy{dummy_trial};
                    timing.target_starts.dummy{dummy_trial} = timing.target_starts.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    timing.target_ends.dummy{dummy_trial} = timing.target_ends.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    
                    counter = counter + 1; % increment counter variable.
                end
                timing.response_ends.dummy{dummy_trial} = Screen('Flip', window, timing.response_starts.dummy{dummy_trial} + (waitframes-0.5)*exp.ifi);
            end
            % calculate the duration of the response period.
            timing.response_durations.dummy{dummy_trial} = timing.response_ends.dummy{dummy_trial} - timing.response_starts.dummy{dummy_trial};
            
            % calculate the amount of time remaining until the desired trial length is reached.
            time_remaining = exp.trial_length - (timing.response_ends.dummy{dummy_trial}- timing.trial_starts.dummy{dummy_trial});
            
            timedout = false; % establish timing variable.
            while ~timedout % while the total trial duration has not yet been met:
                
                % if the length of this trial is exceeds the specified trial duration.
                if((GetSecs - timing.response_ends.dummy{dummy_trial}) >= time_remaining)
                    timedout = true; % exit the timing loop.
                    
                    % record the time at the end of the trial.
                    timing.trial_ends.dummy{dummy_trial} = GetSecs - timing.trial_starts.dummy{dummy_trial};
                    
                    % calculate relative timings for the response data.
                    timing.response_starts.dummy{dummy_trial} = timing.response_starts.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                    timing.response_ends.dummy{dummy_trial} = timing.response_ends.dummy{dummy_trial} - timing.trial_starts.dummy{dummy_trial};
                end
            end
            
            if dummy_trial == exp.dummy_trials_desired % if the dummy trial is equal to the total number of dummy trials desired:
                
                counter = 0; % establish counter variable.
                Screen('DrawTexture', window, stim.fixations.ib_fix); % draw the black inter-block fixation cross.
                
                timing.interblock_starts.dummy = Screen('Flip', window); % flip to screen and record time.
                for frame = 1:interblockDurationTimeFrames % repeat same timing process as before.
                    
                    while counter < 1 % to ensure we only complete this loop once:
                        
                        % get the upcoming block name (orientation, contrast or shape).
                        feature_task = exp.feature_cues{exp.trial_sequence(1,1)};
                        current_stim.fix = stim.fixations.(sprintf('%s_cue_fix', lower(feature_task))); % get the corresponding block cue.
                        
                        counter = counter + 1; % increment counter variable.
                    end
                    Screen('DrawTexture', window, stim.fixations.ib_fix);
                    timing.interblock_ends.dummy = Screen('Flip', window, timing.interblock_starts.dummy + (waitframes-0.5)*exp.ifi);
                end
                
%                 Screen('DrawTexture', window, current_stim.fix); % draw the upcoming block cue.
                Screen('DrawTexture', window, stim.fixations.ib_fix);
                
                timing.cue_starts.dummy = Screen('Flip', window); % present for the specified cue-duration time.
                for frame = 1:cueDurationTimeFrames
%                     Screen('DrawTexture', window, current_stim.fix);
                    Screen('DrawTexture', window, stim.fixations.ib_fix);
                    timing.cue_ends.dummy = Screen('Flip', window, timing.cue_starts.dummy + (waitframes-0.5)*exp.ifi);
                end
            end
        end
        clear current_stim % clear the current_stim information to ensure new stimuli are created each trial.
        counters.dummy = counters.dummy + 1; % increment the counter variable to esure only one set of dummy volumes.
    end
    
    % ----------------- RUN A COLOUR X FEATURE-SPECIFIC BLOCK --------------- %
    
    for block_trial = 1:exp.trials_desired % for each trial:
        
        % increment the trial counter variable and record the trial start time.
        trial = trial + 1; trial_data.counter{trial} = trial;
        timing.trial_starts.exp{trial} = GetSecs;
        
        if block_trial == 1 % if this is the first trial of the block:
            % get the block-specific feature and colour conditions.
            if block_number == 1
                feature_task = exp.feature_cues{exp.trial_sequence(1,block_trial)};
                colour_condition = exp.colour_cues{exp.trial_sequence(2,block_trial)};
            elseif block_number ~= 1
                feature_task = exp.feature_cues{exp.trial_sequence(1,((block_number-1)*10)+block_trial)};
                colour_condition = exp.colour_cues{exp.trial_sequence(2,((block_number-1)*10)+block_trial)};
            end
        end
        
        % select the relevant colour and feature-data and store some trial-specific data.
        [current_stim, trial_data] = get_trial_specific_data(colour_condition, thresholds,...
            rg, by, lum, trial, phi, A, feature_task, stim, trial_data);
        
        % draw the reference stimulus and condition-specific fixation.
        Screen('DrawTexture', window, current_stim.ref); Screen('DrawTexture', window, current_stim.fix);
        counter = 0; % establish counter variable.
        
        [timing.ref_starts.exp{trial}] = deal(Screen('Flip', window)); % repeat the same timing procedure as previously.
        for frame = 1:stimDurationTimeFrames
            
            while counter < 1 % to ensure this loop is repeated only once:
                
                trial_data.features{trial} = lower(feature_task); % store the trial feature condition.
                trial_data.colours{trial} = lower(colour_condition); % store the trial colour condition.
                
                % get the trial-specific target values for orientation, contrast and shape.
                if block_number == 1
                    inp=[exp.trial_sequence(3,block_trial) exp.trial_sequence(4,block_trial) exp.trial_sequence(5,block_trial)];
                elseif block_number ~= 1
                    inp=[exp.trial_sequence(3,(block_trial+((block_number-1)*10))) exp.trial_sequence(4,(block_trial+((block_number-1)*10))) exp.trial_sequence(5,(block_trial+((block_number-1)*10)))];
                end
                
                % calculate the relative timings for the dummy interblock and cue stimuli.
                timing.interblock_durations.dummy = timing.interblock_ends.dummy - timing.interblock_starts.dummy;
                timing.interblock_starts.dummy = timing.interblock_starts.dummy - timing.trial_starts.dummy{dummy_trial};
                timing.interblock_ends.dummy = timing.interblock_ends.dummy - timing.trial_starts.dummy{dummy_trial};
                
                timing.cue_durations.dummy = timing.cue_ends.dummy - timing.cue_starts.dummy;
                timing.cue_starts.dummy = timing.cue_starts.dummy - timing.trial_starts.dummy{dummy_trial};
                timing.cue_ends.dummy = timing.cue_ends.dummy - timing.trial_starts.dummy{dummy_trial};
                
                % set default response parameters.
                trial_data.RT{trial} = 'no_response';           trial_data.keyLabel{trial} = 'N/A';
                trial_data.type_changes{trial} = 'no_change';   trial_data.response_needed{trial} = 'N';
                trial_data.marked_response{trial} = 'no_response';     trial_data.keyCode{trial} = [];
                [trial_data.response{trial}, trial_data.conjunction{trial}] = deal('no');
                [trial_data.orientation.track{trial}, trial_data.contrast.track{trial}, trial_data.shape.track{trial},...
                    trial_data.num_changes{trial}] = deal(0);
                
                counter = counter + 1; % increment counter variable.
            end
            Screen('DrawTexture', window, current_stim.ref);
            Screen('DrawTexture', window, current_stim.fix);
            timing.ref_ends.exp{trial} = Screen('Flip', window, timing.ref_starts.exp{trial} + (waitframes-0.5)*exp.ifi);
        end
        timing.ref_durations.exp{trial} = timing.ref_ends.exp{trial} - timing.ref_starts.exp{trial};
        
        Screen('DrawTexture', window, current_stim.fix); % draw the condition-specific fixation cross. 
        counter = 0; % establish counter variable. 
        
        timing.interstimulus_starts.exp{trial} = Screen('Flip', window);
        for frame = 1:stimDurationTimeFrames
            while counter < 1
                changes=inp-ref~=0; % get a logical of nonzero locations (where input and reference were unequal).                
                num_change=sum(changes); % sum these logicals together.   
                
                trial_data.num_changes{trial} = num_change;

                if num_change > 0 % if some change from the reference occured:
                    fmt=[repmat('%s_&_',1,num_change-1) '%s']; % create a dynamic format string. 
                    trial_data.type_changes{trial} = sprintf(fmt, exp.feature_cues{changes}); % store the feature-changes as a string. 
                    trial_data.type_changes{trial} = lower(trial_data.type_changes{trial});
                    trial_data.response_needed{trial} = 'Y'; % store a response as needed. 
                end
                
                % if changes in multiple features occured, store this trial as a conjunction. 
                if num_change > 1, trial_data.conjunction{trial} = 'yes'; end 
                
%                 changes = inp ~= ref; % get a logical of nonzero locations (where input and reference were unequal). 
                
                % select and draw the relevant pre-generated stimulus based on the trial-specific feature changes.
                [current_stim.target_draw] = get_target_RFP(current_stim.all, changes); % 
                [current_stim.target_tex] = Screen('MakeTexture', window, current_stim.target_draw);
                [trial_data] = get_target_stim_values(colour_condition, trial, thresholds, phi, A, changes, lum, rg, by, trial_data);
                
                counter = counter + 1; % increment counter variable. 
            end
            Screen('DrawTexture', window, current_stim.fix);
            timing.interstimulus_ends.exp{trial} = Screen('Flip', window, timing.interstimulus_starts.exp{trial} + (waitframes-0.5)*exp.ifi);
        end
        timing.interstimulus_durations.exp{trial} = timing.interstimulus_ends.exp{trial} - timing.interstimulus_starts.exp{trial};
        
        Screen('DrawTexture', window, current_stim.target_tex); % draw target RFP.
        Screen('DrawTexture', window, current_stim.fix); % draw condition-specific fixation cross. 
        
        timing.target_starts.exp{trial} = Screen('Flip', window); % present the target RFP and fixation. 
        for frame = 1:stimDurationTimeFrames
            Screen('DrawTexture', window, current_stim.target_tex);
            Screen('DrawTexture', window, current_stim.fix);
            timing.target_ends.exp{trial} = Screen('Flip', window, timing.target_starts.exp{trial} + (waitframes-0.5)*exp.ifi);
        end
        timing.target_durations.exp{trial} = timing.target_ends.exp{trial} - timing.target_starts.exp{trial};
        
        Screen('DrawTexture', window, current_stim.fix); % draw the response condition-specific fixation. 
        
        timedout = false; counter = 0; % establish counter variables.  
        
        [keyTime, timing.response_starts.exp{trial}] = deal(Screen('Flip',window));
        
        % calculate time remaining until desired trial length is reached. 
        time_remaining = exp.trial_length - (keyTime- timing.trial_starts.exp{trial}); 
        while ~timedout
            
            % if the total trial length time has been reached: 
            if((GetSecs - timing.response_starts.exp{trial}) >= time_remaining)
                timing.response_ends.exp{trial} = GetSecs; % store the trial end timing and exit the loop.
                timedout = true;
            end
            
            [pressed,firstPress,~,~,~] = KbQueueCheck; % otherwise, wait for a keypress. 
            if pressed % if the participant makes a keypress:
                if any(firstPress(targetKeys)) % check if they have pressed one of the response keys. 
                    trial_data.keyLabel{trial} = KbName(firstPress); % if so, store the keyname of this keypress. 
                    
                    % if they have pressed two keys at once, store the keyname of the first key. 
                    if iscell(trial_data.keyLabel{trial}), trial_data.keyLabel{trial} = trial_data.keyLabel{trial}{1}; end
                    
                    trial_data.RT{trial} = GetSecs - timing.response_starts.exp{trial}; % calculate and store the keypress reaction time. 
                    
                    while counter < 1 % ensure the loop is only completed once
                        
                        % store the participant response ('change'/'no_change') based on their keypress. 
                        if isequal(trial_data.keyLabel{trial}, change1) == 1 || isequal(trial_data.keyLabel{trial}, change2) == 1          
                            trial_data.response{trial} = 'change';
                        elseif isequal(trial_data.keyLabel{trial}, no_change1) == 1 || isequal(trial_data.keyLabel{trial}, no_change2) == 1    
                            trial_data.response{trial} = 'no_change';  
                        end
                        
                        % classify participant response as hit, miss, correct rejection or false alarm.
                        [marked_response] = d_prime_classification(trial_data.response{trial},...
                            trial_data.type_changes{trial}, exp.trial_sequence(1,block_trial));
                        trial_data.marked_response{trial} = marked_response;
                        trial_data.marked_response{trial}
                        counter = counter + 1; % increment counter variable. 
                    end
                end
            end
        end
        
        KbQueueFlush; % clear out keypress information so we don't fill up the buffer
        
        % store trial end timings and calculate relative timings. 
        timing.response_durations.exp{trial} = timing.response_ends.exp{trial} - timing.response_starts.exp{trial};
        timing.trial_ends.exp{trial} = timing.response_ends.exp{trial} - timing.trial_starts.exp{trial};
        
% ------------- PRESENT INTER-BLOCK INTERVAL AND BLOCK CUE -------------- %

        if block_trial == exp.trials_desired % if this is the last trial of a block    
                        
            Screen('DrawTexture', window, stim.fixations.ib_fix); % draw the black inter-block fixation. 
            counter = 0; % establish counter variable.
            
            timing.interblock_starts.exp{trial} = Screen('Flip', window);
            for frame = 1:interblockDurationTimeFrames -1
                Screen('DrawTexture', window, stim.fixations.ib_fix);
                
                while counter < 1
                    
                    % get the name of the upcoming condition cue and load the relevant cue texture. 
                    if block_number == length(exp.block_order)
                        current_stim.fix = stim.fixations.ib_fix;
                    elseif block_number < length(exp.block_order)
                        feature_task = exp.feature_cues{exp.trial_sequence(1,(10*(block_number))+1)};
                        stim_name = sprintf('%s_cue_fix', lower(feature_task));
                        current_stim.fix = stim.fixations.(stim_name);
                    end
                    counter = counter + 1; % increment counter variable. 
                end
                timing.interblock_ends.exp{trial} = Screen('Flip', window, timing.interblock_starts.exp{trial} + (waitframes-0.5)*exp.ifi);
            end
            timing.interblock_durations.exp{trial} = timing.interblock_ends.exp{trial} - timing.interblock_starts.exp{trial};
            
%             Screen('DrawTexture', window, current_stim.fix); % draw the upcoming cue. 
            Screen('DrawTexture', window, stim.fixations.ib_fix);
            
            timing.cue_starts.exp{trial} = Screen('Flip', window);
            for frame = 1:cueDurationTimeFrames -1
                
                Screen('DrawTexture', window, stim.fixations.ib_fix);
                timing.cue_ends.exp{trial} = Screen('Flip', window, timing.cue_starts.exp{trial} + (waitframes-0.5)*exp.ifi);
            end
            timing.cue_durations.exp{trial} = timing.cue_ends.exp{trial} - timing.cue_starts.exp{trial};
            timing.trial_ends.exp{trial} = timing.cue_ends.exp{trial} - timing.trial_starts.exp{trial};
        end
        
        if trial == 1 % if this is the first trial:
            % set the time of the first trial end as the trial duration the length of the dummy volumes (9s).
            timing.total_ends{trial} = 9 + round(timing.trial_ends.exp{trial},1);
        elseif trial ~= 1 % otherwise:
            
            % specify the cumulative trial start and trial end timings. 
            timing.total_starts{trial} = timing.total_starts{trial-1} + round(timing.trial_ends.exp{trial-1},1);
            timing.total_ends{trial} = timing.total_ends{trial-1} + round(timing.trial_ends.exp{trial-1},1);
        end
    end
end
toc
KbQueueRelease; sca; ShowCursor; % end the keyboard cue, clear the screen and show the mouse cursor. 
ListenChar; % renable keyboard presses.

% create log files of the trial timing datas for fsl, freesurfer and mrvista. 
timing = classify_eventfile_data(timing, trial_data, trial, wpid, run_number, storage.fmri_data, storage.pp, storage.run);

cd(storage.fmri_data); cd(storage.pp); cd(storage.run); % index into participant and run-specific storage directory. 

% create a cell array of participant response and timing data to output to .csv file. 
excel_output = [trial_data.counter; trial_data.features; trial_data.colours; trial_data.orientation.thresh; trial_data.contrast.thresh;...
    trial_data.shape.thresh; trial_data.orientation.track; trial_data.contrast.track; trial_data.shape.track; trial_data.num_changes;...
    trial_data.conjunction; trial_data.type_changes; trial_data.response_needed; trial_data.response; trial_data.keyLabel;...
    trial_data.marked_response; trial_data.RT; timing.total_starts; timing.ref_durations.exp; timing.ref_ends.final;...
    timing.interstimulus_durations.exp; timing.interstimulus_ends.final; timing.target_durations.exp; timing.target_ends.final;...
    timing.response_durations.exp; timing.response_ends.exp; timing.trial_ends.exp; timing.interblock_durations.exp;...
    timing.interblock_ends.final; timing.cue_durations.exp; timing.cue_ends.final; timing.total_ends]';

% convert this data to a table with column headings. 
data_table = cell2table(excel_output, 'VariableNames', {'trial_number', 'feature', 'colour', 'othresh', 'cthresh', 'sthresh', 'ochange', 'cchange',...
    'schange', 'num_changes', 'conjunction', 'type_changes', 'response_required', 'response', 'key_pressed', 'marked_response', 'RT',...
    'trial_start', 'ref_duration','ref_end', 'isi_duration', 'isi_end', 'target_duration', 'target_end', 'response_duration',....
    'response_end', 'trial_end', 'ib_duration', 'ib_end', 'cue_duration', 'cue_end', 'cum_end'});

% write this data to a .csv file.
writetable(data_table, sprintf('%s_%s_run_0%d_data.csv', wpid, r_number, run_number));

% clear all; clc; % clear the workspace and command window ready for the next run. 

% ----------------------------------------------------------------------- %