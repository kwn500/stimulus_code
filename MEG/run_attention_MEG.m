%% ACHROMATIC ATTENTION EXPERIMENT

% fMRI experiment.
% Displays blocks of trials in which RFP stimuli change in visual features
% at a specified probability. Participants are cued to the visual feature
% of attention.

% Records keypresses and output response and timing data to a .csv file.

% KWN 4/7/19

%%
close all;
clear;
KbName('UnifyKeyNames');
addpath('/home/k/kwn500/Desktop/ppdev-mex-master');

%% Collect participant information
response = inputdlg({'Enter Participant Number:', 'Enter Scan Number:'},...
    'Participant Information', 1);

%% Specify scan parameters
setup.ppd = 46.4; % specify pixels/degree (correct for MRI ProPixx).

%-----TRIALS
setup.trials = 10; % number of trials per block.
setup.block_repeats = 4; % number of block repeats (per condition).
setup.probability = 0.8; % probability of feature change.
setup.block_labels = {'ORIENTATION', 'CONTRAST', 'SHAPE', 'PASSIVE'};
setup.fixation_labels = {'O', 'C', 'S', 'P'};

% pseudo-randomise block order;
setup.block_order = [];
for thisBlock = 1:setup.block_repeats
    setup.block_order = [setup.block_order, randperm(length(setup.block_labels))];
end

%-----TIMINGS
setup.baseline_time = 7.5; % inter-block interval duration (s).
setup.cue_time = 1.5; % block attention cue duration(s).
setup.stim_time = 0.2; % stimulus duration (s).
setup.ISI_time = 0.2; % inter-stimulus duration (s).
setup.response_time = 0.8; % response duration (s).
setup.trial_time = 1.5; % total trial duration (s).

%-----KEYBOARD
setup.keys.space = KbName('space');
setup.keys.escape = KbName('escape');
setup.keys.trigger = KbName('5%');
setup.keys.change1 = '1!';
setup.keys.change2 = '2@';
setup.keys.nochange1 = '3#';
setup.keys.nochange2 = '4$';

%andre's additional keypress specification for use in MRI.
setup.keys.allowedKeys = zeros(1,256);
setup.keys.targetKeys = [KbName({'1','2','3','4'}), 11:14, 49:52];
setup.keys.quitKey = KbName('q');
setup.keys.allowedKeys([setup.keys.targetKeys setup.keys.quitKey]) = 1;

KbQueueCreate([], setup.keys.allowedKeys); % create keyboard queue.

%-----RFP PARAMETERS
% parameters for reference RFP.
setup.rfp.C = 0.5; % contrast
setup.rfp.A = 0.5; % shape.
setup.rfp.phi = 0; % orientation.

%-----THRESHOLD VALUES
% usually we load in participant- and feature-specific threshold values
% from prior psychophysics testing, but for now, we'll set some default
% values.
setup.othresh = 0.1;
setup.cthresh = 0.1;
setup.sthresh = 0.05;

setup.output_dir = ('\Users\Ryan Maloney\Documents\Kirstie\attn_experiment\');

ppdev_mex('Open', 1);
ppdev_mex('Write', 1, 0);

%% Screen setup
PsychDefaultSetup(2)

screens = Screen('Screens');
screenNumber = max(screens);

black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);

% open grey-screen window.
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, (128/255), [0 0 500 500]);

% hide mouse cursor and disable keypress output;
SetMouse(5000,0);
HideCursor(window);
ListenChar(0);

%% Introduction screen
Screen('TextSize', window, 25);
Screen('TextFont', window, 'Arial');

setup.introtext = 'Waiting for trigger\n\n1 or 2 = change\n3 or 4 = no change';
DrawFormattedText(window, setup.introtext, 'center', 'center', white);
Screen('Flip', window);

% wait for trigger (5)
keyPress = 'not_five';
while keyPress ~= '5'
    [~,key] = KbWait;
    keyPress = KbName(key);
    if iscell(keyPress)
        keyPress = keyPress{1};
    end
end

KbQueueStart; % start recording keypresses.

%% Trial structure

% pre-generate feature changes across all changes for more efficient code.
for thisBlock = 1:length(setup.block_order)
    for thisTrial = 1:setup.trials
        
        trials.o_seq(thisTrial) = setup.rfp.phi;
        trials.c_seq(thisTrial) = setup.rfp.C;
        trials.s_seq(thisTrial) = setup.rfp.A;
        
        % alter feature change on approximately 20% of trials.
        if rand() > setup.probability
            trials.o_seq(thisTrial) = setup.othresh + setup.rfp.phi;
        end
        if rand() > setup.probability
            trials.c_seq(thisTrial) = setup.cthresh + setup.rfp.C;
        end
        if rand() > setup.probability
            trials.s_seq(thisTrial) = setup.sthresh + setup.rfp.A;
        end
    end
    trials.block{thisBlock} = [trials.o_seq; trials.c_seq; trials.s_seq];
end

%% Run experiment
Priority(MaxPriority(window));

% pre-generate reference RFP for efficiency.
stim.ref_rfp = RFP_generator(setup.rfp.phi, setup.rfp.C, setup.rfp.A);
stim.ref_rfp_tex = Screen('MakeTexture', window, stim.ref_rfp);

total_trials = 1; % keep track of all trials.

for thisBlock = 1:length(setup.block_order)
    for thisTrial = 1:setup.trials
        
        % store important trial information.
        info.orientation(total_trials) = trials.block{thisBlock}(1,thisTrial);
        info.contrast(total_trials) = trials.block{thisBlock}(2,thisTrial);
        info.shape(total_trials) = trials.block{thisBlock}(3,thisTrial);
        info.block{total_trials} = setup.block_labels{setup.block_order(thisBlock)};
        
        %%-----ATTENTION CUE
        
        % if this is the first trial of the block, display the attention cue.
        if thisTrial == 1
            t0 = GetSecs; Screen('Flip', window);
            timing.exp_start = t0;
            timing.cue(total_trials) = GetSecs();
            
            time_diff = timing.cue(total_trials) - t0;
            
            waitTime = setup.cue_time - time_diff;
            endTime = waitTime + timing.cue(total_trials);
            
            while GetSecs < endTime
                Screen('DrawTexture', window, stim.ref_rfp_tex);
                DrawFormattedText(window,...
                    setup.fixation_labels{setup.block_order(thisBlock)},...
                    'center', 'center', white);
            end
            % send trigger matching block condition.
            ppdev_mex('Write', 1, setup.block_order(thisBlock));
        else
            Screen('DrawTexture', window, stim.ref_rfp_tex);
            DrawFormattedText(window,...
                setup.fixation_labels{setup.block_order(thisBlock)},...
                'center', 'center', white);
        end
        
        test.trial_start = GetSecs;
        
        %%-----REFERENCE RFP
        
        % display the reference radial frequency stimulus and generate the
        % ISI fixation.
        t0 = GetSecs; Screen('Flip', window);
        timing.reference(total_trials) = GetSecs();
         
        % turn off block trigger after single screen flip. 
        ppdev_mex('Write', 1, 0);
       
        time_diff = timing.reference(total_trials) - t0;
        
        waitTime = setup.stim_time - time_diff;
        endTime = waitTime + timing.reference(total_trials);
        
        while GetSecs < endTime
            DrawFormattedText(window,...
                setup.fixation_labels{setup.block_order(thisBlock)},...
                'center', 'center', white);
        end
        
        %%-----INTER-STIMULUS INTERVAL
        
        % display the ISI fixation and generate the target RFP.
        t0 = GetSecs; Screen('Flip', window);
        timing.ISI(total_trials) = GetSecs();
        
        time_diff = timing.ISI(total_trials) - t0;
        
        waitTime = setup.ISI_time - time_diff;
        endTime = waitTime + timing.ISI(total_trials);
        
        target_draw = 0;
        while GetSecs < endTime
            if target_draw == 0
                stim.target_rfp = RFP_generator(trials.block{thisBlock}(1,thisTrial),...
                    trials.block{thisBlock}(2,thisTrial),trials.block{thisBlock}(3,thisTrial));
                stim.target_rfp_tex = Screen('MakeTexture', window, stim.target_rfp);
                Screen('DrawTexture', window, stim.target_rfp_tex);
                
                DrawFormattedText(window,...
                    setup.fixation_labels{setup.block_order(thisBlock)},...
                    'center', 'center', white);
                
                target_draw = target_draw + 1;
                
                % extract event trigger code.
                [event_trigger] = get_event_trigger(trials.block, setup.rfp, thisTrial, thisBlock);
            end
        end
        
        %%-----TARGET RFP
        
        % display the target RFP and generate the response fixation.
        t0 = GetSecs; Screen('Flip', window);
        
        % display event-specific trigger code.
        ppdev_mex('Write', 1, event_trigger);
        
        timing.target(total_trials) = GetSecs();
        
        time_diff = timing.target(total_trials) - t0;
        
        waitTime = setup.stim_time - time_diff;
        endTime = waitTime + timing.target(total_trials);
        
        while GetSecs < endTime
            DrawFormattedText(window, 'X', 'center', 'center', black);
            
            info.num_changes(total_trials) = sum([info.orientation(total_trials) == setup.rfp.phi,...
                info.contrast(total_trials) == setup.rfp.C, info.shape(total_trials) == setup.rfp.A]);
        end
        
        %%-----RESPONSE
        
        % display the response fixation.
        t0 = GetSecs; Screen('Flip', window);
        
        % turn off event trigger code.
        ppdev_mex('Write', 1, 0);
        
        timing.response(total_trials) = GetSecs;
        
        endTime = timing.reference(total_trials) + setup.trial_time + 0.02;
        
        % specify some defaults in the case of no response.
        timing.RT(total_trials) = NaN;
        keyLabel = NaN;
        
        while GetSecs < endTime
            
            % if this is the end of a block, prep the inter-block interval
            % fixation stimulus.
            if thisTrial == setup.trials
                DrawFormattedText(window,'X','center','center', black);
            end
            
            % check to see if key has been pressed
            [pressed, firstPress,~,~,~] = KbQueueCheck;
            
            if pressed
                
                % check if they have pressed any of the response keys
                if any(firstPress(setup.keys.targetKeys))
                    keyLabel = KbName(firstPress);
                    
                    % store the reaction time of the response.
                    timing.RT(total_trials) = GetSecs - timing.response(total_trials);
                    
                    if iscell(keyLabel)
                        keyLabel = keyLabel{1};
                    end
                end
                
                % check if they have quit the programme
                if firstPress(setup.keys.quitKey)
                    finished = true;
                    currMessage = 'Finished';
                end
                
                KbQueueFlush; % clear keys to prevent filling the buffer.
            end
            
            % classify keypress into d' categories.
            if isequal(keyLabel, setup.keys.nochange1) || isequal(keyLabel, setup.keys.nochange2)
                info.keypress{total_trials} = 'No_Change';
                
                if ~isequal(info.block{total_trials}, 'PASSIVE')
                    if isequal(info.block{total_trials}, 'ORIENTATION') && info.orientation(total_trials) ~= setup.rfp.phi
                        info.response{total_trials} = 'Miss';
                    elseif isequal(info.block{total_trials}, 'CONTRAST') && info.contrast(total_trials) ~= setup.rfp.C
                        info.response{total_trials} = 'Miss';
                    elseif isequal(info.block{total_trials}, 'SHAPE') && info.shape(total_trials) ~= setup.rfp.A
                        info.response{total_trials} = 'Miss';
                    else
                        info.response{total_trials} = 'Correct_Rejection';
                    end
                else
                    info.response{total_trials} = 'N/A';
                end
            elseif isequal(keyLabel, setup.keys.change1) || isequal(keyLabel, setup.keys.change2)
                info.keypress{total_trials} = 'Change';
                
                if ~isequal(info.block{total_trials}, 'PASSIVE')
                    if isequal(info.block{total_trials}, 'ORIENTATION') && info.orientation(total_trials) ~= setup.rfp.phi
                        info.response{total_trials} = 'Hit';
                    elseif isequal(info.block{total_trials}, 'CONTRAST') && info.contrast(total_trials) ~= setup.rfp.C
                        info.response{total_trials} = 'Hit';
                    elseif isequal(info.block{total_trials}, 'SHAPE') && info.shape(total_trials) ~= setup.rfp.A
                        info.response{total_trials} = 'Hit';
                    else
                        info.response{total_trials} = 'False_Alarm';
                    end
                else
                    info.response{total_trials} = 'N/A';
                end
            else
                info.keypress{total_trials} = 'No_Response';
                info.response{total_trials} = 'No_Response';
            end
        end
        
        timing.exp_end = GetSecs;
        
        %%-----INTER-BLOCK BASELINE
        
        % if this is the end of a block, display the inter-block interval
        % fixation and prep the attention cue.
        if thisTrial == setup.trials
            t0 = GetSecs; Screen('Flip', window);
            timing.baseline(total_trials) = GetSecs();
            
            time_diff = timing.baseline(total_trials) - t0;
            
            waitTime = setup.baseline_time - time_diff;
            endTime = waitTime + timing.baseline(total_trials);
            
            while GetSecs < endTime
                
                if thisBlock < length(setup.block_order)
                    DrawFormattedText(window,setup.block_labels{setup.block_order(thisBlock+1)},...
                        'center', 'center', white');
                end
            end
        end
        
        total_trials = total_trials + 1;
    end
end

Screen('CloseAll')
ppdev_mex('Close', 1);
KbQueueRelease;
ShowCursor;

%% Format timing data
edit_timing.baseline = timing.baseline;
edit_timing.baseline(edit_timing.baseline == 0) = [];
edit_timing.baseline = edit_timing.baseline-timing.exp_start;

edit_timing.cue = timing.cue;
edit_timing.cue(edit_timing.cue == 0) = [];
edit_timing.cue = edit_timing.cue-timing.exp_start;

edit_timing.reference = timing.reference-timing.exp_start;
edit_timing.ISI = timing.ISI-timing.exp_start;
edit_timing.target = timing.target-timing.exp_start;
edit_timing.response = timing.response-timing.exp_start;

edit_timing.reference_duration = edit_timing.ISI-edit_timing.reference;
edit_timing.ISI_duration = edit_timing.target-edit_timing.ISI;
edit_timing.target_duration = edit_timing.response-edit_timing.target;

edit_timing.reference_temp = edit_timing.reference;
edit_timing.reference_temp(1) = []; edit_timing.reference_temp(end+1) = (timing.exp_end-timing.exp_start);
edit_timing.response_duration = edit_timing.reference_temp-edit_timing.response;

%% Output data
timing_output.headers = {'Trial', 'Condition', 'Orientation', 'Contrast', 'Shape',...
    'Number_Of_Changes', 'RT', 'Keypress', 'Response', 'Reference_Onset', 'Reference_Duration', 'ISI_Onset', 'ISI_Duration',...
    'Target_Onset', 'Target_Duration', 'Response_Onset', 'Response_Duration'};

timing_output.data = table((1:total_trials-1)', info.block', info.orientation', info.contrast', info.shape',...
    info.num_changes', timing.RT', info.keypress', info.response',...
    edit_timing.reference', edit_timing.reference_duration',...
    edit_timing.ISI', edit_timing.ISI_duration', edit_timing.target',....
    edit_timing.target_duration', edit_timing.response', edit_timing.response_duration',...
    'VariableNames', timing_output.headers);

writetable(timing_output.data, strcat(setup.output_dir, sprintf('pp_%s_scan_%s_%s.csv', response{1}, response{2}, datestr(now, 'mm-dd-yyyy_HH-MM'))));

%% 
