function [t, sd] = run_feature_thresholds(pp, run, feature, colour)

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

exp.pp = pp; exp.run = run; exp.feature = feature; exp.colour = colour;

% --- RFP ORIENTATION, CONTRAST AND SHAPE QUEST THRESHOLD STAIRCASES ---- %
        
% Participants must judge whether a target RFP is different in a cued
% feature with respect to a reference RFP. Participants respond using with
% the U\N keys to indicate direction of change.

% ----------------------------------------------------------------------- %

% Ensure display.ppd is correct for your screen in 'RFP_generator' before
% running.

% Remove 'Screen('Preference', 'SkipSyncTests', 1)' when running as an
% experiment (prevents macbook timing issue), also remove
% 'PsychDebugWindowConfiguration' which toggles window transparency.

% KWN (21\02\2018)

% ----------------------------------------------------------------------- %
%                       VITAL EXPERIMENT INFORMATION                      %
% ----------------------------------------------------------------------- %

% number of practice trials is subtracted from trials desired- e.g.
% practice_trials: 10, trials_desired: 50 = 40 genuine trials.
exp.practice_trials = 10;
exp.experiment_trials = 75;
exp.total_trials = exp.practice_trials + exp.experiment_trials;

% Screen('Preference', 'SkipSyncTests', 1); % for macbook PTB timing issue.
% PsychDebugWindowConfiguration; % for debugging (transparent screen).
HideCursor; ListenChar(2);

% ----------------------------------------------------------------------- %
%                            GENERAL SETUP                                %
% ----------------------------------------------------------------------- %

KbName('UnifyKeyNames'); % recommended, uses standard keyboard naming.
space = KbName('space'); escape = KbName('escape'); % set keynames for potential keypresses.

start_directory = pwd; % store origin directory location.

% specify base directory for use with strcat (makes switching PCs easier).
base_directory = 'C:\Users\PS-RMLab\Documents\Kirstie\colour';

% add function directory (and subfolders) path.
addpath(genpath(strcat(base_directory, '\code\psychophysics_code\functions\')));

% specify data output directory.
data_storage = strcat(base_directory,'\psychophysics_data\staircase_thresholds\');

exp.file_pp = sprintf('%d', exp.pp);
if exp.pp < 10
    exp.file_pp = sprintf('0%d', exp.pp);
end

% check that within the participant, run and feature-specific directory,
% that a file of the specified colour condition does not already exist.
cd(data_storage);
pp_dir = sprintf('pp_%s', exp.file_pp); [~,~] = mkdir(pp_dir); cd(pp_dir);
run_dir = sprintf('run_0%d', exp.run); [~,~] = mkdir(run_dir); cd(run_dir);
feature_dir = lower(exp.feature); [~,~] = mkdir(feature_dir); cd(feature_dir);

if ~isempty(dir(sprintf('%s*', lower(exp.colour))))
    fprintf('already completed this colour/feature run, please increment run number! \n');
    return
end
cd(start_directory); % return to starting directory.

% ----------------------------------------------------------------------- %
%                     LOAD IN PARTICIPANT THETAS                          %
% ----------------------------------------------------------------------- %
% load in the participant theta value from the isoluminance task.
cd(strcat(base_directory, '\psychophysics_data\minimum_motion\', sprintf('pp_%s', exp.file_pp)));
load(sprintf('rg_isoluminance_theta_average_pp_%s.mat', exp.file_pp)); % rg_theta
load(sprintf('by_isoluminance_theta_average_pp_%s.mat', exp.file_pp)); % by_theta

% load in the 75% correct minimum contrast detection thresholds for the
% participant, for the luminance, rg and by conditions seperately.
% cd(strcat(base_directory, '\psychophysics_data\contrast_detection\', sprintf('pp_%s', exp.file_pp), '\averaged_thresholds'));
% load(sprintf('lum_averaged_contrast_detection_threshold_pp_%s.mat', exp.file_pp)); % lum_min_contrast_mean;
% load(sprintf('rg_averaged_contrast_detection_threshold_pp_%s.mat', exp.file_pp)); % rg_min_contrast_mean;
% load(sprintf('by_averaged_contrast_detection_threshold_pp_%s.mat', exp.file_pp)); % by_min_contrast_mean;

% ----------------------------------------------------------------------- %
%            LOAD STOCKMAN SHARPE 10-DEGREE CONE FUNDAMENTALS             %
% ----------------------------------------------------------------------- %

% load the stockman sharpe 10-degree cone fundamentals.
sensorInit=load('T_cones_ss10'); 

% resample these cone fundamentals so they fit the dimensions needed by 
% getLMS2RGB_gamma and associated functions. 
sensors=interp1(390:1:830,sensorInit.T_cones_ss10',380:1:740);

% replace the NaNs created when we ask to resample from 390 to 380 with zeros). 
sensors(isnan(sensors))=0; 

% ----------------------------------------------------------------------- %
%                  GET LMS TO RGB COLOR SPECIFICATIONS                    %
% ----------------------------------------------------------------------- %

% specify screen-specific calibration file.
dpy.calib_file = strcat(base_directory, '\calibration_files\','B116_VIEWPixx_14032018');

% specify background in RGB space (usually mid-gray [0.5 0.5 0.5]);
dpy.backRGB.dir = [1 1 1]';
dpy.backRGB.scale = 0.5;

% specify stimulus in LMS space, dependent on colour condition.
% scales will eventually be based on individual participant's minimum
% colour-specific contrast detection thresholds.
if isequal(lower(exp.colour),'luminance')
    dpy.stimLMS.dir = [1 1 1]';
    
    % set the contrast level to double the 75% minimum luminance contrast
    % detection threshold for that participant.
    dpy.baseline_contrast = (0.05)*1.5;
%     lum_min_contrast_mean*4;
    dpy.stimLMS.scale = dpy.baseline_contrast;
    
elseif isequal(lower(exp.colour), 'rg')
    % load in pp-specific RG theta values from minimum motion task.
    dpy.theta = mean_rg_theta;
    
    dpy.stimLMS.dir = [sin(dpy.theta) cos(dpy.theta) 0]'; 
    
    % set the contrast level to double the 75% minimum RG contrast detection
    % threshold for that participant.
    dpy.baseline_contrast = (0.04)*1.5;
%     rg_min_contrast_mean*5;
    dpy.stimLMS.scale = dpy.baseline_contrast;
    
elseif isequal(lower(exp.colour), 'by')
    % load in pp-specific BY theta values from minimum motion task.
    dpy.theta = mean_by_theta;
    
    dpy.stimLMS.dir = [cos(dpy.theta) cos(dpy.theta) sin(dpy.theta)]';
    
    % set the contrast level to double the 75% minimum BY contrast detection
    % threshold for that participant.
    dpy.baseline_contrast = (0.15)*1.5;
%     by_min_contrast_mean*15;
    dpy.stimLMS.scale = dpy.baseline_contrast;
end

% get the LMS-RGB conversions and compute inverse gamma table.
[dpy] = getLMS2RGB_gamma(dpy,1, sensors);

% ----------------------------------------------------------------------- %
%                            INITALISE QUEST                              %
% ----------------------------------------------------------------------- %

% set the feature-specific threshold estimate.
if isequal(lower(exp.feature), 'orientation')
    exp.tGuess = 0.15;
elseif isequal(lower(exp.feature), 'contrast')
    if isequal(lower(exp.colour), 'luminance')
%         exp.tGuess = dpy.baseline_contrast-lum_min_contrast_mean;
        exp.tGuess = dpy.baseline_contrast-(0.05);
    elseif isequal(lower(exp.colour), 'rg')
%         exp.tGuess = dpy.baseline_contrast-rg_min_contrast_mean;
        exp.tGuess = dpy.baseline_contrast-(0.04);
    elseif isequal(lower(exp.colour),'by')
%         exp.tGuess = dpy.baseline_contrast-by_min_contrast_mean;
        exp.tGuess = dpy.baseline_contrast-(0.15);
    end
elseif isequal(lower(exp.feature), 'shape')
    exp.tGuess = 0.10;
end

exp.tGuessSd = 0.5; % set standard deviation estimate.

exp.pThreshold=0.75; % set threshold for % correct responses.
exp.beta=3.5;exp.delta=0.01;exp.gamma=0.5; % range=1;

% create struct for measuring threshold.
q = QuestCreate(log10(exp.tGuess),log10(exp.tGuessSd),exp.pThreshold,exp.beta,exp.delta,exp.gamma);
q.normalizePdf=1;

% ----------------------------------------------------------------------- %
%                              SCREEN SETUP                               %
% ----------------------------------------------------------------------- %

% call default PTB settings & draw to max screen.
PsychDefaultSetup(2); screens = Screen('Screens');  screenNumber = max(screens);

% open gray screen window.
[window] = PsychImaging('OpenWindow', screenNumber, (128/255));

% specify alpha-blending for anti-aliased lines.
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% hide mouse cursor & disable keypress output to the command window.
HideCursor; ListenChar(2);

% load the gamma table of the specified screen (gamma correction).
Screen('LoadNormalizedGammaTable', window, dpy.gamma.inverse);

% set script to maximum processing priority.
topPriorityLevel = MaxPriority(window); Priority(topPriorityLevel);

% ----------------------------------------------------------------------- %

% specify stimulus timing in seconds and frames.
exp.stim_duration = 0.15;  % seconds
exp.sound_duration = 0.6;     % seconds
exp.response_duration = 0.9;  % seconds
exp.ib_duration = 0.5; % seconds

exp.ifi = Screen('GetFlipInterval', window);

stimDurationTimeFrames = round(exp.stim_duration / exp.ifi);
soundDurationTimeFrames = round(exp.sound_duration / exp.ifi);
ibDurationTimeFrames = round(exp.ib_duration / exp.ifi);

waitframes = 1; % Numer of frames to wait before re-drawing

% ----------------------------------------------------------------------- %
%                          INSTRUCTION SCREEN                             %
% ----------------------------------------------------------------------- %

% load the relevant feature-specific fixation letter and instructions.
[fixation_image, instructions, reference] = get_feature_task_info(exp.feature,...
    base_directory, '\', 'psychophysics');

% load the response-period fixation cross.
response_image = imread('C:\Users\PS-RMLab\Documents\Kirstie\colour\fixation_letters\interblock_affinity_white.png');

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial'); % set text parameters.

% draw text to screen and wait for spacebar keypress.
DrawFormattedText(window, instructions, 'center', 'center', [1 1 1]);
Screen('Flip',window); [~,KeyCode] = KbWait;

while KeyCode(space) == 0  % while spacebar is not pressed.
    [~,KeyCode] = KbWait; % wait for spacebar response.
    if KeyCode(escape)  % if escape key pressed:
        
        % close screen, renable keypresses and make mouse cursor visible.
        sca; ListenChar(); ShowCursor;
    end
end

% ----------------------------------------------------------------------- %
%                        GENERATE REFERENCE RFP                           %
% ----------------------------------------------------------------------- %

% set desired parameters for RFP.
% orientation (angular rotation), contrast & radial amplitude modulation (shape).
exp.phi = 0; exp.C = 1; exp.A = 0.5;

% generate RFP [500 x 500] image, with these specified parameters.
distances_reference = RFP_generator(exp.phi,exp.C,exp.A);

% modulate image matrix with the specified LMS->RGB colour and scale (contrast).
[reference_RFP_to_draw] = RFP_colour(distances_reference, dpy.stimRGB.dir*dpy.stimRGB.scale);

% convert RFP, feature-specific fixation and response fixation to textures.
RFP_reference = Screen('MakeTexture', window, reference_RFP_to_draw);
fix = Screen('MakeTexture', window, fixation_image);
response_fix = Screen('MakeTexture', window, response_image);

% ----------------------------------------------------------------------- %
InitializePsychSound(1);  % Initialize sounddriver
nrchannels = 2;  freq = 48000;   repetitions = 1; % Set key variables. 

% Wait for device to start & start sound immediately when called. 
startCue = 0; waitForDeviceStart = 1; 

pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels); % Initialise PsychPort audio. 

stimBeep = MakeBeep(1000, 0.02);

% ----------------------------------------------------------------------- %
%                   RUN THRESHOLD STAIRCASE EXPERIMENT                    %
% ----------------------------------------------------------------------- %

for i = 1:exp.total_trials % for each trial in turn:
    
    rsp.trial_number{i} = i;
    
    timing.trial_start{i} = GetSecs; % record time at start of trial.
    
%     PsychPortAudio('FillBuffer', pahandle, [stimBeep; stimBeep]); % Fill sound buffer with relevant beep frequency.
%     PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);   % Play beep.
    
     sound(stimBeep); % play first stimulus beeps (this method seems to resolve the 'multiple-beeps' issue.
     
    timing.first_sound_start{i} = Screen('Flip', window);

    for frame = 1:soundDurationTimeFrames 
        timing.first_sound_end{i} = Screen('Flip', window, timing.first_sound_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.first_sound_duration{i} = timing.first_sound_end{i} - timing.first_sound_start{i};
    
    Screen('DrawTexture', window, RFP_reference);
    Screen('DrawTexture', window, fix);
    
    timing.first_stim_start{i} = Screen('Flip', window);
    for frame = 1:stimDurationTimeFrames 
        Screen('DrawTexture', window, RFP_reference);
        Screen('DrawTexture', window, fix);
        timing.first_stim_end{i} = Screen('Flip', window, timing.first_stim_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.first_stim_duration{i} = timing.first_stim_end{i} - timing.first_stim_start{i};
   
    timing.ib_stim_start{i} = Screen('Flip', window);
    for frame = 1:ibDurationTimeFrames
        timing.ib_stim_end{i} = Screen('Flip', window, timing.ib_stim_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.ib_stim_duration{i} = timing.ib_stim_end{i} - timing.ib_stim_start{i};
    
    %     PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);   % Play beep.
    sound(stimBeep); % play second stimulus beep.
    
    counter = 0;
 
    timing.second_sound_start{i} = Screen('Flip', window);
    for frame = 1:soundDurationTimeFrames - 1
        timing.second_sound_end{i} = Screen('Flip', window, timing.second_sound_start{i} + (waitframes-0.5)*exp.ifi);
        
        while counter < 1
            
            % unless this is the first trial, or a practice trial:
            if i > 1 && i <= exp.practice_trials
                tTest = tTest;
            else
                tTest=QuestQuantile(q); % get the recommended threshold value from Quest (log10 units).
            end
            
            % calcnulate feature-value for target RFP, and categorise this feature
            % change (e.g. clockwise, lower contrast).
            [change_value, test_value, ii, ~] = calculate_value(dpy, dpy.baseline_contrast,...
                tTest, reference, exp.feature, 'psychophysics');
            
            rsp.test_values{i} = test_value;
            rsp.thresholds{i} = abs(change_value);
            
            % create a target RFP, with a change in the relevant feature
            % dimension.
            if isequal(lower(exp.feature), 'orientation')
                distances_target = RFP_generator(test_value,exp.C,exp.A);
            elseif isequal(lower(exp.feature), 'contrast')
                
                rsp.test_values{i} = dpy.baseline_contrast-change_value;
                
                distances_target = RFP_generator(exp.phi,1,exp.A);
                
                dpy.stimLMS.scale = dpy.baseline_contrast - change_value;
                
                % get the LMS-RGB conversions and compute inverse gamma table.
                [dpy] = getLMS2RGB_gamma(dpy,0,sensors);
            elseif isequal(lower(exp.feature), 'shape')
                distances_target = RFP_generator(exp.phi,exp.C,test_value);
            end
            
            % set LMS->RGB colour of target RFP and prepare texture.
            [target_RFP_to_draw] = RFP_colour(distances_target, dpy.stimRGB.dir*dpy.stimRGB.scale);
            RFP_target = Screen('MakeTexture', window, target_RFP_to_draw);
            
            counter = counter + 1;
        end
    end
    timing.second_sound_duration{i} = timing.second_sound_end{i} - timing.second_sound_start{i};
    
    % draw target RFP and feature-specific fixation letter to screen.
    Screen('DrawTexture', window, RFP_target);
    Screen('DrawTexture', window, fix);
    
    timing.second_stim_start{i} = Screen('Flip',window);
    
    for frame = 1:stimDurationTimeFrames 
        Screen('DrawTexture', window, RFP_target);
        Screen('DrawTexture', window, fix);
        timing.second_stim_end{i} = Screen('Flip', window, timing.second_stim_start{i} + (waitframes-0.5)*exp.ifi);
        
    end
    timing.second_stim_duration{i} = timing.second_stim_end{i} - timing.second_stim_start{i};
    
    Screen('DrawTexture', window, response_fix);
    
    timing.response_stim_start{i} = Screen('Flip',window);
    
    timedout = false;

    rsp.response{i} = 'no_response';
    rsp.answer{i} = 'no_answer';
    
    test_counter = 0;
    rsp.RT{i} = 'NaN'; rsp.keyCode{i} = []; rsp.keyName{i} = 'N\A';
    rsp.reference{i} = reference;
    
    reference_value{i} = reference;
    if isequal(lower(exp.feature), 'contrast')
        reference_value{i} = 1\dpy.baseline_contrast;
    end
    
    while ~timedout
        [keyisDown, keyTime, keyCode] = KbCheck;
        
        if(keyisDown)
            if test_counter == 0
                test_counter = 1;
                rsp.RT{i} = keyTime - timing.response_stim_start{i};
                rsp.keyCode{i} = keyCode;
                rsp.keyName{i} = KbName(rsp.keyCode{i});
                
                if isequal(rsp.keyName{i}, 'u') || isequal(rsp.keyName{i}, 'n')
                    [response_value, response, myBeep, answer] =...
                        response_classify(rsp.keyName{i}, ii, exp.feature);
                    
                    rsp.response{i} = response;
                    rsp.answer{i} = answer;
                    
%                     for sound_loop = 1:1
%                         PsychPortAudio('FillBuffer', pahandle, [myBeep; myBeep]);                       % Fill sound buffer with relevant beep frequency.
%                         PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);   % Play beep.
%                         PsychPortAudio('DeleteBuffer');
%                     end
                    sound(myBeep);
                    
                    if i > exp.practice_trials  % if this is not a practice trial:
                        % update quest with presented threshold and participant
                        % response.
                        q=QuestUpdate(q,log10(abs(change_value)),response_value);
                    end
                    
                elseif isequal(rsp.keyName{i}, 'ESCAPE')
                    sca; ListenChar(); ShowCursor; return
                end
            end
            timedout = true; WaitSecs(0.5);
        end
        if((keyTime - timing.response_stim_start{i}) > 3), timedout = true; end
    end
    if(~timedout)
        rsp.RT{i} = keyTime - timing.response_stim_start{i};
        rsp.keyCode{i} = keyCode;
        rsp.keyName{i} = KbName(rsp.keyCode{i});
    end
    
    timing.response_stim_end{i} = GetSecs;
    timing.response_stim_duration{i} = timing.response_stim_end{i} - timing.response_stim_start{i};
    timing.trial_end{i} = GetSecs - timing.trial_start{i};
    
    timing.first_sound_start{i} = timing.first_sound_start{i} - timing.trial_start{i};
    timing.first_sound_end{i} = timing.first_sound_end{i} - timing.trial_start{i};
    timing.first_stim_start{i} = timing.first_stim_start{i} - timing.trial_start{i};
    timing.first_stim_end{i} = timing.first_stim_end{i} - timing.trial_start{i};
    timing.ib_stim_start{i} = timing.ib_stim_start{i} - timing.trial_start{i};
    timing.ib_stim_end{i} = timing.ib_stim_end{i} - timing.trial_start{i};
    timing.second_sound_start{i} = timing.second_sound_start{i} - timing.trial_start{i};
    timing.second_sound_end{i} = timing.second_sound_end{i} - timing.trial_start{i};
    timing.second_stim_start{i} = timing.second_stim_start{i} - timing.trial_start{i};
    timing.second_stim_end{i} = timing.second_stim_end{i} - timing.trial_start{i};
    timing.response_stim_start{i} = timing.response_stim_start{i} - timing.trial_start{i};
    timing.response_stim_end{i} = timing.response_stim_end{i} - timing.trial_start{i};
    
    if i == 1
        timing.total_trial_ends{i} = timing.trial_end{i};
    else
        timing.total_trial_ends{i} = timing.trial_end{i} + timing.total_trial_ends{i-1};
    end
    PsychPortAudio('DeleteBuffer');
end

t=QuestMean(q); sd=QuestSd(q); % calculate final threshold mean and standard deviation.
t = 10^t;  % convert mean threshold from log10 units.

% print final threshold estitame and standard deviation to command window.
fprintf('Final threshold estimate (mean+-sd) is %.2f +- %.2f\n',t,sd)

sca; ListenChar(); ShowCursor; % close screen, renable keypresses and mouse input.

% ----------------------------------------------------------------------- %
%                              DATA OUTPUT                                %
% ----------------------------------------------------------------------- %

% index into relevant storage directories.
cd(data_storage); cd(pp_dir); cd(run_dir); cd(feature_dir);

% plot participants threshold estimations across the experiment.
thresholds = cell2mat(rsp.thresholds);
plot(abs(thresholds(exp.practice_trials+1:end)));
ylim([0 max(thresholds)]);ylabel(sprintf('%s_%s',exp.feature, exp.colour));
xlim([1 (exp.experiment_trials)]); xlabel('Trials');

% save this figure as a .png file.
saveas(gcf,sprintf('%s_%s_feature_thresholds_pp_%s_run_0%d.png', lower(exp.colour), lower(exp.feature), exp.file_pp, exp.run));

final_data = [rsp.trial_number', rsp.reference', rsp.test_values', rsp.thresholds', rsp.response', rsp.answer',...
    rsp.keyName', rsp.RT', timing.trial_start', timing.first_sound_start', timing.first_sound_end', timing.first_sound_duration',...
    timing.first_stim_start', timing.first_stim_end', timing.first_stim_duration',...
    timing.second_sound_start', timing.ib_stim_start', timing.ib_stim_end', timing.ib_stim_duration',...
    timing.second_sound_end', timing.second_sound_duration', timing.second_stim_start', timing.second_stim_end',...
    timing.second_stim_duration', timing.response_stim_start', timing.response_stim_end', timing.response_stim_duration',...
    timing.trial_end', timing.total_trial_ends'];

data_table = cell2table(final_data, 'VariableNames', {'TrialNumber', 'Reference', 'TestValue', 'Threshold',...
    'ResponseRequired', 'PPResponse','KeyPressed', 'RT', 'TrialStart', 'Sound1Start', 'Sound1End', 'Sound1Duration',...
    'Stim1Start', 'Stim1End', 'Stim1Duration', 'IBStart', 'IBEnd' ,'IBDuration', 'Sound2Start', 'Sound2End', 'Sound2Duration', 'Stim2Start', 'Stim2End',...
    'Stim2Duration','ResponseStart', 'ResponseEnd', 'ResponseDuration', 'TrialEnd', 'CumulativeTrialEnd'});

% write this data to a .csv file.
writetable(data_table, sprintf('%s_%s_feature_threshold_data_pp_%s_run_0%d.csv',lower(exp.colour), lower(exp.feature), exp.file_pp, exp.run));

% write participant feature-specific thresholds to .mat files for later combination.
if isequal(lower(exp.feature), 'orientation') % if an orientation task:
    
    % rename the mean and standard deviation for later processing.
    orientationT = t; orientationSD = sd;
    
    % save the participant-specific orientation thresholds to .mat file.
    save(sprintf('%s_%s_feature_thresholds_pp_%s_run_0%d.mat', lower(exp.colour), ...
        lower(exp.feature), exp.file_pp, exp.run), 'orientationT', 'orientationSD');
    
    % repeat a similar process if a contrast or shape feature-task.
elseif isequal(lower(exp.feature), 'contrast')
    contrastT = t; contrastSD = sd;
    save(sprintf('%s_%s_feature_thresholds_pp_%s_run_0%d.mat', lower(exp.colour), ...
        lower(exp.feature), exp.file_pp, exp.run), 'contrastT', 'contrastSD');
elseif isequal(lower(exp.feature), 'shape')
    shapeT = t; shapeSD = sd;
    save(sprintf('%s_%s_feature_thresholds_pp_%s_run_0%d.mat', lower(exp.colour), ...
        lower(exp.feature), exp.file_pp, exp.run), 'shapeT', 'shapeSD');
end

% save the workspace to a .mat file. 
save(sprintf('%s_%s_feature_threshold_workspace_pp_%s_run_0%d', lower(exp.colour),...
    lower(exp.feature), exp.file_pp, exp.run));

cd(start_directory); % return to starting directory.
% ----------------------------------------------------------------------- %
