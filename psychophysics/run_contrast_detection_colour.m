% function [t, sd] = run_contrast_detection(pp, run, colour)

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

exp.pp = 1; exp.run = 3; exp.colour = 'Luminance';

exp.practice_trials = 10;
exp.experiment_trials = 50;
exp.total_trials = exp.practice_trials + exp.experiment_trials;

% PsychDebugWindowConfiguration; % for debugging (transparent screen).

% ----------------------------------------------------------------------- %
%                            GENERAL SETUP                                %
% ----------------------------------------------------------------------- %

KbName('UnifyKeyNames'); % recommended, uses standard keyboard naming.
space = KbName('space'); escape = KbName('escape'); % set keynames for potential keypresses.

start_directory = pwd; % store origin directory location.

% specify base directory for use with strcat (makes switching PCs easier).
base_directory = 'C:\Users\PS-RMLab\Documents\Kirstie\colour';

% add function directory (and subfolders) path.
addpath(genpath(strcat(base_directory, '\code\psychophysics_code\functions')));

% specify data output directory.
data_storage = strcat(base_directory,'\psychophysics_data\contrast_detection\');

if exp.pp < 10
    exp.file_pp = sprintf('0%d', exp.pp);
end

% check that within the participant, run and feature-specific directory,
% that a file of the specified colour condition does not already exist.
cd(data_storage);
pp_dir = sprintf('pp_%s', exp.file_pp); [~,~] = mkdir(pp_dir); cd(pp_dir);
run_dir = sprintf('run_0%d', exp.run); [~,~] = mkdir(run_dir); cd(run_dir);

if ~isempty(dir(sprintf('%s*', lower(exp.colour))))
    fprintf('already completed this colour run, please increment run number! \n');
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
    dpy.stimLMS.scale = (.05*0.3);
elseif isequal(lower(exp.colour), 'rg')
    % load in pp-specific RG theta values from minimum motion task.
    dpy.theta = mean_rg_theta;
    dpy.stimLMS.dir = [sin(dpy.theta) cos(dpy.theta) 0]';
    dpy.stimLMS.scale = (.04*0.6);
elseif isequal(lower(exp.colour), 'by')
    % load in pp-specific BY theta values from minimum motion task.
    dpy.theta = mean_by_theta;
    dpy.stimLMS.dir = [cos(dpy.theta) cos(dpy.theta) sin(dpy.theta)]';
    dpy.stimLMS.scale = (.15*0.6);
end

% get the LMS-RGB conversions and compute inverse gamma table.
[dpy] = getLMS2RGB_gamma(dpy,1, sensors);

mask_image = (rand(300,300,3));
% ----------------------------------------------------------------------- %
%                            INITALISE QUEST                              %
% ----------------------------------------------------------------------- %

% set the feature-specific threshold estimate.
exp.tGuess = dpy.stimLMS.scale;
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
exp.stim_duration = 0.0833;  % seconds
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
[~, instructions, ~] = get_feature_task_info('Contrast',...
    base_directory, '\', 'contrast_detection');

% load the response-period fixation cross.
response_image = imread('C:\Users\PS-RMLab\Documents\Kirstie\colour\fixation_letters\isi_affinity.png');
first_image = imread('C:\Users\PS-RMLab\Documents\Kirstie\colour\fixation_letters\isi_affinity.png');
second_image = imread('C:\Users\PS-RMLab\Documents\Kirstie\colour\fixation_letters\isi_affinity.png');

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
exp.phi = 0; exp.C =0; exp.A = 0.5;

% generate RFP [500 x 500] image, with these specified parameters.
distances_reference = RFP_generator(exp.phi,exp.C,exp.A);

% modulate image matrix with the specified LMS->RGB colour and scale (contrast).
[reference_RFP_to_draw] = RFP_colour(distances_reference, dpy.stimRGB.dir*dpy.stimRGB.scale);

% convert RFP, feature-specific fixation and response fixation to textures.
RFP_reference = Screen('MakeTexture', window, reference_RFP_to_draw);
response_fix = Screen('MakeTexture', window, response_image);
first_fix = Screen('MakeTexture', window, first_image);
second_fix = Screen('MakeTexture', window, second_image);
mask = Screen('MakeTexture', window, mask_image);

orientations = [deg2rad(15), deg2rad(10), deg2rad(5), deg2rad(0), deg2rad(-5), deg2rad(-10), deg2rad(-15)];

% ----------------------------------------------------------------------- %
InitializePsychSound(1);  % Initialize sounddriver
nrchannels = 1;  freq = 48000;   repetitions = 1; % Set key variables. 

% Wait for device to start & start sound immediately when called. 
startCue = 0; waitForDeviceStart = 1; 

pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels); % Initialise PsychPort audio. 

stimBeep = MakeBeep(1000, 0.02);
% ----------------------------------------------------------------------- %
%                   RUN THRESHOLD STAIRCASE EXPERIMENT                    %
% ----------------------------------------------------------------------- %

for i = 1:exp.total_trials % for each trial in turn:
    
    timing.trial_start{i} = GetSecs; % record time at start of trial.
    
%         PsychPortAudio('FillBuffer', pahandle, [stimBeep]);  % Fill sound buffer with relevant beep frequency.
%         PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);  % Play first stimulus beep.
%         PsychPortAudio('Stop', pahandle, 1,1);

    sound(stimBeep); % play first stimulus beeps (this method seems to resolve the 'multiple-beeps' issue.
    
    counter = 0; % intialise counter variable.
    
    timing.first_sound_start{i} = Screen('Flip', window); % record time at start of first beep presentation. 
    for frame = 1:soundDurationTimeFrames -1  
       
        timing.first_sound_end{i} = Screen('Flip', window, timing.first_sound_start{i} + (waitframes-0.5)*exp.ifi);
        
        while counter < 1
           
            rsp.trial_number{i} = i;
          
            % unless this is the first trial, or a practice trial:
            if i > 1 && i <= exp.practice_trials
                tTest = tTest;
            else
                tTest=QuestQuantile(q); % get the recommended threshold value from Quest (log10 units).
            end
            
            change_value=10^tTest;    %convert threshold value to radians.
            
            % if quest's suggested threshold testing value is greater than this, set the threshold to the maxmimum possible contrast.
            change_value=min(change_value,dpy.stimLMS.maxScale);
            
            rsp.threshold{i} = abs(change_value);
            
            orientation_index = randi(7,1);
            
            % generate the target RFP stimulus with a random orientation selected from
            % the pre-defined possible orientations. 
            distances_target = RFP_generator(orientations(orientation_index),1,exp.A);
            
            % set LMS->RGB colour of target RFP and prepare texture.
            [target_RFP_to_draw] = RFP_colour(distances_target, dpy.stimRGB.dir*change_value);
            RFP_target = Screen('MakeTexture', window, target_RFP_to_draw);
            
            ii = rand(); % generate a random number between 0 and 1.
            
            if ii >= 0.5 % if this random number is more than or equal to 0.5:
                stim1 = RFP_reference; stim2 = RFP_target;
            else
                stim1 = RFP_target; stim2 = RFP_reference;
            end
            counter = counter + 1;
        end
    end
    timing.first_sound_duration{i} = timing.first_sound_end{i} - timing.first_sound_start{i};
    
    % draw the reference RFP and feature-specific fixation image to screen.
    Screen('DrawTexture', window, stim1);
    Screen('DrawTexture', window, first_fix);
    
    timing.first_stim_start{i} = Screen('Flip', window);
    for frame = 1:stimDurationTimeFrames 
        Screen('DrawTexture', window, stim1);
        Screen('DrawTexture', window, first_fix);
        timing.first_stim_end{i} = Screen('Flip', window, timing.first_stim_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.first_stim_duration{i} = timing.first_stim_end{i} - timing.first_stim_start{i};
   
    Screen('DrawTexture', window, mask);
    timing.ib_stim_start{i} = Screen('Flip', window);
    for frame = 1:ibDurationTimeFrames
        
        Screen('DrawTexture', window, mask);
        timing.ib_stim_end{i} = Screen('Flip', window, timing.ib_stim_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.ib_stim_duration{i} = timing.ib_stim_end{i} - timing.ib_stim_start{i};
    
%     PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);   % Play beep.
    sound(stimBeep); % play second stimulus beep.

    timing.second_sound_start{i} = Screen('Flip', window);
    for frame = 1:soundDurationTimeFrames
        
        timing.second_sound_end{i} = Screen('Flip', window, timing.second_sound_start{i} + (waitframes-0.5)*exp.ifi);
    end
    timing.second_sound_duration{i} = timing.second_sound_end{i} - timing.second_sound_start{i};
    
    Screen('DrawTexture', window, stim2);
    Screen('DrawTexture', window, second_fix);
    timing.second_stim_start{i} = Screen('Flip',window);
    
    for frame = 1:stimDurationTimeFrames 
        Screen('DrawTexture', window, stim2);
        Screen('DrawTexture', window, second_fix);
        timing.second_stim_end{i} = Screen('Flip', window, timing.second_stim_start{i} + (waitframes-0.5)*exp.ifi);
        
    end
    timing.second_stim_duration{i} = timing.second_stim_end{i} - timing.second_stim_start{i};
    
    Screen('DrawTexture', window, response_fix);
    timing.response_stim_start{i} = Screen('Flip',window);
    
    timedout = false;

    rsp.response{i} = 'no_response';
    rsp.answer{i} = 'no_answer';

    test_counter = 0;
    rsp.RT{i} = 'NaN'; rsp.keyCode{i} = []; rsp.keyName{i} = 'N/A';
    while ~timedout
        [keyisDown, keyTime, keyCode] = KbCheck;
        
        if(keyisDown)
            if test_counter == 0
                test_counter = 1;
                rsp.RT{i} = keyTime - timing.response_stim_start{i};
                rsp.keyCode{i} = keyCode;
                rsp.keyName{i} = KbName(rsp.keyCode{i});
                
                if isequal(rsp.keyName{i}, 'l') || isequal(rsp.keyName{i}, 'a')
                    [response, response_value, answer, myBeep] =...
                        response_classify_contrast_detection(ii,rsp.keyName{i});
                    
                    rsp.response{i} = response;
                    rsp.answer{i} = answer;
%                     
%                     for sound_loop = 1:1
%                         PsychPortAudio('FillBuffer', pahandle, [myBeep]);                       % Fill sound buffer with relevant beep frequency.
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
end

t=QuestMean(q); sd=QuestSd(q); % calculate final threshold mean and standard deviation.
t = 10^t;  % convert mean threshold from log10 units.

% print final threshold estimate and standard deviation to command window.
fprintf('Final threshold estimate (mean+-sd) is %.2f +- %.2f\n',t,sd)

sca; ListenChar(); ShowCursor; % close screen, renable keypresses and mouse input.

% ----------------------------------------------------------------------- %
%                              DATA OUTPUT                                %
% ----------------------------------------------------------------------- %

% index into relevant storage directories.
cd(data_storage); cd(pp_dir); cd(run_dir);

% plot participants threshold estimations across the experiment.
thresholds = cell2mat(rsp.threshold);
plot(abs(thresholds(exp.practice_trials+1:end)));
ylim([0 max(thresholds)]);ylabel(sprintf('%s',exp.colour));
xlim([1 (exp.experiment_trials)]); xlabel('Trials');

% save this figure as a .png file.
saveas(gcf,sprintf('%s_contrast_detection_threshold_pp_%s_run_0%d.png', lower(exp.colour), exp.file_pp, exp.run));

final_data = [rsp.trial_number', rsp.threshold', rsp.response', rsp.answer', rsp.keyName', rsp.RT', timing.trial_start',...
    timing.first_sound_start', timing.first_sound_end', timing.first_sound_duration',...
    timing.first_stim_start', timing.first_stim_end', timing.first_stim_duration',...
    timing.ib_stim_start', timing.ib_stim_end', timing.ib_stim_duration',...
    timing.second_sound_start',timing.second_sound_end', timing.second_sound_duration',...
    timing.second_stim_start', timing.second_stim_end', timing.second_stim_duration',...
    timing.response_stim_start', timing.response_stim_end', timing.response_stim_duration',...
    timing.trial_end', timing.total_trial_ends'];


data_table = cell2table(final_data, 'VariableNames', {'TrialNumber', 'Threshold', 'ResponseRequired', 'PPResponse',...
    'KeyPressed', 'RT', 'TrialStart', 'Sound1Start', 'Sound1End', 'Sound1Duration', 'Stim1Start', 'Stim1End',...
    'Stim1Duration', 'IBStart', 'IBEnd' ,'IBDuration','Sound2Start', 'Sound2End', 'Sound2Duration', 'Stim2Start', 'Stim2End', 'Stim2Duration',...
    'ResponseStart', 'ResponseEnd', 'ResponseDuration', 'TrialEnd', 'CumulativeTrialEnd'});

% write this data to a .csv file.
writetable(data_table, sprintf('%s_contrast_detection_data_pp_%s_run_0%d.csv',lower(exp.colour), exp.file_pp, exp.run));

% write participant feature-specific thresholds to .mat files for later combination.
contrastdetectionT = t; contrastdetectionSD = sd;
save(sprintf('%s_contrast_detection_thresholds_pp_%s_run_0%d.mat', lower(exp.colour), exp.file_pp, exp.run),...
    'contrastdetectionT', 'contrastdetectionSD');

% save the workspace to a .mat file.
save(sprintf('%s_threshold_workspace_pp_%s_run_0%d', lower(exp.colour), exp.file_pp, exp.run));

cd(start_directory); % return to starting directory.
% ----------------------------------------------------------------------- %