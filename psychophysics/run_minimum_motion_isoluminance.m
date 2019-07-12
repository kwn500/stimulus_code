function [mean_theta, std_theta, all_theta] = run_minimum_motion_isoluminance(pp, colour_condition)

% OLD CODE FOR REFERENCE, SHOULD BE CLEANED IF RE-RUNNING %%

% ---------------- MINIMUM MOTION ISOLUMINANCE TASK --------------------- %

% Participants must shift the theta of 'RG' or 'BY' LMS-modulated stimuli
% (circles with a sine-wave border) to find their point of isoluminance-
% where the motion in the stimulus slows\ stops.

% Participants increase\decrease theta using (the 'U' and 'N' keys in the
% psychophysics lab, and the '1' and '3' buttons in the fMRI scanner).
% Participant confirm their choice pressing the spacebar in the
% psychophysics lab and with the '2' or '4' buttons in the scanner.

% INPUTS;
% pp: participant number (a number, not string).
% colour condition: 'RG' or 'BY' (not case-sensitive).

% ----------------------------------------------------------------------- %

% Remove 'Screen('Preference', 'SkipSyncTests', 1)' when running as an
% experiment (prevents macbook timing issue), also remove
% 'PsychDebugWindowConfiguration' which toggles window transparency.

% KWN (28\02\2018)

% ----------------------------------------------------------------------- %
%                            GENERAL SETUP                                %
% ----------------------------------------------------------------------- %

% Screen('Preference', 'SkipSyncTests', 1); % for macbook PTB timing issue.
% PsychDebugWindowConfiguration; % screen transparency toggle for debug mode.

start_dir = pwd; % get location of starting directory.

% add function directory and subdirectories to path.
addpath(genpath('C:\Users\PS-RMLab\Documents\Kirstie\colour\code\psychophysics_code\functions\'));

wpid = num2str(pp); % convert the participant number to a string for file output.

if pp < 10 % if the participant number is less than 10, format with a leading 0.
    wpid = sprintf('0%d', pp);
end

num_repetitions = 3; % specify number of repetitions desired.

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
%                  GET LMS TO RGB COLOUR SPECIFICATIONS                   %
% ----------------------------------------------------------------------- %

% specify a starting theta for the stimulus.
dpy.theta = 1.7;

% specify screen-specific calibration file.
dpy.calib_file = strcat('C:\Users\PS-RMLab\Documents\Kirstie\colour\calibration_files\B116_VIEWPixx_14032018');

% specify background in RGB space (usually mid-gray [0.5 0.5 0.5]);
dpy.backRGB.dir = [1 1 1]';
dpy.backRGB.scale = 0.5;

% specify stimulus in LMS space, dependent on colour condition.
if isequal(lower(colour_condition), 'rg')
    dpy.stimLMS.dir = [sin(dpy.theta) cos(dpy.theta) 0]';
    dpy.stimLMS.scale = 0.04*(1.5);
    
elseif isequal(lower(colour_condition), 'by')
    dpy.stimLMS.dir = [cos(dpy.theta) cos(dpy.theta) sin(dpy.theta)]';
    dpy.stimLMS.scale = 0.15*(1.5);
end

% specify possible theta rang (90 degrees). 
dpy.maxtheta = pi; dpy.mintheta = pi/2;

% get the LMS-RGB conversions and compute inverse gamma table.
[dpy] = getLMS2RGB_gamma(dpy, 1,sensors);

% ----------------------------------------------------------------------- %
%                         SCREEN\SOUND SETUP                              %
% ----------------------------------------------------------------------- %

% call default PTB settings & draw to max screen.
PsychDefaultSetup(2); screens = Screen('Screens');  screenNumber = max(screens);

% open gray screen window.
[window] = PsychImaging('OpenWindow', screenNumber, (128/255));  

% specify alpha-blending for anti-aliased lines.
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% load the gamma table of the specified screen (gamma correction).
Screen('LoadNormalizedGammaTable', window, dpy.gamma.inverse);

% specify high and low theta beeps.
highBeep = MakeBeep(2000, 0.05); 
lowBeep = MakeBeep(200, 0.1);

ListenChar(2); % disable keypresses to command window.

SetMouse(5000,0); % move mouse to top right corner
HideCursor(window); % make sure the mouse is hidden

% load in the relevant fixation-cross image.
fix_image = imread('C:\Users\PS-RMLab\Documents\Kirstie\colour\fixation_letters\isi_affinity.png');
fix = Screen('MakeTexture', window, fix_image);

% ----------------------------------------------------------------------- %
%                     PRESENT INTRODUCTION SCREEN                         %
% ----------------------------------------------------------------------- %

Screen('TextSize', window, 25); Screen('TextFont', window, 'Arial'); % Set text parameters.

% Present relevant instructions.
line1 = 'Minimum Motion Task\n\n';
line2 = 'Press ''U'' to increase theta\n';
line3 = 'Press ''N'' to decrease theta\n\n';
line4 = 'Press ENTER to start\n\n';
line5 = 'Press SPACEBAR to confirm your selection and advance trial\n';

DrawFormattedText(window, [line1 line2 line3 line4 line5], 'center', 'center', [1 1 1]);

% draw text to screen and wait for keypress.
Screen('Flip',window);
KbWait;

% ----------------------------------------------------------------------- %
%                  PRESENT MINIMIUM-MOTION STIMULUS                       %
% ----------------------------------------------------------------------- %

in = linspace(0,pi*2,15); % specify the phase values of the stimulus.
in_tracker = 1:15; % create range of numbers 1-> length phase values specified.

for phase = 1:length(in) % for each of the phase values specified:
    % create a sine-wave circle with that phase. 
    circle = sine_circle_generator(in(phase));
    
    % store this in the phase-specific position in the storage cell array. 
    circles{phase} = circle;
end

% flip these phase values so we can modulate the shape to appear to move towards and away.
out_tracker = fliplr(in_tracker);

% repeat these matrices, to allow the circle to move more than one phase
% cycle before switching direction.
in_tracker = repmat(in_tracker,1,2);
out_tracker = repmat(out_tracker,1,2);

% combine the matrices to specify the phase cycle of the stimulus.
phase_tracker = [in_tracker, out_tracker];

% ----------------------------------------------------------------------- %
%                  PRESENT MINIMIUM-MOTION STIMULUS                       %
% ----------------------------------------------------------------------- %

i = 1; % initialise repetition counter variable.

% while the number of repetitions is less than the number of repetitions desired:
while i <= num_repetitions
    
    % for each phase desired, up to 500 cycles (this is arbitrary, just to ensure 
    % we definitely present the stimulus for a long enough!)
    for thisPhase = repmat(phase_tracker, 1,500)
        
        circle = circles{thisPhase}; % get the phase-specific sine-wave circle.
        
        if isequal(lower(colour_condition), 'rg') % if a RG run:
            
            % update the RG-LMS definition with the new theta value.
            dpy.stimLMS.dir = [sin(dpy.theta) cos(dpy.theta) 0]';
            
        elseif isequal(lower(colour_condition), 'by') % else, if BY:
            
            % udpdate the BY-LMS definition with the new theta value.
            dpy.stimLMS.dir = [cos(dpy.theta) cos(dpy.theta) sin(dpy.theta)]';
        end
        
        % get the LMS->RGB transformations for the desired theta value.
        [dpy] = getLMS2RGB_gamma(dpy,0, sensors);
        
        % modulate this image with the relevant LMS->RGB modulation specification.
        [circle_draw] = RFP_colour(circle, dpy.stimRGB.dir*dpy.stimRGB.scale);
        
        % create a texture and draw this stimulus to the screen along with a
        % fixation cross.
        sine_circle = Screen('MakeTexture', window, circle_draw);
        Screen('DrawTexture', window, sine_circle);
        Screen('DrawTexture', window, fix);
        
        Screen('Flip', window); % flip the stimulus to the screen.
        
        key_press_counter = 0; % set the current keypress status to zero.
        
        [keyisDown, ~, KeyCode] = KbCheck; % check for a keypress.
        if key_press_counter == 0
            if keyisDown % if a key have been pressed:
                key_press_counter = 1; % update the key press status.
                keyLabel = KbName(KeyCode); % store the keyboard name corresponding to the keypress.
  
                % if the participant pressed 'u' (increase theta):
                if isequal(keyLabel,'u') || isequal(keyLabel, '1!')
                    
                    % increase the current theta value.
                    dpy.theta = dpy.theta + 0.005;
                    
                    % if this theta value exceeds the maximum specified,
                    % just present the maximum.
                    dpy.theta = min(dpy.maxtheta, dpy.theta);
                    
                    % if the theta value has reached the maximum specified:
                    if dpy.theta == dpy.maxtheta
                        
                        % play a high-pitched beep to let the participant
                        % know they have reached the maximum theta.
                        Snd('Play', highBeep);
                    end
                    
                    % else, if the participant pressed 'n' (decrease theta):
                elseif isequal(keyLabel,'n') || isequal(keyLabel, '3#')
                    
                    % decrease the current theta value.
                    dpy.theta = dpy.theta - 0.005;
                    
                    % if this theta value exceeds the minimum specified,
                    % just present the minimum.
                    dpy.theta = max(dpy.mintheta, dpy.theta);
                    
                    % if the theta value has reached the minimum specified:
                    if dpy.theta == dpy.mintheta
                        
                        % play a low-pitched beep to let the participant
                        % know they have reached the minimum theta.
                        Snd('Play',lowBeep);
                    end
                    
                    % else, if the participant pressed the spacebar (end trial):
                elseif isequal(keyLabel, 'space') || isequal(keyLabel, '2@') || isequal(keyLabel, '4$')
                    
                    % if this repetition is less than the number of repetitions desired:
                    if i < num_repetitions
                        
                        % store the theta value in the colour-specific cell array.
                        if isequal(lower(colour_condition),'rg')
                            rg_theta{i} = dpy.theta;
                        elseif isequal(lower(colour_condition),'by')
                            by_theta{i} = dpy.theta;
                        end
                        
                        i = i + 1; % increment repetition number.
                        
                        % randomly generate a new starting theta value between
                        % the minimum and maximum theta values.
                        dpy.theta = dpy.mintheta+(dpy.maxtheta-dpy.mintheta)*rand(1,1);
                        dpy.theta = round(dpy.theta,2);
                        
                        % present the fixation stimulus for 2 seconds as a
                        % break between the repetitions.
                        Screen('DrawTexture', window, fix);
                        Screen('Flip', window);
                        WaitSecs(2);
                        
                        % continue to the next iteration of the loop (i.e.
                        % start the next repetition).
                        continue
                        
                        % else, if this repetition is equal to the number of repetitions desired:
                    elseif i == num_repetitions
                        
                        % close the screen, show the mouse and re-enable keypresses to command window.
                        sca; ShowCursor; ListenChar(0);
                        
                        % index into the relevant data storage directory.
                        cd('C:\Users\PS-RMLab\Documents\Kirstie\colour\psychophysics_data\minimum_motion\');
                        
                        % make and index into the relevant participant directory.
                        pp_dir = sprintf('pp_%s', wpid); [~,~] = mkdir(pp_dir); cd(pp_dir);
                        
                        % save out the final lms modulation-specific theta values in .mat files
                        % (both the repetition-specific thetas and the average across repetitions).
                        if isequal(lower(colour_condition),'rg')
                            rg_theta{i} = dpy.theta;
                            save(sprintf('rg_isoluminance_theta_runs_pp_%s', wpid), 'rg_theta');
                            
                            mean_rg_theta = mean(cell2mat(rg_theta));
                            save(sprintf('rg_isoluminance_theta_average_pp_%s', wpid), 'mean_rg_theta');
                            
                            mean_theta = mean_rg_theta; % output average theta.
                            std_theta = std(cell2mat(rg_theta)); % output standard deviation.
                            all_theta = rg_theta; % output individual theta guesses. 

                        elseif isequal(lower(colour_condition), 'by')
                            by_theta{i} = dpy.theta;
                            save(sprintf('by_isoluminance_theta_runs_pp_%s',wpid), 'by_theta');
                            
                            mean_by_theta = mean(cell2mat(by_theta));
                            save(sprintf('by_isoluminance_theta_average_pp_%s', wpid), 'mean_by_theta');
                            
                            mean_theta = mean_by_theta; % output average theta.
                            std_theta = std(cell2mat(by_theta)); % output standard deviation.
                            all_theta = by_theta; % output individual theta guesses. 
                        end
                        
                        cd(start_dir); % return to starting directory.
                        return % exit the loop (hence end the script).
                    end
                    
                elseif isequal(keyLabel, 'ESCAPE')
                    % close screen, renable keypresses and make mouse cursor visible.
                    sca; ListenChar(); ShowCursor;
                    [mean_theta, std_theta, all_theta] = deal('N\A');
                    cd(start_dir); % return to starting directory.
                    return
                end
            end
        end
    end
end
cd(start_dir); % return to starting directory.
end

