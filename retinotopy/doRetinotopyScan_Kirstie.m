function doRetinotopyScan_Kirstie(params, site, WaitMRITrigger)
% doRetinotopyScan - runs retinotopy scans
%
% doRetinotopyScan(params, site, WaitMRITrigger)
%
% Runs any of several retinotopy scans
%
% 99.08.12 RFD wrote it, consolidating several variants of retinotopy scan code.
% 05.06.09 SOD modified for OSX, lots of changes.
% 11.09.15 JW added a check for modality. If modality is ECoG, then call
%           ShowScanStimulus with the argument timeFromT0 == false. See
%           ShowScanStimulus for details.

%  16.11.17 AG/BM added site label to control site specific trigger/response
%           code
% 16.11.17 AG/BM added option to wait for trigger or not - 0 or 1 (no or
%           yes)

% defaults
if ~exist('params', 'var'), error('No parameters specified!'); end

% make/load stimulus
stimulus = retLoadStimulus(params);

% loading mex functions for the first time can be
% extremely slow (seconds!), so we want to make sure that
% the ones we are using are loaded.
KbCheck;GetSecs;WaitSecs(0.001);


assignin('base','myParams6',params)


% try
    % check for OpenGL
    AssertOpenGL;
    
    % to skip annoying warning message on display (but not terminal)
    Screen('Preference','SkipSyncTests', 1);
    
    % Open the screen
    params.display                = openScreen(params.display);
    params.display.devices        = params.devices;
    
    % to allow blending
    Screen('BlendFunction', params.display.windowPtr, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Store the images in textures
    stimulus = createTextures(params.display,stimulus);
    
    % If necessary, flip the screen LR or UD  to account for mirrors
    % We now do a single screen flip before the experiment starts (instead
    % of flipping each image). This ensures that everything, including
    % fixation, stimulus, countdown text, etc, all get flipped.
    retScreenReverse(params, stimulus);
    
    % If we are doing ECoG, then add photodiode flash to every other frame
    % of stimulus. This can be used later for syncing stimulus to electrode
    % outputs.
    stimulus = retECOGtrigger(params, stimulus);
    
    assignin('base','myStimulus',stimulus)

    
    for n = 1:params.repetitions,
        % set priority
        %%%%Priority(params.runPriority);
        Priority(1);
        % reset colormap?
        retResetColorMap(params);
        
        % wait for go signal
        onlyWaitKb = false;
        %         pressKey2Begin(params.display, onlyWaitKb, [], [], params.triggerKey); % bdk edit
        
        % bdk edit - wait for mouse click then MRI pulse
        % create and specify properties for parallel port (MRI pulse in)
        WaitMRIPulse = 0;
        if WaitMRIPulse == 1;
            ioObj = io64;
            % call the object to initiate the driver
            status = io64(ioObj);
            %seleft the base parallel port address; usually LPT1 =  hex 379
            address = hex2dec('379'); % NB status register not data pins (378)
        end
        
        clicks = 0;
        Screen('FillRect',params.display.windowPtr, 128); % to clear last texture of stimulus
        Screen('DrawText',params.display.windowPtr,'Waiting for mouse click...', 100, 100);
        Screen('Flip',params.display.windowPtr);
        while clicks < 1
            clicks = GetClicks(params.display.windowPtr);
        end
        
        % wait for scanner trigger if doing that
        % if WaitMRIPulse = 1, waits for MRI trigger before starting
          if WaitMRITrigger == 1
            Screen('DrawText',params.display.windowPtr,'Waiting for mri pulse...', 100, 100);
            Screen('Flip',params.display.windowPtr);
            
            WaitForStartTrigger(site); %module containing site specific trigger code
            %once this fucntion returns, the code moves on
        end
        
        % If we are doing eCOG, then signal to photodiode that expt is
        % starting by giving a patterned flash
        retECOGdiode(params);
        
        % countdown + get start time (time0)
        [time0] = countDown(params.display,params.countdown,params.startScan, params.trigger);
        time0   = time0 + params.startScan; % we know we should be behind by that amount
        
        
        % go
        if isfield(params, 'modality') && strcmpi(params.modality, 'ecog')
            timeFromT0 = false;
        else timeFromT0 = true;
        end
%         keyboard
        [response, timing, quitProg] = showScanStimulus(params.display,stimulus,time0, timeFromT0); %#ok<ASGLU>
        
        % reset priority
        Priority(0);
        
        % get performance
        [pc,rc] = getFixationPerformance(params.fix,stimulus,response);
        fprintf('[%s]: percent correct: %.1f %%, reaction time: %.1f secs',mfilename,pc,rc);
        
        % save
        if params.savestimparams,
            % Save in the current working directory
            filename = [datestr(now,30) '.mat'];
            save(filename);                % save parameters
            fprintf('[%s]:Saving in %s.',mfilename,filename);
        end;
        
        % don't keep going if quit signal is given
        if quitProg, break; end;
        
    end;
    
    % Close the one on-screen and many off-screen windows
    closeScreen(params.display);
    
% catch ME
%     % clean up if error occurred
%     Screen('CloseAll'); setGamma(0); Priority(0); ShowCursor;
%     warning(ME.identifier, ME.message);
% end;


return;








