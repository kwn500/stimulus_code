
function WaitForStartTrigger(site)


% wait for pulse to start scanner
if strcmp(site, 'YNiC') == 1
% % %  change pulse to SIEMENS set up / BM 2017
% % % % % % %     % at YNiC we are using Windows 64 bit and a parallel port for triggers
% % % % % % %     
% % % % % % %     % instantiate the io64 library
% % % % % % %     ioObj = io64;
% % % % % % %     % call the object to initiate the driver
% % % % % % %     status = io64(ioObj);
% % % % % % %     %seleft the base parallel port address; usually LPT1 =  hex 379
% % % % % % %     address = hex2dec('379'); % NB status register not data pins (378)
% % % % % % %     
% % % % % % %     currState = io64(ioObj,address);
% % % % % % %     stateChanged = currState;
% % % % % % %     
% % % % % % %     % wait for a state change to signal that we need to move on
% % % % % % %     while currState == stateChanged
% % % % % % %         stateChanged = io64(ioObj,address);
% % % % % % %     end
% % % % % % %     %go!
% % % % % % %     
      %also, unify the keyboard names for different OSs
    KbName('UnifyKeyNames');
    
    %wait for state change
    while(1)
        [exKeyIsDown, exSecs, exKeyCode] = KbCheck;
        % in Jerusalem, the Siemens scanner gives a 5 at the start of each
        % TR
        if (exKeyIsDown&&exKeyCode(KbName('5%')))
            %only go if we get a 5
            disp('got a 5 -- starting ...')
            break;
        else
            % give some time back to OS
            WaitSecs(0.001);
        end
    end
    %go!
    
    
elseif  strcmp(site, 'Magdeburg') == 1
    %%% initialise your devide
    
    %%% we dont need this - we probably can copy the same thing as below as
    %%% -1 is a fixed setting and will never change so the only thing that
    %%% gets executed is nevertheless the part down there 
    
     %also, unify the keyboard names for different OSs
    KbName('UnifyKeyNames');
    
    %wait for state change
    while(1)
        [exKeyIsDown, exSecs, exKeyCode] = KbCheck;
        % in Magdeburg, the scanner is triggered by keyboard press (T)
        if (exKeyIsDown&&exKeyCode(KbName('T')))
            %only go if we get a T
            disp('got a T -- starting ...')
            break;
        else
            % give some time back to OS
            WaitSecs(0.001);
        end
    end
    %go!
    
    
elseif  strcmp(site, 'Jerusalem') == 1
    %%% Don't need part below any more
   %%%% [junk time0]=deviceUMC('wait for trigger',params.devices.UMCport); 
    
    %also, unify the keyboard names for different OSs
    KbName('UnifyKeyNames');
    
    %wait for state change
    while(1)
        [exKeyIsDown, exSecs, exKeyCode] = KbCheck;
        % in Jerusalem, the Siemens scanner gives a 5 at the start of each
        % TR
        if (exKeyIsDown&&exKeyCode(KbName('5%')))
            %only go if we get a 5
            disp('got a 5 -- starting ...')
            break;
        else
            % give some time back to OS
            WaitSecs(0.001);
        end
    end
    %go!
    
    
    
else
    disp('WARNING!! Site specific triggering code not available! - strating straight away ...')
end
