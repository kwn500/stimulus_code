function MT_MST_Localiser_Vpixx (SubjCode, runNumber)

% The MT MST localiser version is modified to allow a distinction between MT and MST, based on
% the paper by Fischer, Bulthoff, Logothetis & Bartels (2012):
% 'Visual Motion Responses in the Posterior Cingulate Sulcus: A Comparison to V5/MT and MST'
% We have 3 stimulus conditions, static, and 2 motion: left 1/3 and right 1/3. Ipsilateral responses
% can be considered MST, contralateral only responses can be considered MT.
% R Maloney 7 Nov 2011
%
% Version  _BOLDscreen written 12 September 2013 by R Maloney; for use on the BOLDscreen stimulus display.
% Screen resolution and stimulus parameters scaled for use on the CRS BOLDscreen.
% Stimuli are pre-generated using gen_dot_stim_fMRI_BOLDscreen.m, found in the /stimuli subfolder.
%
% Vpixx version:
% The Vpixx version includes an additional, 4th stimulus condition:
% full-field coherent motion.
% Dot images matrices must be pre-genereated and stored in subfolder /stimuli by gen_dot_stim_fMRI_PROpixx.m
% There were initial problems with the background grey levels of the stimulus overlaid
% on the PTB background screen.
% This was because the dots images are stored as unit8 (probably to save memory/space)
% and were then converted into int16.
% Obviously, a grey level background value of 0.5 is not an integer!
% So now the dots are converted into doubles. 
% Even though this demands more memory this should not cause problems on 
% powerful modern machines with fast processors / high RAM.
% Fixed by RTM 4 June 2016

%first, set default input values if not entered:
if nargin < 2
    runNumber = 1;
end

if nargin < 1
    SubjCode = 'test';
end

%%%%-------------------------------%%%%
%           Define parameters:
%%%%-------------------------------%%%%

TestFRs = true; %to test whether the screen(s) are at 120 Hz

PPD = 46.4;     %there are 46.39 pixels per degree on the PROpixx at 57cm viewing distance (on average over H/W)

% If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

% work out output filenames:
dateString = datestr(now);
dateString(dateString == ' ') =  '_';
dateString(dateString == '-') =  '_';
dateString(dateString == ':') =  '_';
%Set aside information
subjectData.experimentDescriptor = 'MT_MST_Localiser';
subjectData.subjectCode = SubjCode;
subjectData.runNumber = runNumber;

% File name for the data to be saved as:
subjectData.filename = [ SubjCode, '_', ...  %maybe put the /data prefix in later.
    subjectData.experimentDescriptor, ...
    ['_' , dateString, '_'] ...
    num2str( runNumber ), ...
    '_data.mat'];

% Just abort if the file already exists:
if exist(subjectData.filename,'file')
    userResponse = input('WARNING: Run file exists. Overwrite? Enter y or n: ','s');
    if ~strcmp( userResponse, 'y' )
        subjectData = [];
        error('Aborting function! File exists!');
    end
end

jheapcl; % clear the java heap space.

%%%%-------------------------------%%%%
%           Test display/s
%%%%-------------------------------%%%%

% First up, test the frame rate of the display/s if it's the first run.
% It is crucial for 3D/dichoptic presentations to be displayed at 120Hz on the Vpixx device.

if TestFRs && runNumber == 1
    NumScreens = Screen('Screens');
    for ii=1:length(NumScreens)
        ScreenFRs(ii) = Screen('NominalFrameRate', NumScreens(ii)) %print result to command window
    end
    % If any of the displays are not 120 Hz, provide a warning and check them again using the more accurate
    % 'FrameRate' utility. Then abort the program.
    if any(ScreenFRs<120)
        beep
        sprintf('Nominal frame rate is not 120 Hz! \n Obtaining accurate frame rates. Please wait...')
        WaitSecs(1)
        for ii=1:length(NumScreens)
            ScreenFRs(ii) = FrameRate (NumScreens(ii))              %print result to command window
        end
        error ('Check screen display settings!')
    end
end

% Choose the screen: it is usually the max screen no. available.
% Frustratingly, the Shuttle XPC (purchased June 2015) always seems to make the Vpixx display == 1. Not sure why, & can't seem to change it.
% So if we're on that machine, need to -1 from the screen number:
[~, CompName] = system('hostname');                                 % find out the computer name
if strncmpi(CompName, 'pspcawshuttle', length('pspcawshuttle')) ... % normal strcmp not working here, can't figure out why...
        && (length(Screen('Screens'))>1)                            % and there is more than 1 display connected...
    WhichScreen = max( Screen( 'Screens' ) )-1;
else
    WhichScreen = max( Screen( 'Screens' ) );                       % should be right for any other machine!
end

screenRect = Screen('Rect',WhichScreen);                            % get the screen resolution.
centreX = screenRect(3)/2;
centreY = screenRect(4)/2;
RefreshRate = Screen('NominalFrameRate', WhichScreen);

% Save the information about the PC hardware:
subjectData.HostPC = CompName;

% Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

%%%%----------------------------------------%%%%
%           Define timing parameters:
%%%%----------------------------------------%%%%

%Some parameters:
numBlocks = 31; %31 blocks all up
blockLength = 12; %12 seconds. 12 * 31  = 372 s = 6.2 mins 6 mins 12 sec
dotContrast = 1; %contrast of the dots

fixationBlock = 0;
motionBlockL = 1;
motionBlockR = 2;
staticBlock = 3;
motionBlockFF = 4; %full-field, coherent motion.

% fix | {motionL motionR static}! | fix  %% the 24 permutations of the 4 conditions are displayed, separated by fixation blocks
blockOrder = perms([motionBlockL motionBlockR staticBlock motionBlockFF])'; %the 24 possible permutations of the 4 conditions
blockOrder = blockOrder(:,1:6); % just take the first 6 of these permutations, whatever they are.
blockOrder = [blockOrder; zeros(1,6)]; %separate each of the perms with a fixation block
blockOrder = reshape(blockOrder,30,1);
blockOrder = [fixationBlock; blockOrder]; % add in a fixation at the start

%%%%-------------------------------%%%%
%       Set up the fixation task:
%%%%-------------------------------%%%%

% Set up a 'shortcut' list of values to be scanned by KbCheck.
% These are the only keyboard/Current Designs fORP responses we are interested in.
% This should save a little bit of time with each call to 'KbCheck', because it only needs to scan these values, and ignore others
scanListVals = [KbName('1!'), KbName('2@'), KbName('3#'), KbName('4$'), KbName('5%'), KbName('q'), KbName('5')]; % 5 is the trigger, q to quit.
scanList = zeros(1,256); scanList(scanListVals) = 1;

% Set up the fixation cross or spot:
% This is drawn directly to the screen using Screen('FillRect')
% if you're using a cross:
crossWidth = 2;
crossHeight = 10;
fixationCross = fixation_cross(crossWidth,crossHeight,centreX,centreY);

% Make the fixation lock rings:
% We have an inner one around fixation and an outer one right on the edge of screen.
% These could probably be defined as a single texture (rather than 2) but I think that will complicate matters with the alpha-blending settings.
% (they are complicated enough already)
ringRadiusInner = PPD*0.4;                  % ring surrounding fixation
ringRadiusOuter = screenRect(4)/2;          % outer edge (radius) of the ring: the edge of the screen
ringWidthInner = ringRadiusInner - PPD/5;   % 1/5 of a degree thick

% Make the rings. Both are in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
% Define the ring:
xx = (1-screenRect(4))/2:(screenRect(4)-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.

%The fixation cross will randomly alternate between two different shades of grey (above 50% grey)
FixnGreyShift = [0 0 0; 0.7 0.7 0.7];
FixnCol = FixnGreyShift(1,:); % Parsed when fixation cross is drawn: just assign this here now

% generate frames when fixation cross dims / brightens
frameCount = 1;
totalFrames = numBlocks * blockLength * RefreshRate; % all the frames for the whole experiment

% in frames: define parameters of the fixation task
earliestTimeBetweenTasks = 1.5 * RefreshRate;   % minimum time of change
jitterRange = 7.5 * RefreshRate;                % jitter the time between min time & jitterRange+min time

taskFrames = [];
while frameCount < totalFrames
    offset = round((rand(1,1)*jitterRange)+earliestTimeBetweenTasks);
    taskFrames(size(taskFrames,1)+1,1) = frameCount + offset;
    frameCount = frameCount + offset;
end

% To plot a distribution of the fixation task change times:
%FixnChangeDist = diff(taskFrames)
%hist(FixnChangeDist*1000/RefreshRate) %in ms

%%%%-------------------------------%%%%
%           Set up stimuli:
%%%%-------------------------------%%%%

% Set the ratio to adjust stimulus size:
% this is the ratio of the BOLDscreen pixels/deg to the projector pixels/deg; because the BOLDscreen has a greater
% pixel density (79 pixels/deg), we need to increase the number of pixels in the stimuli to match them to the size they would have
% appeared on the projector (a density of about 53 pixels/deg) so 79/53 = 1.4906 or about 1.5

%Propixx has 46.4 PPD, so 46.4/53 = 0.87547

SizeRatio = 46.4/53;

GreyCol = 0.5;

bgColour = [GreyCol GreyCol GreyCol];

textureLoc = [(screenRect(3)-screenRect(4))/2, 0, ...
    screenRect(3)-((screenRect(3)-screenRect(4))/2) screenRect(4)];

% load the dots, stored in a mat file from the output of moving_dots
dotsRaw = load('stimuli\dots_PROpixx.mat');

% Here is were previously, the dots were normalised & converted into int16 format
%dotsRaw = int16(dotsRaw.images / 255); %normalise by 256 to set in range 0-1
%dotsRaw = int16(dotsRaw.images); %normalising seemed to screw up alpha blending

% Now, we make them doubles & normalise, so grey-level background = 0.5, and there is no
% difference between the background of the stimulus & the PTB background outside it :-)
% RTM 4 June 2016
dotsRaw = double(dotsRaw.images)/256;

% apply the contrast to the images
for imageIndex = 1:size(dotsRaw,3)
    dotsRaw(:,:,imageIndex) = dotsRaw(:,:,imageIndex) - GreyCol;
    dotsRaw(:,:,imageIndex) = dotsRaw(:,:,imageIndex) .* dotContrast;
    dotsRaw(:,:,imageIndex) = dotsRaw(:,:,imageIndex) + GreyCol;
end

% ...And we no longer convert back into unit8 format...
%dotsRaw = uint8(dotsRaw);

%make the left and right Alpha channel transparency values:
sz = size(dotsRaw,1);
x = (1-sz)/2:(sz-1)/2;
[x,y] = meshgrid(x,x);
[a,r] = cart2pol(x,y);
z1 = (-a<deg2rad(240)) & (a>deg2rad(120)); %top left half
z2 = (a<deg2rad(240)) & (-a>deg2rad(120)); %bottom left half
LeftAlpha = z1+z2;
%RightAlpha = fliplr(LeftAlpha); %simply flip the left hemifield one to get the right hemifield one

dotsRaw(:,:,end+1) = LeftAlpha; %*255; %now put the alpha channel values at the end of the dots matrix
%dotsRaw(:,:,end+1) = RightAlpha*255;

% pregenerate the random static frames to display during the static condition
staticUpdateRate = 3; %comes to once per second
randomStaticFrames = zeros((blockLength/staticUpdateRate)*4,1); %Not sure what is special about the number 20
minDistance = 20; %?

% generate the first frame to be between 1 and 61
randomStaticFrames(1,1) = round(rand(1,1)*60+1);

for n=2:length(randomStaticFrames)
    randomFrame = round(rand(1,1)*60+1);
    %not sure what this is all about
    invalidValues = mod(randomStaticFrames(n-1,1):randomStaticFrames(n-1,1)+ minDistance,61)+1;
    invalidValues = [invalidValues mod(randomStaticFrames(n-1,1):-1:randomStaticFrames(n-1,1)- minDistance,61)+1];
    
    while ~isempty(find(invalidValues==randomFrame,1))
        randomFrame = round(rand(1,1)*60+1);
    end
    randomStaticFrames(n,1) = randomFrame;
end

% order of motion frames
motionFrameOrder = [1:61 60:-1:2]; %expand, then contract


try
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % Set the color range to be normalised between 0 and 1 (rather than 0-255):
    PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);
    
    % Initialise the Vpixx device:
    if UsingVP        % Enable DATAPixx blueline support, and VIEWPixx scanning backlight for optimal 3D
        
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        % The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses
        
        if Datapixx('IsViewpixx3D')     % If it's the Viewpixx3D
            
            Datapixx('EnableVideoLcd3D60Hz');
            parameters.DisplayType = 'Viewpixx3D'; % set aside the device type for reference
            Datapixx('RegWr');
            
        elseif Datapixx('IsPropixx')    % if it's the Propixx DLP projector
            
            parameters.DisplayType = 'PROpixx';         % set aside the device type for reference
            Datapixx('SetPropixxDlpSequenceProgram',0); % set to normal RGB video processing for driving the LEDs & DLP MMDs
            % Datapixx('RegWr');
            
            % Modify the per-eye crosstalk on the PROpixx.
            % Apparently this cross-talk correction only works when using RB3D video mode,
            % where the red/blue channels contain the left/right eye greyscale images (which we are not using).
            % Datapixx('SetPropixx3DCrosstalkLR', 1);
            % Datapixx('SetPropixx3DCrosstalkRL', 1);
            Datapixx('RegWrRd'); % seem to need to do this after setting 'SetPropixxDlpSequenceProgram' to 0
        end
    else
        parameters.DisplayType = 'Not Vpixx';
    end
    
    % Open an on screen (grey) window and configure the imaging pipeline
    % Info about the 'blueline' mechanism for synching to the 3D glasses:
    % There seems to be a blueline generation bug on some OpenGL systems.
    % SetStereoBlueLineSyncParameters(windowPtr, windowRect(4)) corrects the
    % bug on some systems, but breaks on other systems.
    % We'll just disable automatic blueline, and manually draw our own bluelines!
    if useHardwareStereo
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5, [], [], [], 1); %flag of 1 engages stereomode
        SetStereoBlueLineSyncParameters(win, windowRect(4)+10);
    else
        [win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);
    end
    
    %Define the 'blue line' parameters
    blueRectLeftOn   = [0,                 windowRect(4)-1, windowRect(3)/4,   windowRect(4)];
    blueRectLeftOff  = [windowRect(3)/4,   windowRect(4)-1, windowRect(3),     windowRect(4)];
    blueRectRightOn  = [0,                 windowRect(4)-1, windowRect(3)*3/4, windowRect(4)];
    blueRectRightOff = [windowRect(3)*3/4, windowRect(4)-1, windowRect(3),     windowRect(4)];
    
    HideCursor;
    % Normally we might do the gamma correction here using
    % PsychColorCorrection('SetEncodingGamma', win, [1/R_gamma, 1/G_gamma, 1/B_gamma]);
    % But since the PROpixx has built-in linear channels there is probably no need.
    % raise priority level:
    priorityLevel=MaxPriority(win); Priority(priorityLevel);
    % Query the screen refresh rate:
    ifi = Screen('GetFlipInterval',win); %in sec
        
    % Set the alpha-blending:
    % We want a linear superposition of the dots should they overlap:
    % Just like the Gabors in GarboriumDemo.m (see there for further info).
    % Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    %Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % convert the dot images to PTB textures
    for frame = 1:61
        dotsTexturesL(frame,1) = Screen('MakeTexture',win,dotsRaw(:,:,[frame end])); %left hemifield %end (62) gives the alpha channel
        %dotsTexturesR(frame,1) = Screen('MakeTexture',win,dotsRaw(:,:,[frame end]));  %right hemifield %end (63) gives the alpha channel
        dotsTextures(frame,1) = Screen('MakeTexture',win,squeeze(dotsRaw(:,:,frame))); %static
    end
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    
    %%%%-------------------------------%%%%
    %         Begin the experiment!
    %%%%-------------------------------%%%%
    
    % set some variables that increment over frames:
    userStop = false;
    frame = 1;
    overallFrame = 1;
    block = 1;
    taskCount = 1;
    dimmed = false;
    randomStaticFrameIndex = 1;
    missedFrames = 0;
    
    BlocktTimes = nan(numBlocks,3);              % keep a record of when each event begins/ends
    pulseTimes=[];                              % this is where we will store the time of each scanner pulse
    fixationResponses = zeros(totalFrames,2);   % Store the button presses & colour status of the fixation cross
    
    f = 0; % this value increases with each iteration, but on a per eye basis only
    % And it determines the phase of the motion at the given frame (ie place in the cycle)
    
    % wait for the trigger, obtain a starting time when received
    jheapcl;                % clear the java heap space again.
    KbCheck();              % take a quick KbCheck to load it now & flush any stored events
    startTime = GetSecs();  % load GetSecs and assign startTime into memory: doing this takes some time on the first occasion
    
    %%%%-------------------------------%%%%
    %         Wait for the trigger
    %%%%-------------------------------%%%%
    
    % Wait for either a TTL pulse from the scanner,
    % or the user presses the '5' key to begin.
    Screen('TextFont',win, 'Arial');
    Screen('TextSize',win, 20);
    % Display welcome screen to both eyes:
    if useHardwareStereo
        Screen('SelectStereoDrawBuffer', win, 0);   % flag of 0= left eye
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('SelectStereoDrawBuffer', win, 1);   % flag of 1= right eye
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('Flip', win); % , [], [], [], 1);
    else
        DrawFormattedText(win,['Welcome ' SubjCode '. \n \nWaiting for trigger... \n \nPress '' 5 '' on keyboard to override, \n \nor '' q '' to quit at any time.'], ...
            'center', 'center', 1);
        Screen('Flip', win); %, [], [], [], 1);
    end
    
    % if in the scanner, need to wait for the trigger before beginning
    start = false;
    while ~start
        StartResp = CheckForResponses (scanList);
        if StartResp ~= 0
            % If it's a five, it means we've either received the scanner trigger,
            % or someone has pressed '5' on the keyboard, so begin
            if StartResp == 5
                start = true; %exit this while loop & return 'startTime': continue playback to begin the program!
                %startTime = GetSecs();
            elseif StartResp == 81 %they've pressed 'q' to quit
                ExitGracefully (UsingVP);
                error( 'MT_MST_Localiser_Vpixx: You quit the program!' );
            end
        end
    end
    
    %sync vbl to startTime
    vbl = Screen('Flip',win); %first flip after pulse
    startTime = vbl;
    %runTimer = GetSecs() - startTime; %stimuli are timed off this start point (should be zero, or very close to it!)
    
    
    %%%%-------------------------------%%%%
    %         Run through the blocks
    %%%%-------------------------------%%%%
    
    while block <= numBlocks && ~userStop
        
        % if this is a frame where the fixation dimming task should change
        if overallFrame == taskFrames(taskCount,1)
            % swap between true / false ie 0=black, 1=grey
            dimmed = mod(dimmed+1,2);
            taskCount = taskCount + 1;
            FixnCol = FixnGreyShift(dimmed+1, :); %will be either dark or light grey
        end
        
        % Keep record of what that fixation task status was
        % make the dimmed output either 1 or 2 (rather than 0 or 1)
        fixationResponses(overallFrame,1) = dimmed+1; %1=black, 2=grey
        
        % *** Draw the dots! ***
        
        % Select left-eye image buffer for drawing:
        if useHardwareStereo
            Screen('SelectStereoDrawBuffer', win, 0);
        end
        
        %%%%------------------------------------------------%%%%
        %               Draw left eye stimulus:
        %%%%------------------------------------------------%%%%
        
        Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        switch blockOrder(block)
            
            case fixationBlock    
                % do nothing. No stimuli displayed (just fixation)
            
            case motionBlockFF %full-field coherent motion.
                Screen('DrawTexture',win,dotsTextures(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc);
            
            case motionBlockL
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 20, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);
                Screen('DrawTexture',win,dotsTexturesL(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc);
                
            case motionBlockR
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 16, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                %to draw the Right hemi stimulus, we simply rotate the left hemi stimulus by 180 deg; this means there is 1/3 less textures open
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);
                Screen('DrawTexture',win,dotsTexturesL(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc,180); %rotationAngle=180deg
                
            case staticBlock
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 16, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);              
        end
        
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        %Draw the fixation cross:
        Screen('FillRect',win,FixnCol,fixationCross);
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        %Draw blue lines:
        Screen('FillRect', win, [0, 0, 1], blueRectLeftOn);
        Screen('FillRect', win, [0, 0, 0], blueRectLeftOff);
        
        % Select right-eye image buffer for drawing:
        if useHardwareStereo
            Screen('SelectStereoDrawBuffer', win, 1);
        else %Not sure what would happen if this 'else' was actually true. I suspect something would go wrong, as stim would be presented twice.
            %But this is more or less how the Vpixx people do it in their 'DatapixxImagingStereoDemo'
            Screen('DrawingFinished', win);
            [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
            %frame = frame + 1;
            overallFrame = overallFrame + 1;
            %runTimer = GetSecs() - startTime %print out the timer, synched to clock
            runTimer = vbl - startTime; %print out the timer, synched to vbl
            if missed > 0
                missedFrames = missedFrames + 1;
            end
        end
        
        %%%%------------------------------------------------%%%%
        %         Check for button presses/ scanner pulses
        %%%%------------------------------------------------%%%%
        resp = CheckForResponses (scanList);
        
        % if resp > 0, something has happened
        if resp ~= 0
            % if it was the trigger: store the time
            if resp == 5
                pulseTimes(size(pulseTimes,1)+1,1) = GetSecs() - startTime; %could use round(GetSecs())
                % if the 'q' key is pressed, quit the program
            elseif resp == KbName( 'q' )
                ExitGracefully (UsingVP)
                error( 'MT_MST_Localiser_Vpixx: You quit the program!' );
            else %any other response, should be a button press (or an accident/error)
                fixationResponses( overallFrame, 2 ) = resp; %store the subject's response
            end
        end
        
        % if this is a frame where the fixation dimming task should change
        if overallFrame == taskFrames(taskCount,1)
            % swap between true / false ie 0=black, 1=grey
            dimmed = mod(dimmed+1,2);
            taskCount = taskCount + 1;
            FixnCol = FixnGreyShift(dimmed+1, :); %will be either dark or light grey
        end
        
        % Keep record of what that fixation task status was
        % make the dimmed output either 1 or 2 (rather than 0 or 1)
        fixationResponses(overallFrame,1) = dimmed+1; %1=black, 2=grey
        
        %%%%------------------------------------------------%%%%
        %               Draw right eye stimulus:
        %%%%------------------------------------------------%%%%
        
        Screen('Blendfunction', win, GL_SRC_ALPHA, GL_ONE); %turn on alpha blending for lin superposition: see GarboriumDemo.m
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        switch blockOrder(block)
            
            case fixationBlock
                % do nothing. No stimuli displayed (just fixation)
                
            case motionBlockFF %full-field coherent motion.
                Screen('DrawTexture',win,dotsTextures(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc);
            
            case motionBlockL
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 20, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);
                Screen('DrawTexture',win,dotsTexturesL(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc);
                
            case motionBlockR
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 16, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                %to draw the Right hemi stimulus, we simply rotate the left hemi stimulus by 180 deg; this means there is 1/3 less textures open
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);
                Screen('DrawTexture',win,dotsTexturesL(motionFrameOrder(mod(frame-1,120)+1),1),[],textureLoc,180); %rotationAngle=180deg
                
            case staticBlock
                
                % change to a new random frame for the static regions
                if mod(frame,60*staticUpdateRate) == 0
                    %if the index exceeds 16, take it back to 1:
                    randomStaticFrameIndex = mod(randomStaticFrameIndex,16) + 1; %randomStaticFrameIndex + 1;
                end
                Screen('DrawTexture',win,dotsTextures(randomStaticFrames(randomStaticFrameIndex,1),1),[],textureLoc);              
        end
        
        %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
        Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
        Screen('DrawTexture',win, fixationRingTextureInner);
        %Draw the fixation cross:
        Screen('FillRect',win,FixnCol,fixationCross);
        Screen('BlendFunction', win, GL_ONE_MINUS_DST_ALPHA, GL_DST_ALPHA); %Done drawing so flip alpha back again
        Screen('BlendFunction', win, GL_ONE, GL_ZERO);
        %Draw blue lines:
        Screen('FillRect', win, [0, 0, 1], blueRectRightOn);
        Screen('FillRect', win, [0, 0, 0], blueRectRightOff);
        
        Screen('DrawingFinished', win);
        
        [vbl , ~ , ~, missed] = Screen('Flip', win, vbl + (ifi*0.5)); %, [], [], 1); %update display on next refresh (& provide deadline)
        
        %keep record of any missed frames:
        if missed > 0
            missedFrames = missedFrames + 1;
        end
        
        %%%%------------------------------------------------%%%%
        %         Check for button presses/ scanner pulses
        %%%%------------------------------------------------%%%%
        resp = CheckForResponses (scanList);
        
        % if resp > 0, something has happened
        if resp ~= 0
            % if it was the trigger: store the time
            if resp == 5
                pulseTimes(size(pulseTimes,1)+1,1) = GetSecs() - startTime; %could use round(GetSecs())
                % if the 'q' key is pressed, quit the program
            elseif resp == KbName( 'q' )
                ExitGracefully (UsingVP)
                error( 'MT_MST_Localiser_Vpixx: You quit the program!' );
            else %any other response, should be a button press (or an accident/error)
                fixationResponses( overallFrame, 2 ) = resp; %store the subject's response
            end
        end
        
        if frame == 1
            blockStartVBL = vbl;
            blockEndVBL = blockStartVBL + (blockLength - ifi);
        end
        
        %Increment frame and overallFrame
        frame = frame + 1;
        overallFrame = overallFrame + 1
        
        % if we've reached the end of the block
        % synched to time rather than frames to overcome any slight
        % deviation from 120 Hz refresh cumulating over the run
        if vbl >= blockEndVBL
            frame = 1;
            block = block + 1;
            blockTimes(block,1) = round(vbl - startTime);
        end
        
    end % end of iterations across frames
    
catch MException
    
    % We throw the error again so the user sees the error description.
    rethrow (MException)
    psychrethrow(psychlasterror);
    ExitGracefully (UsingVP)
    error('Error!')
    
end% End of try/catch statement

% Done. Last flip to take end timestamp and for stimulus offset:
% And print diagnostic info to screen:
vbl = Screen('Flip', win); %, [], [], [], 1);
LastFrame = overallFrame
runTimer = vbl - startTime % should be about 1 frame over correct run time

% Store relevant results and parameters, & write to disc:
subjectData.BlockTimes = BlocktTimes;
subjectData.missedFrames = missedFrames;
subjectData.runTime = runTimer;
subjectData.startTime = startTime;
subjectData.LastFrame = LastFrame;
subjectData.pulseTimes = pulseTimes;
subjectData.responses = fixationResponses;

AverageFlipsPerS = overallFrame / (vbl - startTime) %report the average number of flips per sec. Should be close to the refresh rate.
missedFrames
ifi
%AccuEr/missedFrames %should be close to 1 frame, in sec

% Save the data:
save(subjectData.filename, 'subjectData', 'parameters');

ExitGracefully (UsingVP)


end %end of main function

%Now for the sub-functions:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resp = CheckForResponses(scanList)
% resp is either 1,2,3, or 4 if coming from the Current Designs fibre optic response pad (fORP-932),
% or 5 from the scanner, or the keycode from the keyboard (all as numbers (doubles), not chars)
% Note that all inputs from the fORP-932 are received in the same way as the '1!', '2@', '3#', '4$' & '5%' American qwerty keys
% (not the numeric keypad).
% Find out the codes of different keys using 'KbName' eg KbName(53) = '5%'
% On the response pad, buttons are assigned in a clockwise direction beginning at the right (3 o'clock)

resp = 0;
[ keyIsDown, ~, keyCode ] = KbCheck([],scanList);
%if toc > 0.2 %wait at least 200 ms to re-check for response
if keyIsDown
    
    %tic;
    
    % return a trigger input (coded as a 5) or a key press of either of the 2 '5' keys to begin the program
    if find(keyCode,1) == 53; %KbName('5%') %=53 that's the trigger, or someone's pressed the qwerty 5/% key
        resp = 5;
    elseif find(keyCode,1) == 101; %KbName('5') %=101 someone's pressed 5 on the numeric keypad, that's ok
        resp = 5;
    elseif find(keyCode,1) == 49; %KbName('1!') %=49 button 1 = blue/right button has been pressed
        resp = 1;
    elseif find(keyCode,1) == 50; %KbName('2@') %=50 button 2 = yellow/bottom button has been pressed
        resp = 2;
    elseif find(keyCode,1) == 51; %KbName('3#') %=51 button 3 = green/left button has been pressed
        resp = 3;
    elseif find(keyCode,1) == 52; %KbName('4$') %=52 button 4 = red/top button has been pressed
        resp = 4;
    elseif find(keyCode,1) == 81; %KbName('q') %=81 someone has pressed 'q' key to quit the program!
        resp = 81; %to abort/exit
    else
        resp = 999; %if it's anything else, there's been some error/key pressed accidentally
    end
end
%end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ExitGracefully (UsingVP)
%...need to shut everything down here...

%turn off the prioritisation:
Priority( 0 ); %restore priority

if UsingVP        % close down the ViewPixx or ProPixx
    Datapixx('DisableVideoScanningBacklight');
    if Datapixx('IsViewpixx3D')
        Datapixx('DisableVideoLcd3D60Hz');
    end
    Datapixx('RegWr');
    %Datapixx('Close'); %closing it here might cause it to crash?
end

%Close down the screen:
Screen('CloseAll')

Datapixx('Close'); %closing the Datapixx here (after closing the screen) might stop it from crashing

%Bring back the mouse cursor:
ShowCursor();

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rect = fixation_cross(width,height,centrex,centrey)
% $Id: fixation_cross.m 4 2007-06-23 11:13:19Z damienm $

rect = zeros(4,2);

width = width/2;
height = height/2;

rect(1,1) = -height;
rect(2,1) = width;
rect(3,1) = height;
rect(4,1) = -width;

rect(1,2) = -width;
rect(2,2) = height;
rect(3,2) = width;
rect(4,2) = -height;


rect(1:2:4,:) = rect(1:2:4,:) + centrex;
rect(2:2:4,:) = rect(2:2:4,:) + centrey;
end

