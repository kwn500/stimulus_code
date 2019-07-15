
function Run_MID_scanner_colour_1 (SubjCode, runNumber)



%first, set default input values if not entered:
if nargin < 2
    runNumber = 1;
end

if nargin < 1
    SubjCode = 'test';
end

%%%%-------------------------------%%%%
%           Define parameters:
% v%%%-------------------------------%%%%

TestFRs = true; %to test whether the screen(s) are at 120 Hz

% If you're using the PROpixx or Viewpixx
UsingVP = true;
useHardwareStereo = false;

% Switch for whether you want the annulus superimposed over the dots:
DrawAnnulus = true;

% Colour parameters:
UseS_cone = true; %Flag for S cone stimulus or L-M cone stimulus (=1 for S cone)

% If you want the 2 deg cone fundamentals, set to true:
Load2Deg = false; % otherwise, the 10 deg ones will be loaded.

% work out output filenames:
dateString = datestr(now);
dateString(dateString == ' ') =  '_';
dateString(dateString == '-') =  '_';
dateString(dateString == ':') =  '_';
%Set aside information
subjectData.experimentDescriptor = 'MID_EventRelated_Colour_1';
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

if ~ismac
    jheapcl; % clear the java heap space.
end

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

% Determine the file name that contains the theta values for isoluminance.
% These must have been specified previously.
% NOTE: Here we assume that you are using the first run of the isolum task (_1).
% It would be unusual if you were using another run...
% AND we are assuming it's for S cone. We can build in other functionality there later if needed.
subjectData.IsoFileName = fullfile('data', ...
    [SubjCode, '_isolum_S_cone_1.mat']);

%%%%-------------------------------%%%%
%         Load the spectra:           %
%%%%-------------------------------%%%%

% Obviously, we have different spectra for the Viewpixx and PROpixx. We load the appropriate ones
% depending on what machine we're running. If we're not on a Vpixx device (or not using it), we just load
% the Viewpixx spectra anyway.
% It is easier to load these here so they can be parsed to the LMS2RGB function.
% Here we also nominate the no. of pixels per degree, which also depends on the display.
% We also determine the pixels per degree and load the subject's isoluminance values here too.
if UsingVP
    Datapixx('Open');
    if Datapixx('IsViewpixx3D')  % If it's the Viewpixx3D
        load('C:\Users\Ryan Maloney\Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')
        PPD = 37; %At 57cm viewing distance, there are 37 pixels/deg on the Viewpixx
        % Load the isoluminance values:
        load(subjectData.IsoFileName);
    elseif Datapixx('IsPropixx') % If it's the PROpixx projector
        % *** Here, we need to load the PROpixx spectra (don't need gamma because it's already linear):
        % NOTE: obviously, we would never load these when on the Viewpixx.
        load('C:\Users\Ryan Maloney\Documents\Calibration\PROpixx_Processed_cal_data_21_4_2016.mat')
        PPD = 46.4; % Pixels per degree on the Propixx, also at 57 cm viewing distance
        % Load the isoluminance values:
        load(subjectData.IsoFileName);
    end
else % If we're not using a Vpixx device, just load the Viewpixx spectra anyway (need spectra from somewhere)...
    load('C:\Users\Ryan Maloney\Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')
    PPD = 37;
    % Load the isoluminance values:
    load(subjectData.IsoFileName);
end

% Save the information about the PC hardware:
subjectData.HostPC = CompName;

% Unify the keyboard names in case we run this on a mac:
KbName('UnifyKeyNames')

%%%%----------------------------------------%%%%
%           Define timing parameters:
%%%%----------------------------------------%%%%

% Load the Optseq optimisation file here.
% This contains the optimum event related design given our conditions/events of interest.
% These must be generated previously & be within the folder 'Optseq' inside the search path.
% They must also have a consistent naming convention.
% Put together the file name for this current subject/run:
OptSeqFile = ['Optseq\', SubjCode, '_MID_colour_design-', ...
    sprintf('%03d', runNumber), '.par'];
% Load the file. This bit is borrowed from convertOptseqtoParfile.m from the MrVista vistadisp package
fid = fopen(OptSeqFile);
OptseqParams = textscan(fid,'%f%f%f%f%s');          % onset,condNum,duration,??,label
fclose(fid);
eventOrder = [OptseqParams{2}, OptseqParams{3}];    % first col gives condition/event, 2nd gives its duration

% Set a bunch of variables important in determining the timing of stimuli:
DummyTRs = 4;       % dummy volumes at the beginning
EventLengthTRs = 1; % number of TRs in an event
TR = 3;             % length of volume (TR), in sec
EvtLengthSec = EventLengthTRs * TR;

% Set up the (9) event conditions:
% These provide the given condition for each event of the scan, and are 'read in' when displaying the stimuli
% These must be specified in exactly the same way when defining events in Optseq
NullISI = 0;      % null ISIs, as determined by Optseq
CD = 1;
CD_control = 2;
IOVD = 3;
IOVD_control = 4;
CD_Scone = 5;
CD_control_Scone = 6;
IOVD_Scone = 7;
IOVD_control_Scone = 8;
Fixation = 9;
Dummy = -1;       % dummy scans will be same stimulus as fixation/Nulls

% Insert the dummy scans into the beginning of 'eventOrder':
eventOrder = [Dummy, DummyTRs*TR; eventOrder];
numEvents = length(eventOrder);                 % the total number of events, including dummies at the beginning

% set aside for saving:
parameters.EventOrderDurn = eventOrder;
parameters.OptseqParams = OptseqParams;

%%%%-------------------------------%%%%
%           Define stimuli:
%%%%-------------------------------%%%%

% Dot colour parameters:

% Load the Stockman/Sharpe cone fundamentals:
if Load2Deg % for 2 deg fundamentals
    load('C:\Users\Ryan Maloney\Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
    %load('StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
else        % or for the 10 deg fundamentals
    load('C:\Users\Ryan Maloney\Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
    %load('StockmanSharpe_10deg_cone_fundamentals_1nm.mat')
end

% Spatial and temporal parameters:

% the effective frame rate PER EYE: since each eye is stimulated on successive video frames
% Remember that the per-eye frame rate on the Viewpixx/PROpixx is 60 Hz
PerEyeFR = RefreshRate/2; %Per eye f.r. is always half the absolute f.r.

% *** contrast ***
% Specify peak dot contrast:
pdc = 1; % peak dot contrast (used in the temporal windowing).

% define raised cosine temporal window over first & last 300 ms of each stimulus event
% Note that this is done according to the PER EYE frame rate, so that each point is sampled twice now, once by each eye
cont = pdc.*ones(1,EvtLengthSec * PerEyeFR);                        % initialize to peak dot contrast (defined above)
win_length = EvtLengthSec * PerEyeFR/4;                             % define window length (10% of total event duration)
cont(1:win_length) = pdc.*0.5.*(1-cos(pi*(1:win_length)./win_length));  % ramp up
cont(end:-1:end-win_length+1) = cont(1:win_length);                     % ... and down
EventContrast = [cont 0 0 0];                                           % add a few extra zeros in there, just in case it displays the stimulus past the final designated frame,
% before being re-synched to the next event (will crash otherwise). Hopefully these 0s will never be needed!

% Define the dot texture, a square-shaped sheet of dots.
% Make the texture the same size as the height of the screen
% (well really we are using many little dot textures, but 'imsize' defines the region they inhabit)
% Note, this should really be the same for both CD and IOVD...
imsize = screenRect(4);

% specify dot size:
dot_sigma_in_degrees = 0.1;             % size of SD (sigma) of the dot profile in degs/vis angle
dot_sigma = dot_sigma_in_degrees * PPD; % sigma in pixels
dotsize = round(dot_sigma * 10);        % make the dots some multiple of sigma
% NOTE: dotsize is simply half the length of the sides of a square patch of pixels that the dot profile is placed within.
% It is dot_sigma that really determines the size of the dots. Obviously, dotsize must be larger than dot_sigma.

% RM 12/5/16
% Make the dots a cosine profile, rather than Gaussian.
% This is because the S-cone dots are too hard to see.

circ_dot_size = dot_sigma * 2.355 * 2; % this approximates the full dot width of the gaussian dots, since there are 2.355 sigma under the curve at FWHM
dot_cos_smooth = circ_dot_size / 3; % make the smoothing 1/3 width of the dots 

%Set up the cosine edges for the dots:
dot_rad = circ_dot_size/2;
x0 = (dotsize+1)/2;
y0 = (dotsize+1)/2;
env = ones(dotsize);
for (ii=1:dotsize)
    for (jj=1:dotsize)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > dot_rad^2)
            env(ii,jj) = 0;
        elseif (r2 > (dot_rad - dot_cos_smooth)^2)
            env(ii,jj) = cos(pi.*(sqrt(r2)-dot_rad+dot_cos_smooth)/(2*dot_cos_smooth))^2;
        end
    end
end

% % make the Gaussian dot profile:
% x = (1-dotsize)/2:(dotsize-1)/2;
% %x = x/PPD; %rescale x into vis angle
% [x,y] = meshgrid(x,x);
% y = -y;
% [a,r] = cart2pol(x,y);
% 
% % This gives us a white-peaked dot (+ve contrast polarity)
% env = exp(-r.^2/(2*dot_sigma.^2)); % gaussian window
% env = env./max(max(abs(env)));     % normalize peak to +/- 1
img = repmat(env,1,1,3); %pre-allocate the image matrix.
img2 = img; % for the negative polarity dots.

% set up the raised cosine annular window. We need this to be at least as big as the screen,
% so that the texture doesn't poke out from behind it.
% specify parameters for the annulus:
inrad = PPD * 0.5;     % inner radius of annulus (result in pixels), for fixation spot
outrad = PPD * 5;     % outer radius of annulus (result in pixels)
% define extent of spatial raised cosine at edge of aperture (in pixels)
cos_smooth = dotsize;   % make it one dot size wide
% This should plonk the window in the middle of the matrix, which is what we want
imsize2 = imsize*2;     % double the texture size
x0 = (imsize2+1)/2;
y0 = (imsize2+1)/2;
J = ones(imsize2);
for (ii=1:imsize2)
    for (jj=1:imsize2)
        r2 = (ii-x0)^2 + (jj-y0)^2;
        if (r2 > outrad^2)
            J(ii,jj) = 0;
        elseif (r2 < inrad^2)
            J(ii,jj) = 0;
        elseif (r2 > (outrad - cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-outrad+cos_smooth)/(2*cos_smooth))^2;
        elseif (r2 < (inrad + cos_smooth)^2)
            J(ii,jj) = cos(pi.*(sqrt(r2)-inrad-cos_smooth)/(2*cos_smooth))^2;
        end
    end
end

%Set aside some more of the stimulus parameters:
parameters.PerEyeFR = PerEyeFR;
parameters.PeakDotContrast = pdc;
parameters.DotSigmaInDeg = dot_sigma_in_degrees;

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
% Define size of spot in pixels:
%fix_spot_border = 10; %for a spot (if you use it instead of a cross): it's black border
%fix_spot_centre = 6; %for the inner portion, that changes in lum with the task

% if you're using a cross instead:
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
ringWidthOuter = ringRadiusOuter - PPD/3;   % 1/3 of a degree thick

% Make the rings. Both are in a 2*2 checkerboard pattern:
fixationRing = double(checkerboard(screenRect(4)/2,1) > 0.5);
% Define the ring:
xx = (1-imsize)/2:(imsize-1)/2;
[xx,yy] = meshgrid(xx,xx);
[~,r] = cart2pol(xx,yy);
% make the alpha mask for the rings, inner and outer.
ring_alphaInner = ((r>ringWidthInner+1) & (r<ringRadiusInner-1)); % Make the alpha mask a tiny bit thinner than the ring itself.
ring_alphaOuter = ((r>ringWidthOuter+1) & (r<ringRadiusOuter-1));

%The fixation cross will randomly alternate between two different shades of grey (above 50% grey)
FixnGreyShift = [0 0 0; 0.7 0.7 0.7];
FixnCol = FixnGreyShift(1,:); % Parsed when fixation cross is drawn: just assign this here now

% generate frames when fixation cross dims / brightens
frameCount = 1;
totalFrames = sum(eventOrder(:,2)) * RefreshRate; % all the frames for the whole experiment

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
%           Set up the screen:
%%%%-------------------------------%%%%

try % Start a try/catch statement, in case something goes awry with the PTB functions
    
    % initialization of the display
    AssertOpenGL;
    % Open PTB onscreen window: We request a 32 bit per colour component
    % floating point framebuffer if it supports alpha-blending. Otherwise
    % the system shall fall back to a 16 bit per colour component framebuffer:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');
    % required for gamma correction through the PsychImaging pipeline:
    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
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
    end
    % No clue what the RegWr and RegWrRd commands are all about, but they are in the demos, so I've included them.
    
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
    % *** this should be switched off here for the colour to work *** 
    %Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE);
    
    % We also want alpha-blending for smooth (anti-aliased) dots...
    % not sure how this will conflict with the above command
    % about the linear superposition of dots... but it doesn't seem to cause problems
    Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Make the textures for the dots.
    % Work out the RGB values for each type of dot.
    % We have 6 kinds of dots: black and white,
    % violet/yellow left eye and violet/yellow right eye.
    % To make life easier we will just run through this 6 times
    
    thetas = Col_results.MeanThetas; % Mean thetas for left and right eyes.
    LMS = [1 1 1]; % For the white dot (both eyes)
    
    % Place LMS values in for left and right eyes
    if UseS_cone
        LMS(2,:) = [cos(thetas(1))/sqrt(2), cos(thetas(1))/sqrt(2), sin(thetas(1))]; %LE
        LMS(3,:) = [cos(thetas(2))/sqrt(2), cos(thetas(2))/sqrt(2), sin(thetas(2))]; %RE
    else
        LMS(2,:) = [-sin(thetas(1)), cos(thetas(1)), 0];
        LMS(3,:) = [-sin(thetas(2)), cos(thetas(2)), 0];
    end
    
    %Go through each pixel in the dot and assign appropriate RGB values using LMS2RGB_Vpixx
    %Compute RGB values for every pixel in LMS cone activation space.
    % The first entry below is for the black/white dots' contrast.
    % Apparently this needs to be only about 10% of the cone contrast, so we do that here:
    % Although we request full contrast here for the S-cones, we only end up using the max anyway.
    DesiredScale = [ColParams.MaxPossibleContrast * 0.1, ColParams.cone_contrast, ColParams.cone_contrast]; %the scale factor (akin to desired cone contrast, I think);
    LMSidx = 0;
    for MakeDots = 1:2:6 % loop across B/W dots, Y/V for L, Y/V for R
        
        % Re-assign the image matrices.
        img = repmat(env,1,1,3);
        img2 = img;
        LMSidx = LMSidx+1; % Counter to index the LMS values for each colour
        stimLMS.dir = LMS(LMSidx,:);
        tic
        for ii = 1:dotsize
            for jj = 1:dotsize
                
                %Do positive polarity dot first:
                stimLMS.scale = env(ii,jj) * DesiredScale(LMSidx); %Work out cone activation for this pixel
                [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                img(ii, jj, :) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
                
                %Now negative polarity dot:
                stimLMS.scale = env(ii,jj) * -DesiredScale(LMSidx); %Note, negative scale for this dot
                [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra);
                img2(ii, jj, :) = stimRGB.dir * stimRGB.scale + 0.5; %scale & add the background to the obtained RGB triplet
                
            end
        end
        
        %Also plop an alpha channel in to the 4th value of the 3rd dimension.
        %Note that these are identical for the 2 polarity dots and describe the transparency values of each
        %1=opaque, 0=fully transparent.
        img(:,:,4) = env;
        img2(:,:,4) = env;
        
        Dots(MakeDots) = Screen('MakeTexture',win,img,[],[],2);
        Dots(MakeDots+1) = Screen('MakeTexture',win,img2,[],[],2);
    end
    
    %     Dots(1) = Screen('MakeTexture',win,env,[],[],2);                % white dot
    %     Dots(2) = Screen('MakeTexture',win,env2,[],[],2);               % black dot
    %Dots(3) = Screen('MakeTexture',win,zeros(size(env2)),[],[],2);  % blank (grey) dots: to present nothing but keep all other parameters the same (on fixn/blank trials)
    
    % Generate the annulus texture:
    AnnImages = 0.5.*ones(imsize2,imsize2,2); %specify RGB matrix of annulus, grey
    AnnImages(:,:,2) = 1.*J; %specify Alpha channel of annulus
    annulus = Screen('MakeTexture',win,AnnImages,[],[],2);
    
    % Generate the Inner ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaInner;
    fixationRingTextureInner = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Generate the Outer ring/fixation lock texture:
    ringMat(:,:,1) = fixationRing;
    ringMat(:,:,2) = ring_alphaOuter;
    fixationRingTextureOuter = Screen('MakeTexture',win,ringMat,[],[],2);
    
    % Preallocate array with destination rectangles:
    % This also defines initial dot locations
    % for the very first drawn stimulus frame:
    texrect = Screen('Rect', Dots(1));
    
    % Now, load some of the pre-genereated stimuli into memory to define the stimulus variables.
    % We set the first entry, Event, to zero so it knows to update the stimulus parameters
    Event = 0;
    % Load CD
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
        (Event, CD, PPD, Dots, texrect, parameters);
    % Load IOVD
    [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
        (Event, IOVD, PPD, Dots, texrect, parameters);
    
    %%%%-------------------------------%%%%
    %         Begin the experiment!
    %%%%-------------------------------%%%%
    
    % set some variables that increment over frames:
    missedFrames = 0;
    runRecord = [];
    %tic;                                        % needed for the keyboard checking sub-routine (see below)
    frame = 1;
    overallFrame = 1;                           % keep track of the total number of displayed frames.
    Event = 1;
    taskCount = 1;                              % keep track of change in fixation cross task
    dimmed = false;
    EventTimes = nan(numEvents,3);              % keep a record of when each event begins/ends
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
                error( 'Run_MID_scanner_colour_1: You quit the program!' );
            end
        end
    end
    
    %sync vbl to startTime
    vbl = Screen('Flip',win); %first flip after pulse
    startTime = vbl;
    %runTimer = GetSecs() - startTime; %stimuli are timed off this start point (should be zero, or very close to it!)
    
    %%%%-------------------------------%%%%
    %         Run through the events
    %%%%-------------------------------%%%%
    
    while Event <= numEvents % loop through all the events (including Null/ISIs & dummies)
        
        % determine correct event end time (in sec)
        if frame == 1
            eventStartVBL = vbl; %from the last vbl
            % determine the correct length of event (includes null events):
            eventEndVBL = eventStartVBL + eventOrder(Event,2); %  - ifi); % don't subtract one ifi: don't need to! (runs on time on new machines)
            EventTimes(Event,1) = eventStartVBL;
            EventTimes(Event,2) = eventEndVBL;
        end
        
        %Determine and load the sequence parameters for the current event (incl Null/ISIs):
        tic
        [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions ...
            (Event, eventOrder(Event,1), PPD, Dots, texrect, parameters);
        toc
        
        % Determine event contrast. In order to ramp up/down the stimulus contrast using the raised cosine,
        % we need to make sure the length of the contrast vector is correct. This needs to be done with each
        % event because the blank/ISI conditions may have lengths longer than the stimulus event lengths
        switch eventOrder(Event,1)
            case {Fixation, NullISI, Dummy}
                EventContrastL = zeros(PerEyeFR * eventOrder(Event,2) + 3, 1); % add 3 just so there is a little bit extra on the end,
                %& we won't get any 'index out of bounds' errors if there is a presentation deadline or 2 missed...
                EventContrastR = EventContrastL;
            case {CD, CD_control, IOVD, IOVD_control, ...
                    CD_Scone, CD_control_Scone, IOVD_Scone, IOVD_control_Scone}
                EventContrastL = EventContrast;
                EventContrastR = EventContrast;
        end
        
        %Draw each stimulus event:
        eventEnd = false;
        while ~eventEnd %should keep iterating across stim frames until vbl >= eventEndVBL
            
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
            Screen('DrawTextures',win, ...
                DotsIdxL(:,mod(f,FramesFullCycle)+1), ...    % determines the dot texture (colour) to be drawn: ***NOTE: now indexed every frame***
                [], ...                                      % texture subpart (not using)
                dstRectsL(:,:,mod(f,FramesFullCycle)+1), ... % determines dot locations for each frame
                [], ...                                      % rotate texture (not using)
                [], EventContrastL(f+1));                    % Final argument is for contrast (modulates global alpha)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);
            %Draw the fixation cross:
            Screen('FillRect',win,FixnCol,fixationCross);
            %Or, if you want to draw a fixation spot (with a border) instead:
            %Screen('DrawDots', win, [centreX centreY], fix_spot_border, [0 0 0], [], 2); %flag of 2 makes them rounded: high-quality anti-aliasing.
            %Screen('DrawDots', win, [centreX centreY], fix_spot_centre, FixnCol, [], 2);
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
                frame = frame + 1;
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
                    error( 'Run_MID_scanner_colour_1: You quit the program!' );
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
            Screen('DrawTextures',win, ...
                DotsIdxR(:,mod(f,FramesFullCycle)+1), ...    % determines the dot texture (colour) to be drawn: ***NOTE: now indexed every frame***
                [], ...                                      % texture subpart (not using)
                dstRectsR(:,:,mod(f,FramesFullCycle)+1), ... % determines dot locations for each frame
                [], ...                                      % Rotate texture (not using)
                [], EventContrastR(f+1));                    % Final argument is contrast (modulates global alpha)
            
            %Superimpose the annulus:
            if DrawAnnulus
                Screen('Blendfunction', win, GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
                Screen('DrawTexture',win,annulus);
                Screen('Blendfunction', win, GL_ONE, GL_ZERO);
            end
            
            %Now draw the fixation ring/lock: requires some fancy tweaks of the alpha settings:
            Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); %need to flip the alpha around again (anti-aliasing)
            Screen('DrawTexture',win, fixationRingTextureInner);
            Screen('DrawTexture',win, fixationRingTextureOuter);
            %Draw the fixation cross:
            Screen('FillRect',win,FixnCol,fixationCross);
            %Or, if you want to draw a fixation spot (with a border) instead:
            %Screen('DrawDots', win, [centreX centreY], fix_spot_border, [0 0 0], [], 2); %flag of 2 makes them rounded: high-quality anti-aliasing.
            %Screen('DrawDots', win, [centreX centreY], fix_spot_centre, FixnCol, [], 2);
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
                    error( 'Run_EEO_scanner_3: You quit the program!' );
                else %any other response, should be a button press (or an accident/error)
                    fixationResponses( overallFrame, 2 ) = resp; %store the subject's response
                end
            end
            
            %%%%------------------------------------------------%%%%
            %               Prepare for next frame:
            %%%%------------------------------------------------%%%%
            
            %Increment frame and overallFrame
            f = f+1;                     % increment counter for next PER EYE frame
            frame = frame + 1;
            overallFrame = overallFrame + 1;
            %runTimer = GetSecs() - startTime % print out the timer, synched to clock
            runTimer = vbl - startTime;  % print out the timer, synched to vbl
            
            % if we've reached the end of the EVENT
            % synch to time rather than frames to overcome any slight
            % deviation from 120 Hz refresh accumulating over the run
            if vbl >= eventEndVBL
                frame   % print out the frame number it ends on, should be RefreshRate * event durn in sec +1
                frame = 1;
                f = 0; % Reset f, which increments at the PER EYE frame rate
                EventTimes(Event,3) = (vbl - startTime); % store the time the event ended (could round this value)
                eventEnd = true;                         % this should terminate current event execution
            end
            
        end %end of iterations across frames
        
        Event = Event + 1; % increment the event.
        
    end % end of iterations across events
    
catch MException
    
    % We throw the error again so the user sees the error description.
    rethrow (MException)
    psychrethrow(psychlasterror);
    ExitGracefully (UsingVP)
    error('Error!')
    
end % End of try/catch statement

% Done. Last flip to take end timestamp and for stimulus offset:
% And print diagnostic info to screen:
vbl = Screen('Flip', win); %, [], [], [], 1);
LastFrame = overallFrame
runTimer = vbl - startTime % should be about 1 frame over correct run time

% Store relevant results and parameters, & write to disc:
subjectData.EventTimes = EventTimes;
subjectData.missedFrames = missedFrames;
subjectData.runTime = runTimer;
subjectData.startTime = startTime;
subjectData.LastFrame = LastFrame;
subjectData.pulseTimes = pulseTimes;
subjectData.responses = fixationResponses;
save(subjectData.filename, 'subjectData', 'parameters');

AverageFlipsPerS = overallFrame / (vbl - startTime) %report the average number of flips per sec. Should be close to the refresh rate.
missedFrames
AccuEr = runTimer-sum(eventOrder(:,2)) %Accumulated error: how much total time the run was off by
ifi
%AccuEr/missedFrames %should be close to 1 frame, in sec

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
function [parameters, DotsIdxL, DotsIdxR, dstRectsL, dstRectsR, FramesFullCycle] = SetUpDotPositions (Event, EventCode, PPD, Dots, texrect, parameters)

% This is where the pre-generated dot positions are loaded and the disparity added.
% (see GenerateIOVDdots.m & GenerateCDdots.m)
% These are picked at random (with replacement) from those already in the folder.
% Set Event to 0 for pre-loading stimuli / assigning memory and storing a record of parameters
% (ie during set-up before we actually draw the stimuli).
% If Event == 0, it will assign memory and store the relevant stimulus parameters for future reference.
% This is because CD /IOVD parameters will not change within a single run of the experiment.
% If Event > 0, it will load and set up stimuli for the current event, and keep a record of the
% file name of the randomly-drawn stimulus for each event.
% The number of frames in a full cycle of motion are also defined here, just in case we run the 2 cues at different frequencies.

% Indicate how many existing sets of dot stimuli are in the subfolder 'stimuli' for loading later.
StimInFolder = 5;

% Define the max. change in disparity (horiz shift of dots) for the 2 cues:
% the numerator is in arcmin (but converted to pixels) akin to the amplitude of the sine wave
IOVD_disparity = 64/60 * PPD;
CD_disparity = 16/60 * PPD;

% Define the MID modulation frequency for the 2 cues, in Hz.
% Note that these must be precisely the same as those in the pre-generated stimulus files.
% This is more to define the appropriate file names for loading, since the frequency itself will be stored in the parameters struct.
IOVD_freq = 1;
CD_freq = 1;

% Convert frequency into a string so we can load the appropriate stimulus:
IOVDstrfreq = num2str(IOVD_freq);
IOVDstrfreq(IOVDstrfreq == '.') = '_'; % replace decimal point with underscore
CDstrfreq = num2str(CD_freq);
CDstrfreq(CDstrfreq == '.') = '_';

% Draw a random value to load the stimuli with:
RandDraw = ceil(rand*StimInFolder);

% Define codes for the conditions. The timing is determined by the Optseq file loaded earlier
NullISI = 0;
CD = 1;
CD_control = 2;
IOVD = 3;
IOVD_control = 4;
CD_Scone = 5;
CD_control_Scone = 6;
IOVD_Scone = 7;
IOVD_control_Scone = 8;
Fixation = 9;
Dummy = -1;

switch EventCode
    
    % *** ---------------------------------- ***
    %               Set up Blanks/fixation
    % *** ---------------------------------- ***
    case {Fixation, NullISI, Dummy}
        
        % If it's a fixation, blank/null/ISI or Dummy event, we simply draw the grey (invisible) dot texture
        % in exactly the same way we draw other stimuli.
        % However to do so, we still need the parameters of one of the other stimuli (ie the number of dots)
        % So we will just load a random CD stimulus and use those
        % The destination rects (dstRects) do not need to change: they can be whatever was last loaded, because what we
        % are drawing is invisible anyway
        
        FileNameCD = fullfile('stimuli',['CD_dots_', CDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
        CD_stim = load(FileNameCD); % includes parameters file and dot positions matrices
        % Just cut the parameters back to a smaller struct:
        CDparams = CD_stim.DotParams;
        
        % Define the number of frames in a full cycle:
        FramesFullCycle = CDparams.FramesFullCycle;
        
        % Set up the dot texture indices: they index invisible dots to both eyes so can be all zeros
        DotsIdxL = repmat(Dots(1),1,CDparams.NumDots);
        DotsIdxR = repmat(Dots(1),1,CDparams.NumDots);
        
        % assign dstRects: destination rect matrices.
        % These can just be zeros because we are drawing blanks anyway
        dstRectsL = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
        dstRectsR = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
        
        % Set aside the random stimulus name for this event:
        if Event > 0
            parameters.RandomStimulusNo{Event} = FileNameCD;
        end
        
        % *** ---------------------------------- ***
        %               Set up CD
        % *** ---------------------------------- ***
    case {CD, CD_control, CD_Scone, CD_control_Scone}
        
        % *** Note that file names for CD chromatic & CD achromatic stimuli are the same
        FileNameCD = fullfile('stimuli',['CD_dots_', CDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
        CD_stim = load(FileNameCD); % includes parameters file and dot positions matrices
        % Just cut the parameters back to a smaller struct:
        CDparams = CD_stim.DotParams;
        
        % Define the number of frames in a full cycle:
        FramesFullCycle = CDparams.FramesFullCycle;
        
        % Assign the pre-generated dot position matrices.
        % First, we duplicate the left eye dot positions to make the right eye dot positions.
        % At this stage all dots are identical for the two eyes (disparity has not been added yet).
        dot_posL = CD_stim.dot_posL;
        dot_posR = CD_stim.dot_posL;
        
        % set aside the disparity
        CDparams.disparity = CD_disparity;
        % For reference also store the disparity in Arcmin
        CDparams.DisparityInArcmin = (CD_disparity/PPD)*60;
        
        % set up a few extra parameters of the sine waves for each cue:
        % the period, in sec
        CDparams.period = 1/CDparams.frequency;
        
        % the angular frequency of the sine waves:
        CDparams.angFreq = 2 * pi * CDparams.frequency;
        
        % Length of the sine wave:
        % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
        % end of the cycle and the beginning of the next cycle.
        % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
        CDparams.t = linspace(CDparams.period/CDparams.FramesFullCycle, CDparams.period, CDparams.FramesFullCycle);
        
        % Now make one full cycle of the sine wave, no matter the frequency
        % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
        % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
        CDparams.SineWave =  CDparams.disparity * sin(CDparams.angFreq * CDparams.t); %disparity * sin(angFreq * t)
        
        % assign dstRects: destination rect matrices
        dstRectsL = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
        dstRectsR = zeros(4, CDparams.NumDots,CDparams.FramesFullCycle);
        
        % Shift dot trajectory: add disparity
        for fr = 1:CDparams.FramesFullCycle
            
            % determine dot coordinates: remember, y position does not change: horiz disparity only
            % Update dot position according to sinsoidal trajectory on each (per-eye) frame
            % Left eye: -sin
            dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ... % the size of the dot texture
                dot_posL(:,1,fr) - CDparams.SineWave(fr), ...  % the x positions
                dot_posL(:,2,fr))';                            % the y positions
            % Right eye: +sin
            dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
                dot_posR(:,1,fr) + CDparams.SineWave(fr), ...
                dot_posR(:,2,fr))';
        end
        
        % If it is the CD control condition, we simply suffle the frames temporally.
        % This gives us a stimulus with the same distribution of disparities but no consistent/smooth
        % change in disparity over time
        % This should be all we need to change for the CD control
        if EventCode == CD_control
            ShuffledFrames = randperm(CDparams.FramesFullCycle);
            dstRectsL = dstRectsL(:,:,ShuffledFrames);
            dstRectsR = dstRectsR(:,:,ShuffledFrames);
        end
        
        % Now we have the dot position indices (the dstRects), define the dot texture indices.
        % These are simply half white (Dots(1) and half black (Dots(2)) in each eye.
        
        % NOTE: since fixing the Right eye IOVD flicker bug (March 2016), DotsIdx matrices are now 2D with num_dots * FramesFullCycle dimensions.
        % So we need to set up the matrix accordingly. We won't do it like the below any more:
        % DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
        % DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2), repmat(Dots(2),1,CDparams.NumDots/2)];
        
        % Because in CD the dots only last one frame & then move, it doesn't really matter how their colour is assigned (as long as there is equivalent
        % numbers of black and white AND they are the same across the eyes.
        % So we can set them up in the same way as IOVD, except they are not later sorted by colour (using MinDist), and they are identical in the 2 eyes:
        
        % ALSO, here is the switch where we assign whether the dots are coloured (S-cone) or achromatic.
        % All spatial and temporal parameters for coloured/achromatic dots are the same.
        
        % For the achromatic dots:
        if EventCode == CD || EventCode == CD_control
            
            %Left eye:
            DotsIdxL = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
            DotsIdxL = repmat(reshape(DotsIdxL,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
            %Right eye:
            DotsIdxR = [repmat(Dots(1),1,CDparams.NumDots/2); repmat(Dots(2),1,CDparams.NumDots/2)]; %Dots(1) = white; Dots(2)=black
            DotsIdxR = repmat(reshape(DotsIdxR,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
            
            % Or for the coloured dots:
            % Note that we measure isoluminance separately for the 2 eyes so each eye has its own pair of dot textures
        elseif EventCode == CD_Scone || Ev3ntCode == CD_control_Scone
            
            %Left eye:
            DotsIdxL = [repmat(Dots(3),1,CDparams.NumDots/2); repmat(Dots(4),1,CDparams.NumDots/2)]; %Dots(3) = purple; Dots(4)=yellow
            DotsIdxL = repmat(reshape(DotsIdxL,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
            %Right eye:
            DotsIdxR = [repmat(Dots(5),1,CDparams.NumDots/2); repmat(Dots(6),1,CDparams.NumDots/2)]; %Dots(5) = purple; Dots(6)=yellow
            DotsIdxR = repmat(reshape(DotsIdxR,1,CDparams.NumDots)', 1, CDparams.FramesFullCycle); % Replicate for each frame
            
        end
        
        % Set aside details of the stimulus parameters.
        % Most of these are already stored in the pre-generated stimulus files,
        % but it is probably worthwhile to have them here again in a single output file
        if Event > 0
            parameters.RandomStimulusNo{Event} = FileNameCD;
        else
            parameters.CDparameters = CDparams; % Store the CD parameters once, during set-up
        end
        
        % *** ---------------------------------- ***
        %               Set up IOVD
        % *** ---------------------------------- ***
    case {IOVD, IOVD_control, IOVD_Scone, IOVD_control_Scone}
        
        FileNameIOVD = fullfile('stimuli',['IOVD_dots_colour_', IOVDstrfreq, '_Hz_stim_', num2str(RandDraw), '.mat']);
        IOVD_stim = load(FileNameIOVD); % includes parameters file and dot positions matrices
        % Just cut the parameters back to a smaller struct:
        IOVDparams = IOVD_stim.DotParams;
        
        % Define the number of frames in a full cycle:
        FramesFullCycle = IOVDparams.FramesFullCycle;
        
        % set aside the disparity
        IOVDparams.disparity = IOVD_disparity;
        % For reference also store the disparity in Arcmin
        IOVDparams.DisparityInArcmin = (IOVD_disparity/PPD)*60;
        
        % set up a few extra parameters of the sine waves for each cue:
        % the period, in sec
        IOVDparams.period = 1/IOVDparams.frequency;
        
        % the angular frequency of the sine waves:
        IOVDparams.angFreq = 2 * pi * IOVDparams.frequency;
        
        % Length of the sine wave:
        % The sine wave should not begin at 0, because then we get issues with it wrapping back to the zero point at both the
        % end of the cycle and the beginning of the next cycle.
        % So begin at the time (in sec), that comes just after zero; this will be the period divided by the no. of frames needed to give a full cycle.
        IOVDparams.t = linspace(IOVDparams.period/IOVDparams.FramesFullCycle, IOVDparams.period, IOVDparams.FramesFullCycle);
        
        % Now make one full cycle of the sine wave, no matter the frequency
        % Of course faster frequencies must 'jump' through a full cycle in fewer frames (this is akin to 'undersampling')
        % Because our screen frame rate is always a fixed function of time, we can't change this (ie by increasing the screen frame rate)
        IOVDparams.SineWave =  IOVDparams.disparity * sin(IOVDparams.angFreq * IOVDparams.t); %disparity * sin(angFreq * t)
        
        % assign dstRects: destination rect matrices
        dstRectsL = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
        dstRectsR = zeros(4, IOVDparams.NumDots,IOVDparams.FramesFullCycle);
        
        % Shift dot trajectory: add disparity.
        % Remember that with IOVD the dot positions for each eye are NOT identical (before the horizontal shifts are added)
        % as they are with CD, that is why we have 2 matrices for the 2 eyes
        % Dot lifetime is determined by how many dot positions are re-randomised with each frame (see GenerateIOVDdots.m)
        
        % Determine a value to modulate the sine wave by.
        % In this way the sinusoidal disparity is either subtracted or added to the dot positions with each frame.
        % If both eyes are added (as in the IOVD_control), then the dots move in the same direction in the 2 eyes.
        % The two columns in SineModulate are for all dots in each eye: col 1 for LE, col 2 for RE.
        
        if EventCode == IOVD_control || EventCode == IOVD_control_Scone
            % Make vector of 1s & -1s
            SineModulate = ones(IOVDparams.NumDots,2);
            SineModulate(1:2:end,:) = -1; % this bit will make stimuli the IOVD control: both eyes have dots moving both directions (and the same pattern)
            
            %**** need random 1/2 to be black, 1/2 white; because color is assigned in alternate manner
            % motion direction cannot also.
            % If motion direction is assigned as above, in an alternate manner, the control has dots of opposite colours moving in
            % opposite directions. We can't have that.
            % So, re-assign SineModulate. Now it is a 2 * num_dots vector, upper half are -1s, lower half +1s
            % Remember that there is no relationship between dot position on the screen & its location in the matrix, so making the SineModulate
            % matrix in this way won't matter; the motion directions will still be assigned randomly in the IOVD control.
            % And because color is assigned in an alternate manner, that means equal numbers of black/white dots will
            % be going both monocular directions.
            
            SineModulate = ones(IOVDparams.NumDots,2);
            SineModulate(1:end/2,:) = -1;
            
            % *** Not doing the below any more ...
            % (Also, define indices to select a random sample of NumDots/2 from each eye.
            % All these dots will be presented to BOTH eyes to ensure the same amount of dots per eye.
            % These indices should stay constant across frames to ensure continuous motion of the dots (so they are assigned here).
            % Equal numbers of dots are sampled from each eye.
            % To do this, we simply re-assign the random sample of dots from the L eye to the R eye.
            % The resulting R eye matrix then has 50% of dots from the L & R eyes (because 50% of dots in the R eye matrix have not been re-assigned).
            % Then we simply make LeyeMatrix = ReyeMatrix so that the same pattern is displayed to both eyes.
            %LERandDotSample = Shuffle([true(IOVDparams.NumDots/2,1); false(IOVDparams.NumDots/2,1)]); )... ***
            
        elseif EventCode == IOVD || EventCode == IOVD_Scone
            
            % subtract left eye, add right eye for MID
            SineModulate = ones(IOVDparams.NumDots,2);
            SineModulate(:,1) = -1; %-1 * sin for the left eye only in actual IOVD; +1 * sin for right eye
        end
        
        for fr = 1:IOVDparams.FramesFullCycle
            
            % determine dot coordinates: remember, y position does not change: horiz disparity only
            % Update dot position according to sinsoidal trajectory on each (per-eye) frame
            % Left eye: -sin
            dstRectsL(:,:,fr) = CenterRectOnPoint(texrect, ...                                    % the size of the dot texture
                IOVD_stim.DotPosBinoc{1}(:,1,fr) + SineModulate(:,1) * IOVDparams.SineWave(fr), ... % x positions
                IOVD_stim.DotPosBinoc{1}(:,2,fr))';                                               % y positions
            
            %Right eye: +sin
            dstRectsR(:,:,fr) = CenterRectOnPoint(texrect, ...
                IOVD_stim.DotPosBinoc{2}(:,1,fr) + SineModulate(:,2) * IOVDparams.SineWave(fr), ...
                IOVD_stim.DotPosBinoc{2}(:,2,fr))';
            
            % *** Don't do the below any more:
            % AS usual we present 2 sets of dots to L/R eyes,
            % But the direction of motion is alternately left or right within a single eye, & this is determined by
            % 'SineModulate' parameter....
            % If it's the IOVD control, we take a random sample of NumDots dots from BOTH eyes,
            % these dots are displayed to BOTH eyes simultaneously
            %if EventCode == IOVD_control % if it's the control, L & R are the same patterns
            %
            %    dstRectsR(:,LERandDotSample,fr) = dstRectsL(:,LERandDotSample,fr); % reassign 50% of dots from L eye into R eye
            %    dstRectsL(:,:,fr) = dstRectsR(:,:,fr);                             % now make LE pattern = RE pattern
            %end
        end
        
        % Now we have the dot position indices (the dstRects), define the dot texture indices.
        % These are simply half white (Dots(1) and half black (Dots(2)) in each eye, same as CD.
        
        % *** This is where the dot color must be assigned in an alternate fashion for L & R eye
        % to make (mostly) sure that neighbouring dots across the eyes are of a different color
        % However you do it, colour must be assigned oppositely to the 2 eyes for this to work
        
        % Previously we only assigned dot color in alternating fashion for IOVD, not the control.
        % This was because we previously used the same set of dots in the 2 eyes for the control,
        % so if color was assigned in an alternating fashion then the dots would flicker between the 2 colours as the
        % dots were presented in alternating fashion to the 2 eyes.
        % Now, because we use a different set of dots, it is ok for the dots to be assigned in alternating manner
        % across the eyes. This also maintains the method used to reduce binocular correlations in the new IOVD algorithm, applying it to the control also.
        
        % *** To solve the R eye flicker bug associated with colours changing across dot lifetime, the DotIdx matrices
        % must now be a 2D matrix with FrmsFullCycle * num_dots dimensions
        % Hence, both DotsIdxL & DotsIdxR must now be indexed by the current frame in the cycle as the dots are drawn
        % (as was done with the dot position matrices). ***
        
        % Generate the alternating matrices of dot colours. Note that they need to alternate.
        % If we just assigned top half of the matrix to black, bottom half white (for eg.), then that would mean
        % mostly white dots would go one direction in the IOVD control stimulus, mostly black in the other direction, for example.
        %(this is because dots are now positioned in alternating strips across the 2 eyes)
        % We can't have that; in the control we would want equal numbers of white/black dots going in either direction.
        
        %Not doing this now:
        %if  EventCode == IOVD_control % if it's control, assign colour in the usual way (as in CD).
        %
        %    DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2), repmat(Dots(2),1,IOVDparams.NumDots/2)];
        %    DotsIdxR = [repmat(Dots(1),1,IOVDparams.NumDots/2), repmat(Dots(2),1,IOVDparams.NumDots/2)];
        %elseif EventCode == IOVD % if it's the actual IOVD stimulus, then we assign colour in alternating fashion:
        
        % ****
        % ALSO, here is the switch where we assign whether the dots are coloured (S-cone) or achromatic.
        % All spatial and temporal parameters for coloured/achromatic dots are the same.
        % ****
        
        % For the achromatic dots:
        if EventCode == IOVD || EventCode == IOVD_control
            
            % Left eye: WBWBWB ....
            DotsIdxL = [repmat(Dots(1),1,IOVDparams.NumDots/2); repmat(Dots(2),1,IOVDparams.NumDots/2)]; % Dots(1) = white; Dots(2)=black
            DotsIdxL = repmat(reshape(DotsIdxL,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
            % Right eye: BWBWBW....
            DotsIdxR = [repmat(Dots(2),1,IOVDparams.NumDots/2); repmat(Dots(1),1,IOVDparams.NumDots/2)];
            DotsIdxR = repmat(reshape(DotsIdxR,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
            
            % Or for the coloured dots:
            % Note that we measure isoluminance separately for the 2 eyes so each eye has its own pair of dot textures
        elseif EventCode == IOVD_Scone || EventCode == IOVD_control_Scone
            
            % Left eye: WBWBWB ....
            DotsIdxL = [repmat(Dots(3),1,IOVDparams.NumDots/2); repmat(Dots(4),1,IOVDparams.NumDots/2)]; % Dots(3) = purple; Dots(4)=yellow
            DotsIdxL = repmat(reshape(DotsIdxL,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
            % Right eye: BWBWBW....
            DotsIdxR = [repmat(Dots(5),1,IOVDparams.NumDots/2); repmat(Dots(6),1,IOVDparams.NumDots/2)]; %Dots(5) = purple; Dots(6)=yellow
            DotsIdxR = repmat(reshape(DotsIdxR,1,IOVDparams.NumDots)', 1, IOVDparams.FramesFullCycle); % Replicate for each frame
            
        end
        
        % Now sort these dot colour matrices according to the indices established earlier that account for dot lifetime
        % Because dot colour was assigned in an opposite manner above, we sort the matrices using the same indices.
        % This retains the opposite polarity assignment across eyes
        % Note that MinDist is stored as a separate matrix in the stimulus file (not within the IOVDparams struct),
        % so it will become a field in the 'IOVD_stim' struct when the stimulus parameters file is loaded:
        DotsIdxR = DotsIdxR(IOVD_stim.MinDist);
        DotsIdxL = DotsIdxL(IOVD_stim.MinDist);
        
        % Set aside details of the stimulus parameters.
        % Most of these are already stored in the pre-generated stimulus files,
        % but it is probably worthwhile to have them here again in a single output file
        if Event > 0
            parameters.RandomStimulusNo{Event} = FileNameIOVD;
        else
            parameters.IOVDparameters = IOVDparams; % Store the IOVD parameters once, during set-up
        end
        
end

end



