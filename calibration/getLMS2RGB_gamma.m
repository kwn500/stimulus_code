function [dpy] = getLMS2RGB_gamma(dpy, gamma,sensors)

% Loads in a screen calibration file and extracts spectral and gamma 
% measurement information. We then use the spectral data to transform or
% LMS-specified background and stimulus modulations to RGB space for use 
% in PsychToolBox routines. We also compute the inverse gamma table for 
% later screen gamma correction.

% ----------------------------------------------------------------------- %
% ARGUEMENTS:
% cFile- a calibration file containing mean spectral data (with the Jaz,
% 'meanSpectData' [4×16×2048 double], and the gamma measurement values held 
% in 'sumPower' [4×16 double], along with 'nLevels', specifying the number
% of measurement levels for the gamma. 

% dpy is a struct containing the fields:
    % .backRGB.dir- specifying the desired background in LMS space. 
    % .backRGB.scale- specifying the background contrast. 
    
    % .stimLMS.dir- specifying the stimulus modulation in LMS space. 
    % .stimLMS.scale- specifying the stimulus contrast. 
    
% RETURNS:
% fields of dpy:
    % .spectra- [316 x 3] matrix of spectral data from the Jaz calibration.
    
    % .stimLMS.maxScale- the maximum possible contrast for the specified 
        % LMS modulation given the gamut of the screen. 
        
    % .stimRGB.dir- the LMS stimulus modulation in RGB space (-1 1). 
    % .stimRGB.scale- the specified LMS contrast contrast in RGB space. 
    % .stimRGB.maxScale- the maximum possible LMS contrast in RGB space. 

% AUTHOR: KWN (from LW's experiments). 
% DATE: 21/02/2018

% ----------------------------------------------------------------------- %
% SPECTRA 

if gamma == 1;
    % load the screen-specific calibration data.
    cal=load(dpy.calib_file);
    
    % squeeze the mean spectral data into a 2-dimensional array [n x 3].
    sd=(squeeze(cal.meanSpectData(1:3,8,:)))';
    
    % create some variables to specify sizes for the interpolation:
    
    % create a linearly space vector, the same as the size of our original data
    % (this is from 340nm to 1030nm with 2048 points).
    orig=linspace(340,1030,2048);
    
    % create an array of numbers starting at 380nm, increasing in increments of
    % 1 until we reach 361 steps (740nm)- this is the way LW's function
    % requires our data to be formatted.
    dest=380:1:740;
    
    % interpolate our original data into the destination data.
    for t=1:3 % for L, M & S in turn:
        interp_sd(:,t)=interp1(orig,sd(:,t)',dest);
    end
    
    % replace any value less than 0 with 0 (to help remove noise).
    interp_sd(interp_sd<0)=0;
    
    % add this spectral data to the dpy (display) struct.
    dpy.spectra = interp_sd;
end

% calculate the primary (RGB) values needed to create a stimulus defined by
% the specified LMS background and stimulus cone modulations.
[stimLMS, stimRGB] = pry_sensor2primary_KWN(dpy, dpy.stimLMS,sensors);

% add these output arguements to the display struct. 
dpy.stimLMS.maxScale = stimLMS.maxScale;
dpy.stimRGB.dir = stimRGB.dir;
dpy.stimRGB.scale = stimRGB.scale;
dpy.stimRGB.maxScale = stimRGB.maxScale;

% ----------------------------------------------------------------------- %
% GAMMA

% taken from flyTV experiments. 

if gamma == 1;
    % reads in the gamma information (photometer values at different measurement
    % values), computes an inverse gamma function and returns it ready to 
    % pass into Screen('LoadNormalizedGamma');
    igt=computeInverseGammaFromCalibFile(dpy.calib_file);

    % add this inverse gamma table to the dpy struct. 
    dpy.gamma.inverse=igt;
end

end 

