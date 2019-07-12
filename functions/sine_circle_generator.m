function [cImn] = sine_circle_generator(thisPhase)

display.width = 52.4; % width of the screen (cm). 
display.horRes = 1920; % horizontal resolution of the screen (pixels). 
display.viewDist = 57; % viewing distance from the screen (cm). 

% calculate the pixels per degree of visual angle for the specific screen. 
display.ppd = display.horRes/(2*atand((display.width/2)/display.viewDist));

%%
% make a regular 'straight' sine wave grating, then window it, then bend it
% around into a circle. 

% nCyclesPerImage=sqrt(2)/(pi*(7.5/display.ppd)); % spatial frequency- number of cycles. 
nCyclesPerImage = 2.33;
imageSizePix=256; % size of the image (pixels). 
offsetPix=3.8*display.ppd; % offset from the centre (use to set image radius). 
spreadPix=10; % window for the gaussian (how much of the underlying sine wave pattern we 'let through'. 

% make a circle that has a sine-wave border & replicate the matrix so that
% you have a rectangle full of sine wave. 
supportn=linspace(0,2*pi*nCyclesPerImage,imageSizePix)+thisPhase;
sineSpacen=repmat(supportn(:),1,256); 

% now window it with a gaussian: we make a gaussian curve that sits at some
% offset from the center. when we 'curve' the space, the Y axis turns into
% the radius (with zero in the middle) and the X axis turns into polar
% angle. we offset the gaussian from the middle by an amount that
% corresponds to the radius of the thing we want to create. 
expSpacen=linspace(0,256,256);
% eLinn=normpdf(expSpacen,offsetPix,spreadPix);

eLinn=exp((-(expSpacen-offsetPix).^2)/100);

eLinImn=repmat(eLinn(:),1,256);
finalImn=sin(sineSpacen).*eLinImn;

% convert from polar to cartesian co-ordinates. 
cImn=Polar2Im(finalImn,250,'cubic');

end
%% --------------------------------------------------------------------- %%