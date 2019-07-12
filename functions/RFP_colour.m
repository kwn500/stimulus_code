function [RFP_to_draw] = RFP_colour(distances, RGBmodulation)

% Loads in the RFP image file [500 x 500 double], and modulates this image
% with the RGB values reflecting the desired LMS cone excitations. 

% ----------------------------------------------------------------------- %
% ARGUEMENTS:
% distances- a [500 x 500] image matrix containing values (0-255) which 
% define an RFP. We rescale this so that they are the contrasts in a 
% colour space.

% RGB modulation- a previously computed [1 x 3] LMS cone excitation vector
% transformed to RGB space. This is multiplied by a corresponding scale to 
% specify the contrast of the image. 

% RETURNS:
% RFP_to_draw- a [500 x 500 x 3] image matrix containing values 
% specifying a RFP modulated in a specified LMS direction. 

% AUTHOR: KWN 
% DATE: 21/02/2018

% ----------------------------------------------------------------------- %
% get the size of the RFP image.
[iy, ix]=size(distances); 

if min(distances(:)) > -1 && max(distances(:)) < 1
    normDist = distances;
else
    % convert the colours from 0->255 to -1->1.
    normDist=(double(distances)-128)/128;
end


% the modulation vector of values ('RGBmodulation') is also normalised
% between -1 and 1.

% turn this normalised RFP image matrix into a big vector of numbers
% to make the multiplication better.
rgbImage=normDist(:);

% multiply this RFP image matrix by the LMS modulation vector (in RGB
% space), to create an image containing scaled RGB values which define
% a contrast in LMS space.
rgbImageFinal=rgbImage*RGBmodulation';

% reshape this image matrix back into an RGB image (y by x by 3
% planes).
RFP_to_draw=reshape(rgbImageFinal,[iy ix 3]);

% add 0.5 so that the midpoint is gray.
RFP_to_draw = RFP_to_draw+.5;

% this needs fixing, but at .0020 to anything currently set to 0.5 to 
% prevent the slight gray box outline between the stimulus and the background. 
% RFP_to_draw(RFP_to_draw == 0.5) = 0.5020;

RFP_to_draw(round(RFP_to_draw, 3) == 0.5) = 0.5020;
end


