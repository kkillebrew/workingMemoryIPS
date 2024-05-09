function [xyPix] = Degrees2PixelsPTB(xyDeg, xyCenter, pixelsPerDegree)

% Degrees2PixelsPTB.m
%    [xyPix] = Degrees2PixelsPTB(xyDeg, xyCenter, pixelsPerDegree)
%
% DESCRIPTION
%    Convert an xy point array from visual degree/experimental coordinate 
%    system to the pixel/screen coordinate system for use with Psychtoolbox
%    functions.  Basically, we convert visual degrees to pixels, transform
%    (0,0) center into (0,0) upper left, and invert y-axis (positive is up
%    in experimental coordinates and down for Psychtoolbox screen).
%
% ARGUMENTS
%    'xyDeg' is the m x 2 or 2 x n input list of [x y] or [x; y;] pairs in
%            the visual degrees/experimental coordinate system. If xyDeg is
%            2x2, this function will assume [x1 y1; x2 y2;] and produce a
%            warning message on the screen.
%
%            Alternatively, if xyDeg is a 1x4 array, then this function
%            will assume that this is in the Psychtoolbox form of a
%            rectangle coordinates, [xLeft yTop xRight yBottom].  In this
%            case, the output will be in the same format.
%    'xyCenter' is a 2 element array denoting the center of the screen in
%            pixels. [xCenter yCenter] or [xCenter; yCenter;]
%    'pixelsPerDegree' is a scalar denoting the number of pixels per degree
%            of visual angle for the current display setup.
%
% RETURN
%    'xyPix' is the m x 2 or 2 x m output list of [x y] or [x; y;] pairs in
%            the pixel/screen coordinate system.  If xyDeg is 1x4 array,
%            then xyPix will also be 1x4 array.

% 4/28/08   rehbm   Wrote it.
% 8/19/08   rehbm   Small fix for when we need to transpose xyDeg.  need to
%                   recalculate [m n] size matrix.


%__________________________________________________________________________
% validate arguments
if nargin ~= 3
    error('Degrees2PixelsPTB: Usage Degrees2PixelsPTB(xyDeg, xyCenter, pixelsPerDegree).')
end

% validate format of xyDeg matrix
transposeOutput = 0;
rectangleOutput = 0;
[m n] = size(xyDeg);
if (m~=2 && n~=2) && (m~=1 && n~=4)
    error('Degrees2PixelsPTB: Invalid format for xyDeg.  Must be mx2, 2xn, or 1x4 matrix.')
end

if (m == 2) && (n == 2)
   % Then we can't really know what format the xy data is in (columns or rows)
   % We'll assume it is [x1 y1; x2 y2;] and print a warning
   display('_________________________________________________________________________');
   display('WARNING:');
   display('   Degrees2PixelsPTB - input matrix was 2x2. Assuming [x1 y1; x2 y2;] format.');
elseif (m == 1) && (n == 4)
    % Then this appears to be a rectangle coordinate of the form
    %    [xLeft yTop xRight yBottom]
    xyDeg = reshape(xyDeg,2,2)';    % [xLeft yTop; xRight yBottom;]
    [m n] = size(xyDeg);
    rectangleOutput = 1;
elseif (m == 2)
    % Then this appears to be an xy array of the form
    %    [x1 x2 x3; y1 y2 y3;]
    xyDeg = xyDeg';
    [m n] = size(xyDeg);
    transposeOutput = 1;
end

% validate format of xyCenter matrix
if (length(xyCenter) ~= 2)
    error('Degrees2PixelsPTB : invalid xyCenter array.')
elseif (size(xyCenter,1) == 2)
    % Then xyCenter is in the form [xCenter; yCenter;]
    xyCenter = xyCenter';
end


%__________________________________________________________________________
% do the work - scale, invert y-axis and center
xyPix = (xyDeg .* pixelsPerDegree) .* (repmat([1 -1],m,1)) + (repmat(xyCenter,m,1));
% format for output
if transposeOutput
    xyPix = xyPix';
end
if rectangleOutput
   xyPix = xyPix([1 3 2 4]);
end