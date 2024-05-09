function [newXY] = RotatePointsPTB(xy, ang, anchor)

% RotatePointsPTB.m
%    [newXY] = RotatePointsPTB(xy, ang, [,anchor])
%
% DESCRIPTION
%    Rotate the Cartesian coordinates (x,y) around the supplied 'anchor'
%    point by 'ang' (deg).  Degrees are in Psychtoolbox (PTB) format:
%    0 degrees is UP and positive ang go clockwise.
%
% ARGUMENTS
%    'xy' is an mx2 or 2xn matrix defining the Cartesian (x,y) coordinate
%         pairs. In this case 2x2 matrix this function will assume 
%         [x1 y1; x2 y2;] and produce a warning message on the screen.
%    'ang' is the angle, in degrees, to rotate points around anchor point.
%    'anchor' is an optional argument defining the point around which to 
%        rotate the supplied points.  Default is the origin (0,0).
%
% EXAMPLES
%    x = 10*rand(10,1);
%    y = 10*rand(10,1);
%    xy = cat(2,x,y);
%    [newX newY] = RotatePointsPTB(xy, 90, [1 -1]);
%          Rotate the 10 (x,y) points 90 degrees around (1,-1)
%    [newX newY] = RotatePointsPTB(xy, 270);
%          Rotate the 10 (x,y) points 270 degrees around the origin (0,0).
%
% RETURN
%    'newXY' is the mx2 or 2xm output list of [x y] or [x; y;] pairs. The 
%            output matrix will match the format of the input martix.
%
% SEE ALSO
%    Cart2PolPTB(), Pol2CartPTB()

% 5/1/08    rehbm   Wrote it.
% 5/20/08   rehbm   fixed help text.
% 01/02/14  rehbm   Added support for rotating PTB rects (1x4)


%__________________________________________________________________________
% validate arguments
if nargin < 2 || nargin > 3
    error('RotatePointsPTB: Usage RotatePoiuntsPTB(xy, ang, [,anchor]).')
elseif nargin ~= 3
    anchor = [0 0]; % default anchor is the origin
end

% validate format of xyDeg matrix
transposeOutput = 0;
isRect = 0;
[m n] = size(xy);

% special case for a PTB rect (1x4 vector)
if m==1 && n==4
    isRect = 1;
    xy = reshape(xy,2,2)'; % convert from PTB rect to [x1 y1; x2 y2]
    [m n] = size(xy); % update
elseif m~=2 && n~=2
    error('RotatePointsPTB: Invalid format for xy.  Must be mx2 or 2xn.')
elseif (m == 2) && (n == 2)
   % Then we can't really know what format the xy data is in (columns or rows)
   % We'll assume it is [x1 y1; x2 y2;] and print a warning
   display('_________________________________________________________________________');
   display('WARNING:');
   display('   RotatePointsPTB - input matrix was 2x2. Assuming [x1 y1; x2 y2;] format.');
elseif (m == 2)
    % Then this appears to be an xy array of the form
    %    [x1 x2 x3; y1 y2 y3;]
    xy = xy';
    transposeOutput = 1;
end

% check ang
if ~isscalar(ang)
   error('RotatePoiuntsPTB: ang must be scalar.')
end

% validate format of anchor matrix
if (length(anchor) ~= 2)
    error('Degrees2PixelsPTB : invalid anchor.')
elseif (size(anchor,1) == 2)
    % Then anchor is in the form [xAnchor; yAnchor;]
    anchor = anchor';
end



%__________________________________________________________________________
% do the work - account for anchor (translate), convert to polar, adjust
%    angle, convert back to cartesian, account for anchor (re-translate)
anchor = repmat(anchor,size(xy,1),1);
xy = xy - anchor;
[theta rho] = Cart2PolPTB(xy);
theta = theta + ang;
[newX newY] = Pol2CartPTB(theta, rho);
newXY = cat(2,newX,newY) + anchor;

% format for output
if transposeOutput
    newXY = newXY';
end
if isRect
    % convert back to a valid PTB rect (1x4), and make sure that the
    % lowest x/y values come first.  else, PTB can't draw it.
    x1 = min(newXY(:,1));
    x2 = max(newXY(:,1));
    y1 = min(newXY(:,2));
    y2 = max(newXY(:,2));
    newXY = [x1 y1 x2 y2];
end