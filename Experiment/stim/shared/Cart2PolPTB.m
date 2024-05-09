function [theta rho] = Cart2PolPTB(x, y)

% Cart2PolPTB.m
%    [theta rho] = Cart2PolPTB(x,y)
%    [theta rho] = Cart2PolPTB(xy)
%
% DESCRIPTION
%    Convert the Cartesian coordinates (x,y), where (0,0) is the center of 
%    the screen, +x is right and +y is up, to the polar coordinates used by 
%    Psychtoolbox (PTB).  Namely, 0 degrees is UP and positive thetas go 
%    clockwise.  Note that unlike cart2pol(), Cart2PolPTB returns theta in 
%    degrees, not radians.
%
% ARGUMENTS
%    'x' is the Cartesian x coordinate of the point in question. 'x' may
%        also be a vector. 
%    'y' is the Cartesian y coordinate of the point in question. 'y' may
%        also be a vector.
%    length(x) == length(y) for vector input.
%
%    Alternatively, if only a single argument is supplied, it is assumed to
%    supply (x,y) pairs in an mx2 or 2xn matrix. If the sole argument is a 
%    2x2 matrix this function will assume [x1 y1; x2 y2;] and produce a 
%    warning message on the screen.
%
% RETURN
%    'theta' is the CLOCKWISE angular displacement in DEGREES from the
%            POSITIVE Y-AXIS. Note that this is different than Matlab's
%            built-in cart2pol().    0<=theta<360
%    'rho'   is the distance from the origin to a point in the x-y plane.
%    
%    The output vectors 'theta' and 'rho' will be of the same format as the
%    input vectors 'x' and 'y' (i.e., column vectors or row vectors), but
%    will always be two different vectors (i.e., output will not be mx2 or
%    2xn matrix)
%
% SEE ALSO
%    cart2pol()

% 04/30/08   rehbm   Wrote it.


%__________________________________________________________________________
% validate arguments
if ~nargin || nargin > 3
    error('Cart2PolPTB: Usage Cart2PolPTB(x,y).')
end

% validate format of (x,y) point
transposeOutput = 0;
if nargin == 1
    % then check to see if we were given [x y] or [x; y;]
    [m n] = size(x);
    if (m ~= 2) && (n ~= 2)
        error('Cart2PolPTB: Invalid format for singular input.  Must be mx2 or 2xn matrix.')
    elseif (m == 2) && (n == 2)
        % Then we can't really know what format the xy data is in (columns or rows)
        % We'll assume it is [x1 y1; x2 y2;] and print a warning
%        display('_________________________________________________________________________');
%        display('WARNING:');
%        display('   Cart2PolPTB - singular input matrix was 2x2. Assuming [x1 y1; x2 y2;] format.');
    elseif (m == 2)
        % Then this appears to be an xy array of the form
        %    [x1 x2 x3; y1 y2 y3;]
        x = x';
        transposeOutput = 1;
    end
    
    y = x(:,2);
    x = x(:,1);
else
   % we were given two arguments, verify that they are two vectors of equal length
   if (~isvector(x) || ~isvector(y)) || (size(x,1) ~= size(y,1)) || (size(x,2) ~= size(y,2))
       error('Cart2PolPTB: Invalid input format. x and y must be vectors of the same dimensions.')
   end
end


%__________________________________________________________________________
% do the work - transform to polar, adjust thetale to match Psychtoolbox's
%   Screen funciton usage (degrees, 0 is UP, positive is clockwise)
[theta rho] = cart2pol(x,y);
theta = -theta * (180/pi) + 90; % transform to PTB format
theta = mod(theta,360);         % keep it 0<=theta<360

% format for output
if transposeOutput
    theta = theta';
    rho = rho';
end