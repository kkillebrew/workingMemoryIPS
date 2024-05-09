function [x y] = Pol2CartPTB(theta, rho)

% Pol2CartPTB.m
%    [x y] = Pol2CartPTB(theta,rho)
%    [x y] = Pol2CartPTB(thetarho)
%
% DESCRIPTION
%    Convert the polar coordinates (theta,rho) used by Psychtoolbox (PTB)
%    to the Cartesian coordinates where (0,0) is the center of the screen, 
%    +x is right and +y is up.  Namely, 0 degrees is UP and positive 
%    angles go clockwise.  Note that unlike pol2cart(), Pol2CartPTB expects
%    theta in degrees, not radians.
%
% ARGUMENTS
%    'theta' is the CLOCKWISE angular displacement in DEGREES from the
%            POSITIVE Y-AXIS. Note that this is different than Matlab's
%            built-in pol2cart().    0<=theta<360
%    'rho'   is the distance from the origin to a point in the x-y plane
%
%    length(theta) == length(rho) for vector input.
%
% RETURN
%    'x' is the Cartesian x coordinate of the point in question. 'x' may
%        also be a vector. 
%    'y' is the Cartesian y coordinate of the point in question. 'y' may
%        also be a vector.
%
%    The output vectors 'x' and 'y' will be of the same format as the
%    input vectors 'theta' and 'x' (i.e., column vectors or row vectors), 
%    but will always be two different vectors (i.e., output will not be mx2
%    or 2xn matrix)
%
% SEE ALSO
%    pol2cart()

% 5/1/08   rehbm   Wrote it.


%__________________________________________________________________________
% parameters
% if input theta == 0,90,180,or 270, then x or y should be zero, but will
% come back as a very small value (e.g., -2.4493e-014).  So we'll check for
% extremely small values and make them zero.  An EXTREMELY conservative
% lower bound for the size of each pixel in degrees is 1.0e-10.
e = 1.0e-10;

%__________________________________________________________________________
% validate arguments
if ~nargin || nargin > 2
    error('Pol2CartPTB: Usage Pol2CartPTB(theta,rho).')
end

% validate format of (x,y) point
transposeOutput = 0;
if nargin == 1
    % then check to see if we were given [theta rho] or [theta; rho;]
    [m n] = size(theta);
    if (m ~= 2) && (n ~= 2)
        error('Pol2CartPTB: Invalid format for singular input.  Must be mx2 or 2xn matrix.')
    elseif (m == 2) && (n == 2)
        % Then we can't really know what format the theta_rho data is in (columns or rows)
        % We'll assume it is [theta1 rho1; theta2 rho2;] and print a warning
        display('_________________________________________________________________________');
        display('WARNING:');
        display('   Pol2CartPTB - singular input matrix was 2x2. Assuming [theta1 rho1; theta2 rho2;] format.');
    elseif (m == 2)
        % Then this appears to be an theta_rho array of the form
        %    [theta1 theta2 theta3; rho1 rho2 rho3;]
        theta = theta';
        transposeOutput = 1;
    end
    
    rho = theta(:,2);
    theta = theta(:,1);
else
   % we were given two arguments, verify that they are two vectors of equal length
   if (~isvector(theta) || ~isvector(rho)) || (size(theta,1) ~= size(rho,1)) || (size(theta,2) ~= size(rho,2))
       error('Pol2CartPTB: Invalid input format. theta and rho must be vectors of the same dimensions.')
   end
end


%__________________________________________________________________________
% do the work - adjust angle to match Matlab's pol2cart() usage (radians, 
%    0 is RIGHT, positive is counter-clockwise), transform to Cartesian
theta = -1 * (theta-90) * (pi/180); % transorm to Matlab format
[x y] = pol2cart(theta,rho);

% check for extremely small values of x or y.  round to zero.
x = x .* (abs(x)>e);
y = y .* (abs(y)>e);

% format for output
if transposeOutput
    x = x';
    y = y';
end