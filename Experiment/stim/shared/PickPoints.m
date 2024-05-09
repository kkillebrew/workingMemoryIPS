function [newPoints] = PickPoints(nPoints, currentPoints, minDistanceBetween, boundry, nOverrides)
% pickPoints.m
%   [newPoints] = PickPoints(nPoints, currentPoints, minDistanceBetween, boundry [,nOverrides])
%
%   DESCRIPTION
%       Choose nPoints new points [x y] pairs at least minDistanceBetween from each other
%       given any [x y] currenPoints.  All points will bounded by boundry, which can either
%       be an [x y] array defining maximum X and maximum Y dimensions seperatly, or [r] scalar
%       defining the maximum radius.
%
%   ARGS
%       nPoints = number of desired points
%       currentPoints = current points as [x1 y1; x2 y2; ...] matrix
%       minDistanceBetween = [scalar] is minimum distance between any two points.
%                 or {min_dist_at_1_deg_ecc corticalMagFormula} to scale minDistanceBetween
%                 by cortical magnification factor. See ApplyCorticalMagFactor.m for forumla
%                 options.
%       boundry = 1. [maxX maxY] for central rectangle
%                 2. [maxRadius] for central circle
%                 3. [minEcc maxEcc minAng maxAng] for (partial) wedge,
%                                   whose center is defined by the boundry.
%                 boundry should be supplied in Degrees and in reference
%                 frame of Psychtoolbox (0 deg is up, positive angle
%                 clockwise).
%       nOverrides = optional argument specifing number of tries to select 
%               the points given the current parameters.  2 vector input.
%               First element is number of times to re-choose a set of remaining
%               point.  Second element is the number of times to start from
%               scratch (probably shouldn't be more than 2).  Default is
%               [100000 2].
%
%   RETURNS
%       nPoints x 2 matrix of new [x y] point pairs

%   04/08 REHBM
%   06.04.08 rehbm - added tries argument for better control over selection
%                    process
%   06.06.08 rehbm - added partial_wedge boundry criteria option
%   08.18.08 rehbm - fixed restarting so it works.  needed to keep track of
%                    original currentPoints and reset on each restart.
%   03.17.08 rehbm - added support for scaling minDistanceBetween by cortical mag factor,
%                  - cosmetic
%   03.26.09 rehbm - vectoried new point creation/validation for 3-4X speed increase!

if nargin < 4 || nargin > 5
    error('Usage [newPoints] = PickPoints(nPoints, currentPoints, minDistanceBetween, boundry [,nOverrides])');
end
if nargin == 4
   nOverrides = [100000 2]; % DEFAULT VALUE FOR nOverrides
end

%__________________________________________________________________________
% special cases
if nPoints == 0
    newPoints = [];
    return
end

%__________________________________________________________________________
% validate arguments
% nPoints must be a positive integer
if (nPoints < 1 || mod(nPoints,1))
    error('nPoints must be an positive integer (or zero)');
end
% -validate currentPoints structure
if ~isempty(currentPoints)
    % validate structure
    if size(currentPoints,2) ~= 2
        error('invalid sturcture for currentPoints.  Must be n x 2 matric of [x y] pairs (or empty)');
    end
end

% minDistanceBetween
use_cortical_mag = 0;
if iscell(minDistanceBetween) && length(minDistanceBetween)==2
    use_cortical_mag = 1;
    cortical_mag_base = minDistanceBetween{1};
    cortical_mag_formula = minDistanceBetween{2};
elseif ~isscalar(minDistanceBetween)
    error('minDistanceBetween must be a 2-element cell array or a scalar')
end

% validate boundry structure
badBoundry = 0;
switch length(boundry)
    case 1
        % boundry is scalar and defines max radius
        boundryType = 'radius';
    case 2
        % boundry is vector/matrix and defines maxX and maxY
        boundryType = 'xy';
        % validate structure
        [m n] = size(boundry);
        if ~(m == 1 && n == 2)
            badBoundry = 1;
        end
    case 4
        % boundry is vector and defines [minEcc maxEcc minAng maxAng]
        boundryType = 'partial_wedge';
        % validate structure
        [m n] = size(boundry);
        if ~(m == 1 && n == 4)
            badBoundry = 1;
        end
    otherwise
        badBoundry = 1;
end
if badBoundry
    error('invalid sturcture for boundry.  Must be scalar (for maxRadius), 1x2 array for [maxX maxY], or 1x4 array for [minEcc maxEcc minAng maxAng]');
end

if length(nOverrides)==2 || isvector(nOverrides)
    pointTries = nOverrides(1); % maximum times to try to find each point
    totalTries = nOverrides(2); % number of times to start from scratch
    
else
    error('invalid sturcture for nOverrides.  Must be a 2 element vector [triesPerPoint retries]');
end

%__________________________________________________________________________
% parameters
if ~use_cortical_mag
    % we only need to do this once
    minDist2 = minDistanceBetween * minDistanceBetween;
    minDistZero = minDistanceBetween==0;
else
    minDistZero = cortical_mag_base==0;
end
originalPoints = currentPoints;


%__________________________________________________________________________
% do the work

success = 0;
thisTry = 0;
while (thisTry<totalTries && ~success)
    %if thisTry > 0
    %   display(sprintf('RESTARTING %d',thisTry))
    %end

    thisTry = thisTry+1;
        
    totalFailed = 0;
    currentPoints = originalPoints; % in case we are restarting
    newPoints = zeros(nPoints,2);
    newPointsIdx = 1; % start by adding to first slot
    nNeeded = nPoints; % first pass, we need all the new points
    
%    for i = 1:nPoints
%        if totalFailed
%            % kick out to force re-start
%            break
%        end

        tries = 0;
        while (tries < pointTries)
            %%thisPoint = chooseRandomPoint(boundryType,boundry);
            thisPoints = chooseRandomPoints(nNeeded,boundryType,boundry);
            
            pointFailed = zeros(nNeeded,1); % initialize

            if ~minDistZero % skip all this intensive work if we don't care about distances
                
                % update minDist2 if using cortical_mag factor. otherwise, just use minDist2 from initilization
                if use_cortical_mag
                    thisEccs = sqrt(sum(thisPoints.^2,2));
                    minDist2 = ApplyCorticalMagFactor(thisEccs, cortical_mag_base, cortical_mag_formula).^2;
                end

                % check to see if ANY new point is too close to ANY OTHER new point
                d2New = sq_euclidean_distance(thisPoints,thisPoints);
                d2New(d2New==0) = Inf; % don't consider distance to self (along diagonal)
                newFailed = min(d2New)' < minDist2;% check distance^2 to the closest point
                pointFailed = or(pointFailed, newFailed);

                % check to see if ANY new point is too close to ANY existing point
                if ~isempty(currentPoints)
                    d2Current = sq_euclidean_distance(currentPoints,thisPoints);
                    d2Current(d2Current==0) = Inf; % don't consider distance to self (along diagonal)
                    currentFailed = min(d2Current)' < minDist2;% check distance^2 to the closest point
                    pointFailed = or(pointFailed, currentFailed);
                end
                
            end
            
            % keep valid points
            nValid = sum(~pointFailed);
            thisValidPoints = thisPoints(~pointFailed,:);
            newPoints(newPointsIdx:newPointsIdx+nValid-1,:) = thisValidPoints;
            newPointsIdx = newPointsIdx+nValid;
            currentPoints = cat(1,currentPoints,thisValidPoints); % update current points list for next cycle
            
            % continue to select more newPoints if necessary
            if ~any(pointFailed)
                % if we get this far then we have a completely valid set of new points in newPoints
                break
            else
                % try again: increment tries and abort if we haven't suceeded in pointTries tries
                tries = tries+1;
                if (tries == pointTries)
                    totalFailed = 1;
                    break
                end
                
                % update how many we need to select next cycle
                nNeeded = sum(pointFailed);
            end

        end % while
%    end % for nPoints
    
    
    if ~totalFailed
        success = 1;
    end
end % while totalTries

if totalFailed
    error('failed to find a valid set of points after %d tries given %d restarts',pointTries,totalTries);
end
end

%__________________________________________________________________________
% subfunctions
function [p] = chooseRandomPoints(n,boundryType, boundry)
switch boundryType
    case 'xy'
        % boundry is vector/matrix and defines maxX and maxY
        p = bsxfun(@times,boundry,2*(rand(n,2)-0.5));
    case 'radius'
        % boundry is scalar and defines max radius
        th = 2*pi*rand(n,1);       % random theta in radians
        rad = boundry*rand(n,1);   % random radius in units of boundry
        p = [rad.*cos(th) rad.*sin(th)];
    case 'partial_wedge'
        % boundry is vector and defines [minEcc maxEcc minAng maxAng] of a
        % partial wedge region
        
        % convert angle to matlab polar coordinates (0 rad right, positive
        % counter clockwise).
        min_th = deg2rad(-1*(boundry(3)-90));
        max_th = deg2rad(-1*(boundry(4)-90));
        
        th = min_th + (max_th-min_th) * rand(n,1);              % random theta in radians
        rad = boundry(1) + (boundry(2)-boundry(1)) * rand(n,1); % random radius in units of boundry
        p = [rad.*cos(th) rad.*sin(th)];
end
end

%% Euclidean distance
function [d] = sq_euclidean_distance(X,Y)
% The following algorithm was taken from Peter Acklam's
% "Matlab array manipulation tips and tricks"
% http://home.online.no/~pjacklam
% 11.1.4 Euclidean distance matrix 
%
% We return squared Euclidean distance to save time. We'll square minDist once instead,
%
% X = m-by-p matrix representing m points in p-dimensional space
% Y = n-by-p matrix representing another set of points in the same space.
% d(i,j) = the Euclidean distance between X(i,:) and Y(j,:)
m = size(X,1);
n = size(Y,1);
X = permute(X, [1 3 2]); 
Y = permute(Y, [3 1 2]); 
d = sum(abs( X(:,ones(1,n),:) ...
    - Y(ones(1,m),:,:) ).^2, 3);
end
