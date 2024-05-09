function [newX, newY] = AddText(windowPtr, text, xyCenter, textColor, yJust, textSize)

% AddText.m
%    [newX newY] = AddText(windowPtr, text, xyCenter, textColor, yJust [,textSize])
%
% DESCRIPTION
%    Quickly add text to Screen.  'text' is centered horizontally, and can be
%    placed in one of three positions vertically (above, below, center, textSize).
%
% ARGUMENTS
%    'windowPtr' is the pointer for the Screen window to place the text.
%    'text' is a string of text to add to Screen
%    'xyCenter' is a 2-element vector denoting the center of the screen, 
%               [x y], in pixels.
%    'textColor' is the clut index (scalar or [r g b] triplet) that you 
%               want to use for the text.
%    'yJust' is the vertical justification as a proportion of TextHeight (derived).
%               0 for centered, + for above, - for below.
%    'textSize' is the size of the text to be drawn, passed through
%               Screen('textSize').  if left empty, AddText will simply
%               draw at whatever size is currently defined.
%
% RETURN
%   [newX newY] returns the final pen location.  This is the return of
%               Screen('DrawText'...)
%
% SEE ALSO
%    Screen('DrawText')

% 5/6/08   rehbm   Wrote it.
% 01.02.14 rehbm   Added support for changing textSize
% 01.29.14 rehbm   Make yJust more flexible.


%__________________________________________________________________________
% validate arguments
if ~ismember(nargin,[5 6])
    error('AddText: Usage: AddText(windowPtr, text, xyCenter, textColor, yJust, [textSize])')
end

% make sure windowPtr is valid
if ~find(Screen('Windows')==windowPtr)
    error('AddText: Invalid windowPointer.')
end

if ~ischar(text)
    error('AddText: text must be a string.')
end

if ~isvector(xyCenter) || length(xyCenter) ~= 2
    error('AddText: invalid format for xyCenter.  Must be [x y] or [x;y;].')
end

% yJust can not be any value
% if yJust~=-1 && yJust~=0 && yJust~=1
%     error('AddText: yJust must be 0 (centered), 1 (above) or -1 (below).')
% end

%__________________________________________________________________________
% do the work
if ~isempty(textSize)
    old_textSize = Screen('TextSize',windowPtr,textSize);
end

textBounds = Screen('TextBounds',windowPtr,text);
textWidth  = RectWidth(textBounds);
textHeight = RectHeight(textBounds);
x = xyCenter(1) - textWidth/2;
y = xyCenter(2) - yJust*textHeight;
[newX, newY] = Screen('DrawText',windowPtr,text,x,y,textColor);

% revert text size for future drawing
if ~isempty(textSize)
    Screen('TextSize',windowPtr,old_textSize);
end
