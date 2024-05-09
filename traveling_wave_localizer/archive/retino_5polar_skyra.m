%Stim Parameters for Polar Angle Mapping w/ attention dimming task @ stim
clear options

% Datafile --------------------------------
options.subj               =   'test';
options.notes              = '';

% Screen and Display --------------------------------
% change to match your scanner situation.
options.screen_w = 34;                          % screen width in cm. scanner 39, most lcds 50 
options.screen_d = 38.5;                          % screen viewing distance (cm)
options.bg_color = 0*[1 1 1];         % background color on 0-255 [R G B] scale


% Response --------------------------------
options.key.quit_key           =   KbName('Escape'); % Q for quit
options.key.resp_key           =   KbName('1!');  % response button, 44 for Keyboard, 33 for button box
[id names]  = GetKeyboardIndices;
dev_id = id(strcmp(names,'932')); % Trigger key code
options.key.dev_id         =   -1; % device id for checking responses.  -1 for all devices (but might be slow and miss some responses)


% Timing ----------------------------------------------------
options.ncycles            =   8;          % number of cycles
options.cycle_duration     =   40;         % duration (sec) of a single cycle (eg, rotation period of wedge or ring). *** should be an increment of the TR ***
options.init_blank_time    =   20;          % initial blank time (sec)
options.end_blank_time     =   0;          % blank time at end (sec)


% Fixation Point ----------------------------------------------------
options.fix.color          =   [255 0 0];     % color of fixation spot


% Stimulus ----------------------------------------------------
% stimuli are displayed full-screen.  'mask' defines how much and which part of the stimulus is visible at any point in time.
options.stimulus = 'checker'; % type of visual stimulus.  'checkerboard' or 'dots'


% Attention Task ----------------------------------------------------
options.attend_task         =   2;             % 1 for fixation dimming task, 2 for peripheral (ring/wedge/dots)
options.dim_value           =   [137 137 137]; % dim color (0-255 [R G B]) for checkerboard attention task
options.dim_length          =   0.15;          % duration (sec) of fixation and ring/wedge dimming events
options.min_targ            =   3;             % minimum time (sec) between target events
options.max_targ            =   5;             % maximum time (sec) between target events
options.attend_checker_ring  =   2;          % atten ring duty cycle (width) for wedge mask
    
    
% Mask ----------------------------------------------------
% General - shared across all masks
% N.B. the mask type is inherently defined by the direction
options.duty                = 0.125;   % duty cycle (0.125 = 12.5%): fraction of cycle time that stimulus is visible at any one point.  For wedge, angular size is 360 * duty.
options.direction           = {'ccw' 'cw' 'ccw' 'cw' 'ccw'}; % cell array of strings, one for each run.  options are wedge:'cw','ccw'     ring:'in','out'     bar:'bar'
options.transition          = 'smooth';  % 'smooth' for continuous updating of position and smoothly moving mask, 'discrete' for for static mask positions changes (see options.discrete.order)

% additional options for discrete transitions
options.discrete.order       = 'straight'; % defines the order of cycling through all possible cycle-step pairs.  'straight' (default) for straight progression, 'random' for random order, or vector containing values 1:(ncycles*steps_per_cycle) defining index of cycle-step pairs to move through
options.discrete.step_size   = 2; % (sec) if discrete.step_size is true, how long does the mask stay in one position?  in other words, at what time is the mask updated to step to the next position?  if 'calculate' the discrete.step_size will be calculated such that each point on the screen is part of one and only one step
options.discrete.step_offset = 'half_step_size'; % time offset (sec) to shift each discrete step by.  If 0, the first step will not show a simulus (or rather a sliver of one).  By default, it is 'half_step_size', which will be calculated as discrete.step_size/2.

% Rotating Wedge - standard polar angle mapping - ('cw' or 'ccw')
options.wedge.cut_out     =   0.5;        % cut out n deg (radius) from center. set to zero for no cut_out.
options.wedge.cut_in      =   13.5;         % cut in n deg (radius) from outer max; defines extent (max eccen) of wedge


% Blanks ----------------------------------------------------
% Use these parameters to create arbitrarily timed blank periods.
options.blanks = []; % Nx3 element matrix (e.g., [3 0.5 1; 4 0.33 0.66] where each row-triplet defines a blank period.  The 1st column is the cycle during which the blank occurs and the 2nd and 3rd columns defined when the blank starts and stops, respectively, in time_fraction units.  A full cycle goes from 0-1. In the example given, there would be no stimulus during the 2nd half of the 3rd cycle and the middle 3rd of the 4th cycle.


% Extras ----------------------------------------------------
options.debugging         = 0; % if true, will not HideCursor or ListenChar(2) and will keep track of the duration of each frame


TravelingWave(options);
ShowCursor;