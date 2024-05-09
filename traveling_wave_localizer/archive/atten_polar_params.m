%Stim Parameters for Polar Angle Mapping w/ attention dimming task @ stim
%clear options

%Stimulus Timing & Parameters
options.subj        =   ['test'];               % Subjects initials for saving behav data
options.direction   =   {'cw'};          % Direction(s); # of directions = # of runs, choices: cw, ccw in, out
options.rot_period  =   40;                     % Length of individual cycle (in seconds): default 40
options.ncycles     =   2;                      % # of cycles per run: default 8
options.wedge_size  =   45;                     % wedge size (in arc length): default 45
options.duty        =   .125;                   % duty cycle of ring: default .125 (each point on the screen is stimulated for 12.5% of cycle)
options.r_outer_deg_max = 20;                   % Set outer max for ring (in vis. degrees), in degrees:scanner default: 15
options.r_inner_deg_min = 0;                    % Set inner min for ring (in vis. degrees): default 0
options.scaling         = 1;                    % 1 for linear, 2 for logarithmic, 3 log-based cortical magnification ala Boynton & Duncan, 2002: default 3
options.init_blank_time = 0;                    % start of run blank time: default 20?
options.end_blank_time  = 4;                    % end of run blank time: default 0?
options.checker_rot     = 0;                    % rotation period of in background checkboard (only); 0 for no rotation; pos for cw, neg for ccw

%Attention Task Variables
options.atten_task  =   2;                      % 1 for fixation task, 2 for stimulus task: default 2
options.dim_length  =   1;%.05;                    % dim length (in seconds): can change to make task easier
options.dim_value   =   [127 127 127];          % dim color value (RGB)

options.min_dim            =   10;          % minimum time between fixation point dims
options.max_dim            =   10;          % max time between fixation point dims
    


%Button Response
options.key.quit_key    =   20;                 % key for quitting program (q)
%options.key.resp_key    =   33;                 % response key; 33 button box (4 button); 44 space
options.key.resp_key    =   30;                 % response key; index finger at Skyra is '1!' or 30

%change to match your scanner situation.
options.screen_w = 50;                          % screen width in cm. scanner 39, most lcds 50 
options.screen_d = 60;                          % screen viewing distance (cm)
options.savefile = 1;                           % 1 for create savefile

options.debugging = 1; % if true, will not HideCursor or ListenChar(2)

TravelingWave_test(options);
ShowCursor;