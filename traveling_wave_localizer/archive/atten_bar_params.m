%Stim Parameters for Polar Angle Mapping w/ attention dimming task @ stim
clear p

%Stimulus Timing & Parameters
p.subj        =   ['test'];               % Subjects initials for saving behav data
p.direction   =   {'bar'};          % Direction(s); # of directions = # of runs, choices: cw, ccw in, out
p.cycle_duration  =   10;                     % Length of individual cycle (in seconds): default 40
p.ncycles     =   8;                      % # of cycles per run: default 8
%p.wedge_size  =   45;                     % wedge size (in arc length): default 45
p.duty        =   .125;                   % duty cycle of ring or bar: default .125 (each point on the screen is stimulated for 12.5% of cycle)
p.ring.r_outer_deg_max = 20;                   % Set outer max for bar (in vis. degrees), in degrees:scanner default: 15
p.ring.r_inner_deg_min = 0; % Set inner (lower) min for bar (in vis. degrees): default 0

p.wedge.cut_out         =   1;          % cut out n deg (radius) of center (central circular mask, defines minimum eccentricity of wedge)
p.wedge.cut_in          =   15;         % cut in n deg (radius) from outer max (defines maximum eccen of wedge)

%p.ring.scaling         = 1;                    % 1 for linear, 2 for logarithmic, 3 log-based cortical magnification ala Boynton & Duncan, 2002: default 3
p.init_blank_time = 0;                    % start of run blank time: default 20?
p.end_blank_time  = 4;                    % end of run blank time: default 0?
p.checker.checker_rot     = 0;                    % rotation period of in background checkboard (only); 0 for no rotation; pos for cw, neg for ccw

%Attention Task Variables
p.atten_task  =   2;                      % 1 for fixation task, 2 for stimulus task: default 2
p.dim_length  =   1;%.05;                    % dim length (in seconds): can change to make task easier
p.dim_value   =   [127 127 127];          % dim color value (RGB)

p.min_targ            =   1.9;          % minimum time between fixation point dims
p.max_targ            =   2.1;          % max time between fixation point dims
    


%Button Response
p.key.quit_key    =   20;                 % key for quitting program (q)
%p.key.resp_key    =   33;                 % response key; 33 button box (4 button); 44 space
p.key.resp_key    =   30;                 % response key; index finger at Skyra is '1!' or 30

%change to match your scanner situation.
p.screen_w = 50;                          % screen width in cm. scanner 39, most lcds 50 
p.screen_d = 60;                          % screen viewing distance (cm)
p.savefile = 1;                           % 1 for create savefile

p.debugging = 1; % if true, will not HideCursor or ListenChar(2)


TravelingWave_test4(p);
ShowCursor;