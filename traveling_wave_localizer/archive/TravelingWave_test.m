% ver 1, port of J Swisher's Vision Egg code (Swisher et al. 2007,JNeuro), ma 07/17/09
% modified (pure) log scaling of eccentricity, switched timing from frame rate to actual time, converted to function c said 09/20/09
% added cortical mag. function from Swisher et al. 2007 (ala Duncan and Boynton, 2002), ma 08/10/10
% added default value structure, use of paramters wrapper file, ability to do
% multiple runs, save behav data, ma 12/03/10
%
% 05.2011 - REHBM - added a secondary angle for target that is randomized after each target is turned off to avoid
%                   a consistency of the target angle with eccentricity (especially in eccentricity mapping)
%                 - various minor updates
%                 - added moving dot stimulus option

function [] = retinotopy(options)

try

    %clear mex;
    
    AssertOpenGL;
    rand('state',sum(100*clock));
    PsychJavaTrouble;
    Priority(9);
    
    
    % validate arguments
    if ~isfield(options,'direction')
        error('Please define at least one direction');
    end
    
    if ischar(options.direction)
        options.direction = {options.direction};
    end
    
    dir_check = ~ismember(options.direction,{'cw'; 'ccw'; 'in'; 'out'; 'bar'});
    if any(dir_check)
        error(['Invalid experiment type (options.direction) for run(s) ' mat2str(find(dir_check)) '. Must be ''cw'', ''ccw'', ''in'', ''out'', ''bar''.']);
    end
    
    
    
    %Changeable Paramters---------------------------------------------------
    
    %MAIN PARAMETERS--------------------------------
    %change to match your scanner situation.
    defaults.screen_w           =   39;         % screen width, 39 @ scanner (cm), 50 w/ LCD
    defaults.screen_d           =   60;         % screen viewing distance (cm)
    
    %Additional Parameters----------------------------------------------------
    % Stimulus & Run Timing
    defaults.rot_period         =   40;         % rotation period of wedge or ring (in seconds)
    defaults.ncycles            =   9;          % number of rotation cycles
    %NOTE: change rot period to be a increment of your TR
    defaults.init_blank_time    =   0;          % initial blank time (sec)
    defaults.end_blank_time     =   0;          % blank time at end (sec)
    
    %Stimulus Properties
    defaults.color_checker      =   1;          % 0 for B&W, 1 for color
    defaults.flicker_freq       =   4;          % flicker frequency for full black-white cycle (hz)
    defaults.flickstim          =   0;          % Flicker Stimulus 0 for checkerboard, 1 for black
    %Extra Parameter - keep to 0 if unsure
    defaults.checker_rot        =   0;          % rotation period of background checkboard (only); 0 for no rotation; pos for cw, neg for ccw
    
    % moving dots
    %...should go here...
    
    
    
    % Attention task parameters
    defaults.atten_task         =   2;          % 1 for fixation 2 for wedge/ring
    defaults.dim_value          =[137 137 137]; % dim color for attention task
    defaults.dim_length         =   0.08;        % amount dimmed for (in sec)
    defaults.min_dim            =   2;          % minimum time between fixation point dims
    defaults.max_dim            =   5;          % max time between fixation point dims
    %parameters for when atten_task == 2
    defaults.atten_ring         =   3;          % atten ring duty cycle (width) when running polar angle
    defaults.atten_wedge        =   30;         % degress of atten wedge size when runing eccentricity
    
    %Wedge Properties
    defaults.wedge_size         =   45;         % in degrees of a circle,
    defaults.cut_out            =   0;          % cut out n deg (radius) of center
    defaults.cut_in             =   15;         % cut in n deg (radius) from outer max; equals max eccen of wedge
    defaults.start_angle        =   90;         % right horizontal meridian (0 is top, going cw_out), where to start in degree
    
    %Ring Properties
    defaults.duty               =   0.125;       % Set duty cycle .125 = 12.5%
    defaults.r_outer_deg_max    =   15;         % Set outer max for ring
    defaults.r_inner_deg_min    =   0;          % Set inner min for ring
    defaults.scaling            =   3;          % 1 for linear, 2 for pure logarithmic, 3 log-based cortical magnification ala Boynton & Duncan, 2002
    
    %Key Codes
    defaults.key.resp_key       =   20;         % Q for quit
    defaults.resp_key           =   44;         % response button, 44 for Keyboard, 33 for button box
    
    %fixation point
    defaults.dim_color          =   [255 255 255]; % fixation color
    defaults.fp_size            =   4;          % radius in pixels, of fixation point
    
    
    %Additional Variables
    defaults.direction          =   {};
    defaults.subj               =   {};
    defaults.date               =   fix(clock);
    defaults.behav              =   [];
    defaults.savefile           =   1;
    behav.accuracy              =   [];
    
    defaults.debugging          = 0; % if true, will not HideCursor or ListenChar(2)
    
    % Set size of screen to display background images
    % Leave to max if unsure - Better to scale wedge and ring stimuli
    % directly
    % in code below
    stim_r = 'max';             % stimulus radius (deg), max for maximum screen size
    
    
    % merge input arguments (options) with defaults (defaults)
    options = propval(options,defaults);
    
    %End of Parameters-------------------------------------------------------------------
    
    
    % define variables for all of the option fields so that we can refer to them
    % without the option. strucutre.  just to keep the code below cleaner....
    fns = fieldnames(options);
    for i = 1:length(fns)
        fn = fns{i}; % string
        eval(sprintf('%s = options.%s;',fn,fn));
    end
    % angle of target for peripheral dimming task.  position is calculated (based on time) relative to this base angle.
    % base angle is reset after each target to keep target positions random relative to stimulus angle.  Technically,
    % this is only necessary for the ring stimulus since the wedge stimulus has additional randomness to target eccen
    targ_base_angle = 360*rand;
    
    %Create filename to save data: Subj_date_time
    filename = [num2str(options.subj) '_' date '_' num2str(options.date(4)) '_' num2str(options.date(5))];
    
    % Do not change for now -
    stim = 3;                                   % 1 for checkerboard  - Future editions will use other sitmuli
                                                % 2 scenes?
                                                % 3 random dot motion with coherent direction
                                                
    % Initialize Screen
    screens      = (Screen('Screens'));
    screenNumber = max(screens);                % should equal 1 if using my office Mac
    [w, rect]    = Screen('OpenWindow', screenNumber, 32);
    [keyIsDown, secs, keycode] = KbCheck(-1);   % check response
    Screen('FillRect',w, [0 150 150]);
    Screen('Flip', w,0,1);
    
    % Enable alpha blending for contrast manipulations (for drawing circular dots)
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    % Calculate frame rate
    frate = Screen('FrameRate',w);                % get frame rate from 'FrameRate'
    ifi = Screen('GetFlipInterval', w);
    if frate == 0                                 % if FrameRate returns 0....
        frate = 1/ifi;                            % get frame rate from GetFlipINterval
    end;
    
    white_index = WhiteIndex(w);                  % set white RGB values
    black_index = BlackIndex(w);                  % set black RGB values
    
    % Find center of screen, screen viewing, ring scaling, and wedge
    % starting position
    xc = rect(3)/2;                                         % calculate central point on x-axis
    yc = rect(4)/2;                                         % calculate central point on y-axis
    max_viewing = floor(sqrt((rect(4)^2)+(rect(3)^2)));     % maximum screen distance in radius (deg)
    self.nrp_C  = ((r_outer_deg_max))/(exp(2.5)-1);         % scaling constant for cortical mag.
    start_wedge = (start_angle*pi)/180;                     % where to start in radians
    
    % Calculate pixels per degree (ppd).
    ppd = rect(3)/(2*atan(screen_w/(screen_d*2))*180/pi);
    
    if strcmp(stim_r,'max')
        stim_r = max([xc yc])/ppd;
    end
    
    % Convert from ppd into pixels
    StimSize   = stim_r*2;
    stim_r     = stim_r*ppd;
    cut_out    = round(cut_out*ppd);
    cut_in     = round(cut_in*ppd);
    atten_ring = atten_ring*ppd;
    
    % Setup up stim timing variables
    total_time = init_blank_time + rot_period * ncycles + end_blank_time;
    flick_dur = 1/flicker_freq/2;
    count_dur = rot_period;
    
    % Read in background images
    if color_checker == 1
        checker1 = imread('./images/background0.png');
        checker2 = imread('./images/background1.png');
    elseif color_checker == 0
        checker1 = imread('./images/background_bw_0.png');
        checker2 = imread('./images/background_bw_1.png');
    end
    if stim == 2
        bkimage1  = imread('./images/Scene_06.jpg');
        bkimage2  = imread('./images/Scene_07.jpg');
        bkimage3  = imread('./images/Scene_11.jpg');
        bkimage4  = imread('./images/Scene_12.jpg');
        bkimage5  = imread('./images/Scene_14.jpg');
        bkimage6  = imread('./images/Scene_18.jpg');
        bkimage7  = imread('./images/Scene_25.jpg');
        bkimage8  = imread('./images/Scene_28.jpg');
        bkimage9  = imread('./images/Scene_35.jpg');
        bkimage10 = imread('./images/Scene_38.jpg');
        bkimage11 = imread('./images/Scene_40.jpg');
    end
    if flickstim
        checker1 = checker1./checker1;
    end
    
    % Resize images
    ImageSize = StimSize*ppd;
    
    imgscale = ImageSize/length(checker1);
    checker1 = imresize(checker1,imgscale);
    checker2 = imresize(checker2,imgscale);
    
    if stim == 2
        imgscale  = ImageSize/length(bkimage1);
        bkimage1  = imresize(bkimage1,imgscale);
        bkimage2  = imresize(bkimage2,imgscale);
        bkimage3  = imresize(bkimage3,imgscale);
        bkimage4  = imresize(bkimage4,imgscale);
        bkimage5  = imresize(bkimage5,imgscale);
        bkimage6  = imresize(bkimage6,imgscale);
        bkimage7  = imresize(bkimage7,imgscale);
        bkimage8  = imresize(bkimage8,imgscale);
        bkimage9  = imresize(bkimage9,imgscale);
        bkimage10 = imresize(bkimage10,imgscale);
        bkimage11 = imresize(bkimage11,imgscale);
        
    end
    
    checker_tex(1) = Screen('MakeTexture', w, checker1);
    checker_tex(2) = Screen('MakeTexture', w, checker2);
    
    
    %     if stim == 2
    %         bkimage_tex(1)  = Screen('MakeTexture', w, bkimage1);
    %         bkimage_tex(2)  = Screen('MakeTexture', w, bkimage2);
    %         bkimage_tex(3)  = Screen('MakeTexture', w, bkimage3);
    %         bkimage_tex(4)  = Screen('MakeTexture', w, bkimage4);
    %         bkimage_tex(5)  = Screen('MakeTexture', w, bkimage5);
    %         bkimage_tex(6)  = Screen('MakeTexture', w, bkimage6);
    %         bkimage_tex(7)  = Screen('MakeTexture', w, bkimage7);
    %         bkimage_tex(8)  = Screen('MakeTexture', w, bkimage8);
    %         bkimage_tex(9)  = Screen('MakeTexture', w, bkimage9);
    %         bkimage_tex(10) = Screen('MakeTexture', w, bkimage10);
    %         bkimage_tex(11) = Screen('MakeTexture', w, bkimage11);
    %
    %     end
    
    
    % setup for (potential) dot stimulus
    % make box texture for masking non-visible stimuli - textures are the easiest thing to place in an arbitrary coordinate with an arbitrary rotation
    corner_dist_pix = round(sqrt(2*max([xc yc])^2)); % maximum screen dimension along diagonal in pixels
    bar_mask       = repmat(0,256,256); % an image of a rectangle
    bar_mask_texid = Screen('MakeTexture',w,bar_mask,[],1);
    bar_mask_rect = [xc-corner_dist_pix yc-corner_dist_pix xc+corner_dist_pix yc+corner_dist_pix]; % central template for mask's drawing rect
    bar_angles_deg = [NaN 0 45 90]; % bars move perpendicular to the bar_angle vector.  0 is top-to-bottom, 90 is right-to_left, 45 is upperRight-to-lowerLeft, etc. use NaN for fiation-only (no stimulus, null) cycles *** SHOULD BE MOVED TO OPTIONS/DEFAULTS - needs to be linked to NCYCLES***
    

    if stim == 3
        p.norm_speed 	= 7;	% normal dot speed (deg/sec)
        p.ndots 		= 20000;	% number of dots
        %r_outer_deg_max 		= 15;	% maximum radius of  annulus (degrees)
        %p.min_d 		= 1;	% minumum
        %p.th_min = [0; pi];
        %p.th_max = [pi/4; 5*pi/4];
        %p.th0 = -pi/8;
        %p.th0 = 2.77;
        %p.n_prop = length(p.th_min);	% number of propellers
        %p.n_prop = 1;
        p.dotsize = 6;
        p.dotcolor = [255 255 255]; % 255*rand(3,p.ndots);
        
        %p.switch_coh_min = 1;
        %p.switch_coh_max = 1;
        p.fcoh = 1.0;%0.33;	% fraction of coherent dots
        p.fkill = 0.01; 	% fraction of dots to kill each frame
        
        
        % convert speed from deg/sec to pix/frame
        ds_to_pf = ppd / frate;	% conversion factor from deg/sec to pixels/frame
        pfs = p.norm_speed * ds_to_pf;	% normal dot speed (pixels/frame)
        
        
        % -----------------------------------
        % set up dot positions and velocities
        % -----------------------------------
        
        % define maximum radius of annulus (pixels from center)
        if stim == 3
            % bar - maximum dot eccentricity needs to extend to corner of screen (corner_dist_pix)
            rmax = corner_dist_pix;
        else
            % wedge/ring - maximum dot eccentricyt defined by r_outer_deg_max
            rmax = r_outer_deg_max * ppd;
        end
        rmin = 0; 	% minimum
    end
    
    
        
    % --------------------
    % start experiment now: draw fixation point and text and wait for key press to begin
    % --------------------
    
    for current_run = 1: length(options.direction);
        
        if ~debugging
            HideCursor;                         % Hide the mouse cursor
            ListenChar(2);
        end
        
        % Initialize variables that get reset at begining of each run
        fix_dimmed              = 0;
        behav.falarms           = 0;
        behav.hits              = 0;
        behav.total_dim         = 0;
        flicker_type            = 1;
        count                   = 1;
        timing.check_target     = 0;
        timing.check_prior      = 0;
        key.quit                = 0;
        checker_angle           = 0;
        stim_angle              = start_angle;
        targ_angle              = targ_base_angle; % selected randomly above
        timing.start_resp_time  = [];
        timing.resp_interval    = 1;  %Time alloted for correct response (in seconds)
        
        % Initialize Functions
        CheckButtonResponse(timing,key, behav);
        time_frac = 0;
        ring(r_outer_deg_max, r_inner_deg_min, time_frac, duty);
        
        %Display 'Wait for scanner to start' screen
        Screen('FillRect',w, black_index);
        Screen('FillRect',w, white_index, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]);   % set fixation point
        
        txt = 'Wait for experiment to start';                                                  % define text
        txtloc = [xc - length(txt) * 14 / 2, yc + 40];                                         % define text location
        Screen('TextSize', w, 25);                                                             % set text size
        [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),white_index);                % draw text
        
        Screen('Flip', w,0,1);
        
        
        % Select current stimulus -- curently only option is checkerboard
        switch stim
            case 1 % checker board stimulus
                current_stim = [checker_tex];
            case 2 % scene images
                current_stim = [shuffle(bkimage_tex),shuffle(bkimage_tex),shuffle(bkimage_tex),shuffle(bkimage_tex),shuffle(bkimage_tex),shuffle(bkimage_tex),shuffle(bkimage_tex)];
            case 3 % random dot motion
                % set up initial dot positions/motion directions
                % DOTPOS
                r = (rmax-rmin) * sqrt(rand(p.ndots,1)) + rmin;	% dot starting radius
                t = 2*pi*rand(p.ndots,1);	% dot starting theta
                xy = [r r] .* [cos(t), sin(t)];	% dot positions in Cartesian pixel coordinates relative to center
                
                mdir = 2*pi*rand(p.ndots,1);	% angular motion direction for each dot
                coh = rand(p.ndots,1) < p.fcoh;	% the coherently moving dots
                cohdir = 2*pi*rand; % defined randomly here.  for bar stim this will be updated depending on bar direction
                mdir(coh) = cohdir; 	% direction of coherent motion;
                dxdy = pfs * [cos(mdir), sin(mdir)];	% change in x and y (in pixels) per frame
                
                xdots = xy;%[xy ones(p.ndots,1) p.dotsize*ones(p.ndots,1)];				% [x, y, color, size] for each dot to be drawn

            otherwise
                error('invalid stim number provided')
        end
        
        %look for trigger, or button response to start experiment
        while ~keyIsDown
            [keyIsDown, secs, keycode] = KbCheck(-1); %check response
            if keycode(key.quit_key)
                key.quit = 1;
            end
        end
        
        Screen('FillRect',w,black_index);                                                       % reset screen to black
        Screen('FillRect',w, white_index, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]);    % set fixation point
        Screen('Flip', w);
        
        % Wait for initial blank time
        start_init_blank = GetSecs;
        while GetSecs - start_init_blank < init_blank_time
            WaitSecs(0.001); % avoid CPU hogging
        end
        
        Priority(MaxPriority(w));                                                               % reset to max priority
        
        %Timing Variables that get reset at begining of each run
        start_time = GetSecs;
        flicker_time = start_time+flick_dur;
        timing.next_fix = start_time + min_dim + (max_dim-min_dim).*rand;
        last_cycle = 0; % init for variables that change only during cycle change
        
        % Cycle Loop
        while(GetSecs-start_time < ncycles * rot_period) & ~key.quit
            
            %% Update Ring / Wedge
            cur_time = GetSecs - start_time; % current time relative to start of experiment
            time_frac = mod(cur_time, rot_period)/rot_period;           % fraction of cycle we've gone through
            cur_cycle = floor(cur_time/rot_period)+1;
            if any(strcmp(options.direction(current_run), {'ccw', 'in', 'bar'}))
                time_frac = 1-time_frac;                                            % reversing time for these conditions
            end
            
            
            % update variables that change depending on cycle, but are constant within a given cycle
            if last_cycle ~= cur_cycle
                % we've entered the next cycle
                if strcmp(direction{current_run},'bar') % BAR stim/mask
                    % update bar angle for next cycle
                    this_bar_angle_deg = bar_angles_deg(cur_cycle);

                    if stim == 3 && ~isnan(this_bar_angle_deg) % random dot motion (and bar stimulus)
                        % for dot motion with bar stimulus, direction always moves along bars main axis.
                        % update with each new cycle since bar direction will probably change (ignore for NaN too)
                        
                        % N.B. this_bar_angle_deg and Matlab's polar coordinate systems are offset by 90 deg.  but we
                        %      also need our dot motion direction offset by 90 deg from this_bar_angle_deg.  This is
                        %      inherent in the direct conversion of this_bar_angle_deg to radians (+ 50% chance of 180 deg flip)
                        cohdir = this_bar_angle_deg * pi/180 + round(rand)*pi;
                    end
                end
                last_cycle = cur_cycle; % udpate last_cycle so we only enter this once/cycle
            end
            
            
            switch direction{current_run}
                case {'cw' 'ccw' 'out' 'in'}
                    switch scaling
                        case 1 % linear scaling
                            [r_outer_deg r_inner_deg] = ring(r_outer_deg_max, r_inner_deg_min, time_frac, duty);

                        case 2 % log scaling
                            %First convert visual angle range to cortical distance range, only considering RF centers%
                            r_outer_ctx_max = log(r_outer_deg_max+1);                           % maximum stimulation distance from foveal representation in cortex, in mm, only considering RF centers
                            r_inner_ctx_min = log(r_inner_deg_min+1);
                            
                            %Then allow traveling wave to move linearly across cortex, only
                            %considering RF centers.
                            [r_outer_ctx r_inner_ctx] = ring(r_outer_ctx_max, r_inner_ctx_min, time_frac, duty);
                            
                            
                            %then convert back to visual angle. Duty cycle is preserved.%
                            r_outer_deg = exp(r_outer_ctx) - 1;
                            r_inner_deg = exp(r_inner_ctx) - 1;
                            
                        case 3 % cortical mag scalind
                            r_outer_ctx_max = (log((r_outer_deg_max/self.nrp_C)+1))/2.5;        % maximum stimulation distance from foveal representation in cortex, in mm, only considering RF centers
                            r_inner_ctx_min = (log((r_inner_deg_min/self.nrp_C)+1))/2.5;
                            
                            %Then allow traveling wave to move linearly across cortex, only
                            %considering RF centers.
                            [r_outer_ctx r_inner_ctx] = ring(r_outer_ctx_max, r_inner_ctx_min, time_frac, duty);
                            
                            %then convert back to visual angle. Duty cycle is preserved.%
                            r_outer_deg = (self.nrp_C * (exp((2.5*r_outer_ctx))-1));
                            r_inner_deg = (self.nrp_C * (exp((2.5*r_inner_ctx))-1));
                            
                        otherwise
                            error('Not a valid type of scaling')
                    end
                    
                case {'bar'}
                    if scaling ~= 1
                        error('must use linear (1) scaling with bar stimulus')
                    end
                    % i'm using ring function for now, even though variable names are setup for rings and wedges
                    [r_outer_deg r_inner_deg] = ring(r_outer_deg_max, r_inner_deg_min, time_frac, duty);
                    
                otherwise
                    error('invalid direction supplied')
            end
                    
            % now, regadless of scaling, convert visual angle to pixels
            r_outer_pix = r_outer_deg*ppd;
            r_inner_pix = r_inner_deg*ppd;
            
            % Calculate current angle of wedge (and potential target) based on run time
            stim_angle = (pi + (2*pi*time_frac - pi))*(180/pi)+start_angle;
            targ_angle = (pi + (2*pi*time_frac - pi))*(180/pi)+targ_base_angle;
            if checker_rot
                checker_angle = (pi+ (2*pi*mod(GetSecs-start_time, checker_rot)/checker_rot ))*(180/pi);
            end
            
            %% Check for flicker of stimulus; stim==2 is old code for use with scene/object stimuli. removed for now.
            if GetSecs > flicker_time
                flicker_time = flicker_time + flick_dur;
                %if stim ==1
                if flicker_type == 2
                    flicker_type = 1;
                else
                    flicker_type = 2;
                end
                %             elseif stim ==2
                %                 if flicker_type
                %                     flicker_type=0;
                %                 else
                %                     flicker_type=1;
                %                 end
                %             end
            end
            
            % Old code for different background stimuli. removed for now.
            % if stim ==2
            %     if GetSecs > next_switch
            %         next_switch = next_switch + (max_switch-min_switch).*rand;
            %         count=count+1;
            %     end
            % end
            
            
            %% Code for dimming both @ fixation or stimulus
            if ~fix_dimmed & GetSecs > timing.next_fix  % if not dimmed and time for next fixation dimming...
                % dim fixation cross
                if atten_task == 1
                    dim_color = dim_value;
                end
                fix_dimmed = 1;
                fix_start  = GetSecs;
                
                
                if (cut_in - cut_out) < atten_ring
                    dim_outer = cut_in;
                else
                    dim_outer = (cut_in-atten_ring-cut_out) * rand + (atten_ring+cut_out); % provides random target eccentricity for wedge stimulus
                end

                
                if stim == 3 && atten_task == 2
                    % for dot motion, choose new coherent motion direction (at least some amount away from current direction
                    % how direction changes depends on 'direction'
                    switch direction{current_run}
                        case {'cw' 'ccw' 'out' 'in'}
                            % WEDGE or RING
                            cohdir = cohdir + (sign(rand-.5) * pi/180 * (75+(105-75)*rand)); % +/- random number between 75 and 105 degrees (+ 50% sign flip)
                        case {'bar'}
                            % BAR
                            % direction always moves along bars main axis, and thus changes are always 180 deg
                            cohdir = cohdir + pi; % always 180 deg for bar mask
                        otherwise
                            error('cohdir change undefined for direction == %s',direction(current_run))
                    end
                    coh = rand(p.ndots,1) < p.fcoh;	% choose a new set of coherently moving dots
                    mdir(coh) = cohdir; % update motion directions for newly coherently moving dots (some dots will continue to move in old direction, but should be swamped by new direction)
                end
                                
                behav.total_dim = behav.total_dim+1;
                timing.start_resp_time = GetSecs;
                timing.check_target = 1;
            end
            
            % undim fixation cross
            if fix_dimmed & GetSecs > fix_start+dim_length % if fix dimmed and length of fix dim has passed (.35) then undim
                fix_dimmed = 0;
                if atten_task == 1
                    dim_color = [255 255 255]; % white
                end
                timing.next_fix = GetSecs + min_dim + (max_dim-min_dim).*rand;
                
                % reset targ_base_angle so target is shown at a random angle for next target presentation
                % provides random target angle for ring stimulus
                targ_base_angle = 360*rand;
            end
            
            % Check Button Responses
            [timing,key,behav] = CheckButtonResponse(timing,key,behav);
            
            
            %% update stimulus based on the current time
            % Draw Stimulus
            switch stim
                case 1 % checker board stimulus
                    Screen('DrawTexture',w,current_stim(flicker_type),[],[],checker_angle);
                    
                case 2 % scene stimuli
                    Screen('DrawTexture',w,current_stim(count));
                    if flicker_type
                        Screen('DrawTexture',w,checker_tex(1),[],[],[],[],flicker_type);
                    end
                    
                case 3 % random dot motionq
                    % update the dot positions for the next video frame
                    dxdy = pfs * [cos(mdir), sin(mdir)]; % change in x and y (in pixels) per frame
                    xy = xy + dxdy;						 % move dots
                    
                    % --------------------------------------------------------------------------------------------------
                    % ---these calculations are not for this frame, but the following frame (since xy is already set)---
                    % check to see which dots have gone beyond the borders of the annuli
                    r2 = xy(:,1).^2 + xy(:,2).^2;
                    r_out = r2 > rmax^2 | r2 < rmin^2;				% dots past outer border of central annulus
                    L = find(r_out | rand(p.ndots,1) < p.fkill);	% locus of all dots that have gone beyond borders or are killed (including a random selection that will always die each frame)
                    nL = length(L);
                    if nL
                        % choose new dot positions for killed dots
                        
                        r(L) = (rmax-rmin) * sqrt(rand(nL,1)) + rmin;	% dot starting radius
                        t(L) = 2*pi*rand(nL,1);	% dot starting theta
                        xy(L,:) = [r(L) r(L)] .* [cos(t(L)), sin(t(L))];	% dot positions in Cartesian pixel coordinates relative to center
                        
                        mdir(L) = 2*pi*rand(nL,1);	% motion direction for each dot
                        mdir(coh) = cohdir; % maintain coherence in motion
                        dxdy(L,:) = pfs * [cos(mdir(L)), sin(mdir(L))];
                    end
                    % --------------------------------------------------------------------------------------------------
                    
                    xdots = xy;	% new dot positions to be drawn
                    Screen('DrawDots', w, xdots', p.dotsize, p.dotcolor, [xc yc], 2);
                    
                otherwise
                    error('invalid stim number provided')
            end
            
            % draw the solid overlap that defines the stimulus shape (e.g., wedge, ring or bar)
            switch direction{current_run}
                case {'cw' 'ccw'}
                    % WEDGE SHAPE
                    if fix_dimmed && atten_task == 2 && ismember(stim,[1 2])
                        Screen('FrameArc',w, dim_value, [xc - dim_outer, yc - dim_outer, xc + dim_outer, yc + dim_outer], targ_angle-(.5*wedge_size), targ_angle+(.5*wedge_size), atten_ring, atten_ring );
                    end
                    Screen('FillArc', w, black_index, [xc-stim_r, yc-stim_r, xc+stim_r, yc+stim_r], stim_angle+(0.5*wedge_size), 360-wedge_size);
                    Screen('FrameArc',w, black_index, [xc-max_viewing, yc-max_viewing, xc+max_viewing, yc+max_viewing],0,360,max_viewing-(cut_in),max_viewing-(cut_in));
                    Screen('FillArc', w, black_index, [xc-(cut_out), yc-(cut_out), xc+(cut_out), yc+(cut_out)], 0, 360);
                    Screen('FillRect',w, dim_color, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]); % fixation spot
                    Screen('DrawingFinished', w);
            
                case {'out' 'in'}
                    % RING SHAPE
                    if fix_dimmed && atten_task == 2 && ismember(stim,[1 2])
                        Screen('FillArc', w, dim_value, [xc-stim_r, yc-stim_r, xc+stim_r, yc+stim_r], targ_angle+(.5*wedge_size), atten_wedge);
                    end
                    Screen('FrameArc', w, black_index, [xc-max_viewing, yc-max_viewing, xc+max_viewing, yc+max_viewing],0,360,max_viewing-r_outer_pix,max_viewing-r_outer_pix);
                    Screen('FillOval', w, black_index, [xc-r_inner_pix, yc-r_inner_pix, xc+r_inner_pix,yc+r_inner_pix]);
                    Screen('FillRect',w, dim_color, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]); % fixation spot
                    Screen('DrawingFinished', w);
                    
                case {'bar'}
                    % MOVING BAR STIMULUS
                    % "fix_dim" events taken care of above...

                    % mask everything beyond max_viewing (which essentially makes stimulus circular)
                    %Screen('FrameArc',w, black_index, [xc-max_viewing, yc-max_viewing, xc+max_viewing, yc+max_viewing],0,360,max_viewing-(cut_in),max_viewing-(cut_in));
                    
                    if isnan(this_bar_angle_deg)
                        % mask entire stimulus (fixation or null cycle)
                        %Screen('FillRect',w, black_index);
                    else
                        % partial stimulus mask that creates a moving bar
                        bar_theta = -1 * (-this_bar_angle_deg-90) * (pi/180); % we want 0 deg up (top-to-bottom direction), positive cw, Matlab is 0 right, positive ccw.  need to transform for matlab's pol2cart
                        
                        % upper mask
                        upper_ecc_pix = -corner_dist_pix-r_outer_pix;
                        [upper_x upper_y] = pol2cart(bar_theta,upper_ecc_pix);
                        upper_rect = OffsetRect(bar_mask_rect,upper_x,upper_y);
                        Screen('DrawTexture', w, bar_mask_texid, [], upper_rect, this_bar_angle_deg, 0);
                        
                        % lower mask
                        lower_ecc_pix = corner_dist_pix-r_inner_pix;
                        [lower_x lower_y] = pol2cart(bar_theta,lower_ecc_pix);
                        lower_rect = OffsetRect(bar_mask_rect,lower_x,lower_y);
                        Screen('DrawTexture', w, bar_mask_texid, [], lower_rect, this_bar_angle_deg, 0);
                    end
                    
                    Screen('FillRect',w, dim_color, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]); % fixation spot
                    Screen('DrawingFinished', w);%

                otherwise
                    error('invalid direction provided')
            end
            
            % Check Button Responses
            [timing,key,behav]= CheckButtonResponse(timing,key,behav);
            Screen('Flip', w);
            % Check Button Responses
            [timing,key,behav]= CheckButtonResponse(timing,key,behav);
            
        end
        
        
        % Put fix on screen for end fixation period
        Screen('FillRect',w, black_index);
        Screen('FillRect',w, white_index, [xc-fp_size; yc-fp_size; xc+fp_size; yc+fp_size]);
        Screen('Flip', w);
        
        % Wait for end blank time
        start_end_blank = GetSecs;
        while GetSecs - start_end_blank < end_blank_time
            WaitSecs(0.001); % avoid CPU hogging
        end
        
        
        %Calculate total runtime
        options.actual_runtime=GetSecs-start_time;
        behav.accuracy=behav.hits/behav.total_dim;
        options.behav=[options.behav;behav.falarms,behav.hits,behav.total_dim,behav.accuracy];
        
        %Report behav.accuracy after each run
        Screen('FillRect',w, [0 0 0]);
        txt1 = ['Your accuracy is: ' num2str(behav.accuracy)];
        txt2 = ['You had ' num2str(behav.falarms) ' false alarms'];
        txtloc = [xc - length(txt) * 7 / 2, yc];
        [newX newY] = Screen('DrawText',w,txt1,txtloc(1),txtloc(2)-15,white_index);
        [newX newY] = Screen('DrawText',w,txt2,txtloc(1),txtloc(2)+15,white_index);
        Screen('Flip', w,0,1);
        WaitSecs(5);
        
        key.quit=0; keyIsDown=0;                                            % reset keys

        % intermitent save
        if savefile
            save(fullfile('data',filename),'options');
        end
    end
    
    % final save
    if savefile
        save(fullfile('data',filename),'options');
    end
    
    ListenChar(0);
    Priority(0);
    Screen('CloseAll');
    ShowCursor;
    
    
catch
    % this "catch" section executes in case of an error in the "try" section
    % above.  Importantly, it closes the onscreen window if its open.
    
    % it is necessary to save the error message for this try...catch in
    % case some of the clean-up code includes a try...catch that might
    % overwrite the error message (i.e., ListenChar(0); seems to do this)
    thiserror = lasterror();
    
    % Clean up in case of an error in main code
    Screen('CloseAll'); % close any open screens/windows
    ShowCursor;         % restore cursor visibility
    ListenChar(0);      % keystrokes make it to command line/editor (Ctrl-c)
    Priority(0);        % restore normal processing priority
    
    % rethrow error, but print STACK first!
    %    "I'm the boss, need the info..." -Dr Evil
    display(sprintf('\n'));
    for i = 1:length(thiserror.stack)
        display(thiserror.stack(i));
    end
    rethrow(thiserror); % display error message that caused crash
    
    
end %try..catch..
end


function [r_outer r_inner] = ring(r_outer_max, r_inner_min, time_frac, duty) %This generic function can be used to generate linearly scaled rings on screen or on cortex.
r_width = duty*(r_outer_max-r_inner_min)/(1-duty); %fixed width of ring; (can obtain by solving for r_width in duty = r_width/(r_outer_max+r_width-r_inner_min)
r_outer_max_NOBLOCK = r_outer_max + r_width; %how far outer ring would go if there was no blocking ("blocking" is stopping the movement of the outer ring near the end of the cycle)

r_outer = r_inner_min + time_frac*((r_outer_max_NOBLOCK)-r_inner_min);% provisionally placing outer ring somewhere between r_inner_min and the maximum possible value if no blocking was done

r_inner = r_outer - r_width; %provisionally making inner ring a fixed distance inside outer ring;

r_outer = min(r_outer, r_outer_max); %... but blocking outer ring from going past its maximum allowed value.
r_inner = max(r_inner, r_inner_min); %...and not letting inner ring go below an assigned minimum (usually zero).
end



function [timing key behav] = CheckButtonResponse(timing,key, behav)
%% check for canceling, so easy to exit out for debugging

if timing.check_target %get response for last targ letter
    
    if GetSecs > timing.start_resp_time + timing.resp_interval
        timing.check_target=0;
    end
    
    [keyIsDown, secs, keycode] = KbCheck(-1); %check response
    
    if keycode(key.quit_key)
        key.quit=1;
    elseif find(keycode(key.resp_key))%if hit response key
        behav.hits = behav.hits + 1;
        timing.check_target=0;
        timing.check_prior=1;
        
    end
else %check for false alarms
    [keyIsDown, secs, keycode] = KbCheck(-1); %check response
    
    if keycode(key.quit_key)
        key.quit=1;
    elseif find(keycode(key.resp_key))%if hit response key
        if ~timing.check_prior
            behav.falarms=behav.falarms+1;
            timing.check_prior=1; %cause got a response
        end
    else
        %response bottom not down, will allow responses to be recorded agaiin
        timing.check_prior=0;
    end
    %__________________end response check
end
end


function [merged unused] = propval(propvals, defaults, varargin)

% Create a structure combining property-value pairs with default values.
%
% [MERGED UNUSED] = PROPVAL(PROPVALS, DEFAULTS, ...)
%
% Given a cell array or structure of property-value pairs
% (i.e. from VARARGIN or a structure of parameters), PROPVAL will
% merge the user specified values with those specified in the
% DEFAULTS structure and return the result in the structure
% MERGED.  Any user specified values that were not listed in
% DEFAULTS are output as property-value arguments in the cell array
% UNUSED.  STRICT is disabled in this mode.
%
% ALTERNATIVE USAGE:
%
% [ ARGS ] = PROPVAL(PROPVALS, DEFAULTS, ...)
%
% In this case, propval will assume that no user specified
% properties are meant to be "picked up" and STRICT mode will be enforced.
%
% ARGUMENTS:
%
% PROPVALS - Either a cell array of property-value pairs
%   (i.e. {'Property', Value, ...}) or a structure of equivalent form
%   (i.e. struct.Property = Value), to be merged with the values in
%   DEFAULTS.
%
% DEFAULTS - A structure where field names correspond to the
%   default value for any properties in PROPVALS.
%
% OPTIONAL ARGUMENTS:
%
% STRICT (default = true) - Use strict guidelines when processing
%   the property value pairs.  This will warn the user if an empty
%   DEFAULTS structure is passed in or if there are properties in
%   PROPVALS for which no default was provided.
%
% EXAMPLES:
%
% Simple function with two optional numerical parameters:
%
% function [result] = myfunc(data, varargin)
%
%   defaults.X = 5;
%   defaults.Y = 10;
%
%   args = propvals(varargin, defaults)
%
%   data = data * Y / X;
%
% >> myfunc(data)
%    This will run myfunc with X=5, Y=10 on the variable 'data'.
%
% >> myfunc(data, 'X', 0)
%    This will run myfunc with X=0, Y=10 (thus giving a
%    divide-by-zero error)
%
% >> myfunc(data, 'foo', 'bar') will run myfunc with X=5, Y=10, and
%    PROPVAL will give a warning that 'foo' has no default value,
%    since STRICT is true by default.
%

% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
%
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================

% Backwards compatibility
pvdef.ignore_missing_default = false;
pvdef.ignore_empty_defaults = false;

% check for the number of outputs
if nargout == 2
    pvdef.strict = false;
else
    pvdef.strict = true;
end

pvargs = pvdef;

% Recursively process the propval optional arguments (possible
% because we only recurse if optional parameters are given)
if ~isempty(varargin)
    pvargs = propval(varargin, pvdef);
end

% NOTE: Backwards compatibility with previous version of propval
if pvargs.ignore_missing_default | pvargs.ignore_empty_defaults
    pvargs.strict = false;
end

% check for a single cell argument; assume propvals is that argument
if iscell(propvals) && numel(propvals) == 1
    propvals = propvals{1};
end

% check for valid inputs
if ~iscell(propvals) & ~isstruct(propvals)
    error('Property-value pairs must be a cell array or a structure.');
end

if ~isstruct(defaults) & ~isempty(defaults)
    error('Defaults struct must be a structure.');
end

% check for empty defaults structure
if isempty(defaults)
    if pvargs.strict & ~pvargs.ignore_missing_default
        error('Empty defaults structure passed to propval.');
    end
    defaults = struct();
end

defaultnames = fieldnames(defaults);
defaultvalues = struct2cell(defaults);

% prepare the defaults structure, but also prepare casechecking
% structure with all case stripped
defaults = struct();
casecheck = struct();

for i = 1:numel(defaultnames)
    defaults.(defaultnames{i}) = defaultvalues{i};
    casecheck.(lower(defaultnames{i})) = defaultvalues{i};
end

% merged starts with the default values
merged = defaults;
unused = {};
used = struct();

properties = [];
values = [];

% To extract property value pairs, we use different methods
% depending on how they were passed in
if isstruct(propvals)
    properties = fieldnames(propvals);
    values = struct2cell(propvals);
else
    properties = { propvals{1:2:end} };
    values = { propvals{2:2:end} };
end

if numel(properties) ~= numel(values)
    error(sprintf('Found %g properties but only %g values.', numel(properties), ...
        numel(values)));
end

% merge new properties with defaults
for i = 1:numel(properties)
    
    if ~ischar(properties{i})
        error(sprintf('Property %g is not a string.', i));
    end
    
    % convert property names to lower case
    properties{i} = properties{i};
    
    % check for multiple usage
    if isfield(used, properties{i})
        error(sprintf('Property %s is defined more than once.\n', ...
            properties{i}));
    end
    
    % Check for case errors
    if isfield(casecheck, lower(properties{i})) & ...
            ~isfield(merged, properties{i})
        error(['Property ''%s'' is equal to a default property except ' ...
            'for case.'], properties{i});
    end
    
    % Merge with defaults
    if isfield(merged, properties{i})
        merged.(properties{i}) = values{i};
    else
        % add to unused property value pairs
        unused{end+1} = properties{i};
        unused{end+1} = values{i};
        
        % add to defaults, just in case, if the user isn't picking up "unused"
        if (nargout == 1 & ~pvargs.strict)
            merged.(properties{i}) = values{i};
        end
        
        if pvargs.strict
            error('Property ''%s'' has no default value.', properties{i});
        end
        
    end
    
    % mark as used
    used.(properties{i}) = true;
end
end
