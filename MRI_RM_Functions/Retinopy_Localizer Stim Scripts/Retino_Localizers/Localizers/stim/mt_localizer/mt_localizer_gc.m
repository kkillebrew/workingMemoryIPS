
% dot motion demo using Screen('DrawDots') subfunction
% author: Keith Schneider, 12/13/04


AssertOpenGL;

try
    
    %clear everything, set randomization
    clear all;     
    clear mex; 
    rand('state',sum(100*clock));
    [id names]  = GetKeyboardIndices;
    dev_id = id(1); % Trigger key code
    resp_key = KbName('1!');
    trigger_key = KbName('5%');
    %dev_id = -1; % -1 for all, but probably too slow to catch responses at Skyra
    
    %_______
    %set up experimental variables
    %__________________
    
    moving=1;
    static=2;
    
    %order of conditions, each line is a run.
    epoch_order=[2,1,2,1,2,1,2,1,2,1,2
                 2,1,2,1,2,1,2,1,2,1,2];
    
    [num_runs, num_epochs]=size(epoch_order);
    epoch_length=15;%should be multiple of TR
    
    total_time=epoch_length*num_epochs;
    
    stim_dur=1;%how long each direction lasts.
    
    
    %-----------
    %fixation task variables
    %--------------
    min_fix=2; %minimum time between fixatino point dims
    max_fix=5; %max time between fixation point dims
    dim_value=[75 0 0];%amount dimmed to from white
    dim_length=.090;%amount dimmed for in s
    fix_r       = 0.15; % radius of fixation point (deg)
    undim_value=[255 0 0];%color of fixation point when not dimmed
    fcolor=undim_value;
    % ------------------------
    % set dot field parameters
    % ------------------------

%     nframes     = 1000; % number of animation frames in loop
    
    mon_width=34;%skyra: 34; my desktop: 50; %normally 40, screen width in cm, really scanner screen width 32 x height 24
    v_dist=38.5; % viewing distance in cm; skyra=38.5
    
    dot_speed   = 8;    % dot speed (deg/sec)
    ndots       = 1000; % number of dots
    max_d       = 15;   % maximum radius of  annulus (degrees)
    min_d       = .75;    % minumum
    dot_w       = 0.15;  % width of dot (deg)
    f_kill      = 0.05; % fraction of dots to kill each frame (limited lifetime)    
    differentcolors =0; % Use a different color for each point if == 1. Use common color white if == 0.
    differentsizes = 0; % Use different sizes for each point if >= 1. Use one common size if == 0.
    waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.
    
    quit_key=KbName('Escape');
    

    
    if differentsizes>0  % drawing large dots is a bit slower
        ndots=round(ndots/5);
    end
    
    % ---------------
    % open the screen
    % ---------------

    doublebuffer=1
    screens=Screen('Screens');
	screenNumber=max(screens);
    % [w, rect] = Screen('OpenWindow', screenNumber, 0,[1,1,801,601],[], doublebuffer+1);
    [w, rect] = Screen('OpenWindow', screenNumber, 0,[], 32, doublebuffer+1);

    % Enable alpha blending with proper blend-function. We need it
    % for drawing of smoothed points:
%    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    [center(1), center(2)] = RectCenter(rect);
	 fps=Screen('FrameRate',w);      % frames per second
    ifi=Screen('GetFlipInterval', w);
    if fps==0
       fps=1/ifi;
    end;
    
    black = BlackIndex(w);
    white = WhiteIndex(w);
    gray=white/2;
    
    
    
    %%%%%%set up some variables for dots
    
    ppd = pi * (rect(3)-rect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
    pfs = dot_speed * ppd / fps;                            % dot speed (pixels/frame)
    s = dot_w * ppd;                                        % dot size (pixels)
    fix_cord = [center-fix_r*ppd center+fix_r*ppd];

    rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
    rmin = min_d * ppd; % minimum

        
%    Screen('FillRect', w, black)
    HideCursor;	% Hide the mouse cursor
    Priority(MaxPriority(w));
    ListenChar(2);
    
    % Do initial flip...
    vbl=Screen('Flip', w, black);
    
    
    for run=1:num_runs
        
        
        % ---------------------------------------
        % initialize dot positions and velocities
        % ---------------------------------------

        
        r = rmax * sqrt(rand(ndots,1));	% r
        r(r<rmin) = rmin;
        t = 2*pi*rand(ndots,1);                     % theta polar coordinate
        cs = [cos(t), sin(t)];
        xy = [r r] .* cs;   % dot positions in Cartesian coordinates (pixels from center)

        mdir = ones(ndots,1);    % motion direction (in or out) for each dot
        dr = pfs * mdir;                            % change in radius per frame (pixels)
        dxdy = [dr dr].* cs;                       % change in x and y per frame (pixels)
        xymatrix = transpose(xy);
        % Create a vector with different colors for each single dot, if
        % requested:
        if (differentcolors==1)
            colvect = uint8(round(rand(3,ndots)*255));
        else
            colvect=white;
        end;

        % Create a vector with different point sizes for each single dot, if
        % requested:
        if (differentsizes>0)
            s=(1+rand(1, ndots)*(differentsizes-1))*s;        
        end;

       
        
        %opening screens
            
        Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
%         Screen('FillOval', w, [128 128 128], fix_cord);	% draw fixation dot 
        
        txt = 'Press Button when ready to begin';
        Screen('TextSize', w, 25);                                                             % set text size
        txtloc = [center(1) - length(txt) * 7 / 2, center(2) + 40];
        [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),white);
        
        if run>1%if not before first run, display behavior
            txt = sprintf('  total_dim = %s',num2str(total_dim));    
            txtloc = [center(1) - length(txt) * 7 / 2, center(2) + 70];
                [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),white);
         end

        Screen('Flip', w); 

        %look for trigger, or button response
        %keyIsDown = 0;
        while 1
             [~, secs, keycode] = KbCheck(dev_id); % check response  ***check all devices to allow for manual triggering***
             if keycode(resp_key)
                 KbReleaseWait
                break;
             end
            WaitSecs(.001);
        end
    
        Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
%         Screen('FillOval', w, [128 128 128], fix_cord);	% draw fixation dot 
        
        txt = 'Wait for experiment to start';

        txtloc = [center(1) - length(txt) * 7 / 2, center(2) + 40];
        [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),white);
        Screen('Flip', w); 

        %look for trigger, or button response
        keyIsDown = 0;
        while 1
             [keyIsDown, secs, keycode] = KbCheck(dev_id); % check response  ***check all devices to allow for manual triggering***
             if keycode(trigger_key)
                break;
             end
            WaitSecs(.001);
        end
        
        %****WAIT OUT DUMMY SCANS FOR MANUAL TRIGGERING****
%         Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
%         [newX newY] = Screen('DrawText',w,'                      ...',txtloc(1),txtloc(2),white);
%         Screen('Flip', w,0);
%         WaitSecs(3*2.5); % 3 TRs for 2-3s TRs at the Skyra
        %**************************************************
        
        % trigger detected (and dummy scans waited)
        Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
        Screen('DrawDots', w, xymatrix, s, colvect, center,1);  % change 1 to 0 to draw square dots
        Screen('Flip',w,0,1);
        
        
        start_time=GetSecs;
        epoch_num=0;
        
        next_fix = start_time+ min_fix + (max_fix-min_fix).*rand;
        fix_dimmed=0;
        total_dim=0;
        fcolor=undim_value;
        quit=0;
        % --------------
        % animation loop
        % --------------    
        while(GetSecs - start_time < total_time & ~quit)
            
            epoch_start_time=GetSecs;
            epoch_num=epoch_num+1;
            %redraw fix at beginng of each block
%             Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
%             Screen('Flip',w);
        
            if epoch_num <= num_epochs %as long as havve epochs left to do

                %____________start new epoch    

                if epoch_order(run, epoch_num)==static
                    fix=1;
                else
                    fix=0;
                end
                
                firsttime=1;
                next_switch=epoch_start_time+stim_dur;
                
                
                while(GetSecs-epoch_start_time < epoch_length & ~quit)
                    if firsttime
                       Screen('DrawDots', w, xymatrix, s, colvect, center,1);  % change 1 to 0 to draw square dots
                       Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
                       Screen('Flip',w,0,1);  
                       firsttime=0;
                    end;
                    if fix
                       if ~fix_dimmed & GetSecs > next_fix % if not dimmed and time for nex fixation dimming
                            %dimm fixation cross
                            fcolor=dim_value;
                            fix_dimmed=1;
                            fix_start=GetSecs;
                            total_dim=total_dim+1;
                            Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
                             Screen('Flip',w,0,1);
        
                        end

                        %undim fixation cross
                        if fix_dimmed & GetSecs > fix_start+dim_length % if fix dimmed and length of fix dim has passed (.35) then undim
                            fix_dimmed=0;
                            fcolor=undim_value;
                            next_fix = GetSecs+ min_fix + (max_fix-min_fix).*rand;
                            Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot 
                            Screen('Flip',w,0,1);
        
                        end 
                        [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                         if keycode(quit_key)
                            quit=1;
                         end
                         WaitSecs(.001);
                    else%if moving epoch
                        if ~firsttime
                            Screen('FillOval', w, fcolor, fix_cord);	% draw fixation dot (flip erases it)

                            Screen('DrawDots', w, xymatrix, s, colvect, center,1);  % change 1 to 0 to draw square dots
                            Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                        end;
                        firsttime=0;

                        if GetSecs > next_switch
                            
                            if mdir(1) == 1    %swich motion dir
                                mdir=-1*ones(ndots,1);
                            else
                                mdir=ones(ndots,1);
                            end
                            %recalcuate
                            dr = pfs * mdir;                            % change in radius per frame (pixels)
                            dxdy = [dr dr] .* cs;  
                            next_switch=GetSecs+stim_dur;
                        end
                        
                            
                        xy = xy + dxdy;						% move dots
                        r = r + dr;							% update polar coordinates too

                        % check to see which dots have gone beyond the borders of the annuli

                        r_out = find(r > rmax | r < rmin | rand(ndots,1) < f_kill);	% dots to reposition
                        nout = length(r_out);

                        if nout

                            % choose new coordinates

                            r(r_out) = rmax * sqrt(rand(nout,1));
                            r(r<rmin) = rmin;
                            t(r_out) = 2*pi*(rand(nout,1));

                            % now convert the polar coordinates to Cartesian

                            cs(r_out,:) = [cos(t(r_out)), sin(t(r_out))];
                            xy(r_out,:) = [r(r_out) r(r_out)] .* cs(r_out,:);

                            % compute the new cartesian velocities

                            dxdy(r_out,:) = [dr(r_out) dr(r_out)] .* cs(r_out,:);
                        end;
                        xymatrix = transpose(xy);
                        if ~fix_dimmed & GetSecs > next_fix % if not dimmed and time for nex fixation dimming
                            %dimm fixation cross
                            fcolor=dim_value;
                            fix_dimmed=1;
                            fix_start=GetSecs;
                            total_dim=total_dim+1;
        
                        end

                        %undim fixation cross
                        if fix_dimmed & GetSecs > fix_start+dim_length % if fix dimmed and length of fix dim has passed (.35) then undim
                            fix_dimmed=0;
                            fcolor=undim_value;
                            next_fix = GetSecs+ min_fix + (max_fix-min_fix).*rand;
        
                        end 
                        [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                         if keycode(quit_key)
                            quit=1;
                         end
                        
                        if (doublebuffer==1)
                            vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
                         end;
                         %pause(0.001);
                         %pause;
                    end;%if fix
                end;%while time left in epoch
                
            end;%if have epochs left
        end;%while time left in epoch
        total_run_time=GetSecs-start_time
    end;%for each run
    
    Priority(0);
    ShowCursor
    ListenChar(0);
    Screen('CloseAll');
catch
    Priority(0);
    ShowCursor
    ListenChar(0);
    Screen('CloseAll');
    rethrow(lasterror);
end
