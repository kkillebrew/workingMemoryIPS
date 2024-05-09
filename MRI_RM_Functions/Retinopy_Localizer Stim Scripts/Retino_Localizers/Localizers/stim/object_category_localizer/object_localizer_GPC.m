try
    %get name of file
    
    %126 TRs each 2.5sec long = 255sec per run 4min 15sec per run
    % 25 axial slices 2x2x2 + 50% gap
    % TE = 40ms
    % Flip angle = 76º
    
    
    %initialzie random generate so get different trials
    
    clear all;
    clear mex;
    rand('state',sum(100*clock))
    filename = input('Input filename:', 's');
    
    addpath('./stimuli/'); % location of image files
    
    
    %conditions
    faces      = 1;
    houses     = 2;
    scrambled  = 3; % will use flowers and chairs for scrambled.
    genobjects = 4;
    fixation   = 10;

        epoch_order = [...
            10 1 2 3 4 10 2 4 1 3 10 4 3 2 1 10 3 1 4 2 10;
            10 3 1 4 2 10 4 3 2 1 10 2 4 1 3 10 1 2 3 4 10;
            ];

    
    epoch_times = repmat(NaN,size(epoch_order)); % to store ACTUAL time of block starts relative to run start
        
    [num_runs, num_epochs] = size(epoch_order);
    
    num_objcond = 4;
        
    %     resp_key=44;%50; %for space
    resp_key = KbName('1!'); % Skyra='1!'
    trigger_key = KbName('5%'); % Skyra='1!'

    quit_key = KbName('Escape');
    
    [id names]  = GetKeyboardIndices;
    dev_id = id(1); % Trigger key code
    %dev_id = -1; % -1 for all, but probably too slow to catch responses at Skyra
    
    
    resp_interval=.9; %time have to respond in secs
    first_targ=[2,6];%first target can be a repeat of stim2-6, second one a repeat of stim 9-12
    second_targ=[9,12];
    
    
    %set up stim timing
    
    block_length=15;
    
    total_time=block_length*num_epochs;
    
    stim_dur=.35;
    blank_dur=.40;
    
    white=WhiteIndex(0);
    black=BlackIndex(0);
    gray=(white+black)/2;
    Bkcolor=gray;
    lcolor=black;
    
    fcolor= [255 0 0];
    
    which_win = max(Screen('Screens'));
    [w, screenRect]=Screen('OpenWindow',which_win);
    [keyisdown, secs, keyCode] = KbCheck(dev_id);
    
    Xcenter=screenRect(3)/2;
    Ycenter=screenRect(4)/2;
    
    
    
    %initialize both buffers as black
    Screen('FillRect',w, Bkcolor);
    Screen('Flip', w,0,1);
    
    Screen('FillRect',w, Bkcolor);
    Screen('Flip', w,0,1);
    
    
    p.mon_width=34;%skyra: 34; my desktop: 50; %normally 40, screen width in cm, really scanner screen width 32 x height 24
    
    p.v_dist=38.5; % viewing distance in cm; skyra=38.5
    
    ppd = pi * screenRect(3) / atan(p.mon_width/p.v_dist/2) / 360;	% pixels per degree
    
    ppi=screenRect(3)/p.mon_width; %pixels per inch
    
    areappd=400;%12
    area=[Xcenter-areappd*.5, Ycenter-areappd*.5, Xcenter+areappd*.5, Ycenter+areappd*.5];
    offset=.5*ppd;%jitter for pics
    
    fsize=.15*ppd;
    fsq=[Xcenter - .5*fsize, Ycenter-.5*fsize, Xcenter+.5*fsize, Ycenter+.5*fsize];
    %     picsize=400;
    %
    %     osq=[Xcenter - .5*picsize, Ycenter-.5*picsize, Xcenter+.5*picsize, Ycenter+.5*picsize];
    Screen('TextSize',w, 26);
    
    %setup stimuli
    
    index=['01'; '02';'03'; '04'; '05'; '06'; '07'; '08'; '09'; '10'; '11'; '12'; '13'; '14'; '15'; '16'; '17'; '18'; '19'; '20'; ...
        '21'; '22'; '23'; '24'; '25'; '26'; '27'; '28'; '29'; '30'; '31'; '32'; '33'; '34'; '35'; '36'; '37'; '38'; '39'; '40'];
    numstim=40;
    
    %conditions
    %      faces=1;
    %      houses=2;
    %      scrambled=3;%will use flowers and chairs for scrambled.
    %      bodies=4;
    %      genobjects=5
    %      tools=6
    %      animals=7
    
    %body stim
%     name='Body_';
%     bodies_tex=zeros(40,1);
%     for i=1:numstim
%         bodies_tex(i)=Screen('MakeTexture', w, imread(strcat(name,index(i,:)), 'jpg'));
%     end
    
    %faces stim
    name='Faces_';
    faces_tex=zeros(40,1);
    for i=1:numstim
        faces_tex(i)=Screen('MakeTexture', w, imread(strcat(name,index(i,:)), 'jpg'));
    end
    
    %houses stim
    name='House_';
    houses_tex=zeros(40,1);
    for i=1:numstim
        houses_tex(i)=Screen('MakeTexture', w, imread(strcat(name,index(i,:)), 'jpg'));
    end
    
    %scrambles1 stim
    name='ScramGen_';
    scram1_tex=zeros(40,1);
    for i=1:numstim
        scram1_tex(i)=Screen('MakeTexture', w, imread(strcat(name,index(i,:)), 'jpg'));
    end
    
    %generic objects
    name='ObjectsGen_';
    genobjects_tex=zeros(40,1);
    for i=1:numstim
        genobjects_tex(i)=Screen('MakeTexture', w, imread(strcat(name,index(i,:)), 'jpg'));
    end
  
    
    % --------------------
    % start experiment now: draw fixation point and text and wait for key press to begin
    % --------------------
    
    ListenChar(2);
    HideCursor;
    quit=0;
    trial_num=0;
    
    data=[]; %data: trial_num, condition, corr or miss, rt
    for run=1:num_runs
        
        if ~quit
            tstim=zeros(40,num_objcond);
            
            
            for j=1:num_objcond
                
                tempstim=randperm(40)';
                %set targets
                %first targ, range 2-6
                r = ceil(5.*rand)+1;
                tempstim(r+1)=tempstim(r);
                %set secondtarg, range 9-12
                r = ceil(4.*rand)+8;
                tempstim(r+1)=tempstim(r);
                %first targ, range 22-26
                r = ceil(5.*rand)+21;
                tempstim(r+1)=tempstim(r);
                %set secondtarg, range 29-32
                r = ceil(4.*rand)+28;
                tempstim(r+1)=tempstim(r);
                
                tstim(:,j)=tempstim;
            end
            
            %now make one big array of all textures
            %conditions
            %      faces=1;
            %      houses=2;
            %      scrambled=3;%will use flowers and chairs for scrambled.
            %      genobj = 4;
            %      tstim=zeros(40,num_objcond,num_runs);
            counters=ones(1,num_objcond);
            
            stim=[];
  
            for i=1:num_epochs
                switch epoch_order(run,i)
                    case {1 2 3 4} % any image block
                        a=epoch_order(run,i);
                        stim=[stim; tstim(counters(a):counters(a)+19, a)];
                        counters(a)=mod(counters(a)+19,40)+1;
                    case 10 % fixblocks
                        % do nothing
                    otherwise
                        error('unrecognized epoch type.  all possible epochs need to be explicityly defined')
                end
            end
            
            
            %END STIM SETUP
            
            %opening screens
            
            
            Screen('FillRect',w, Bkcolor, []);
            Screen('FillRect',w,fcolor, [fsq(1), fsq(2), fsq(3), fsq(4)]);
            
            txt = 'Ready for next run?';
            txtloc = [Xcenter - length(txt) * 7 / 2, Ycenter + 40];
            [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),lcolor);
            
            if run>1%if not before first run, display behavior
                Screen('DrawText',w, num2str(run_behav), screenRect(3)-45, screenRect(2)+20, fcolor);
            end
            
            Screen('Flip', w,0,1);

            keyIsDown = 0;
            while 1
                [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                if keycode(resp_key)
                    break;
                end
            end
            
            FlushEvents('keydown');
            Screen('FillRect',w, Bkcolor,[]);
            
            Screen('FillRect',w,fcolor, [fsq(1), fsq(2), fsq(3), fsq(4)]);
            
            Screen('Flip', w,0,1);
            
            %waitsecs(1);
            
            txt = 'Wait for experiment to start';
            
            
            txtloc = [Xcenter - length(txt) * 7 / 2, Ycenter + 40];
            [newX newY] = Screen('DrawText',w,txt,txtloc(1),txtloc(2),lcolor);
            
            Screen('Flip', w,0,1);
            FlushEvents('keydown');
            
            keyIsDown = 0;
            while 1
                [keyIsDown, secs, keycode] = KbCheck(dev_id); % check response  ***check all devices (-1) to allow for manual triggering***
                if keycode(trigger_key)
                    break;
                end
            end
            
            % % %             %****WAIT OUT DUMMY SCANS FOR MANUAL TRIGGERING****
            % % %             Screen('FillRect',w,Bkcolor,[ ]);
            % % %             Screen('FillRect',w,fcolor, [fsq(1), fsq(2), fsq(3), fsq(4)]);
            % % %             [newX newY] = Screen('DrawText',w,'          ...',txtloc(1),txtloc(2),lcolor);
            % % %             Screen('Flip', w,0);
            % % %             WaitSecs(3*2.5); % 3 TRs for 2-3s TRs at the Skyra
            % % %             %**************************************************
            
            Screen('FillRect',w,Bkcolor,[ ]);
            Screen('FillRect',w,fcolor, [fsq(1), fsq(2), fsq(3), fsq(4)]);
            Screen('Flip', w,0,1);
            
            start_time=GetSecs;
            epoch_num=0;
            stimcount=1;
            
            
            %____________________Start new run
            
            while(GetSecs - start_time < total_time & ~quit)
                
                epoch_start_time=GetSecs;
                epoch_num=epoch_num+1;
                epoch_times(run,epoch_num) = epoch_start_time - start_time;
                
                if epoch_num <= num_epochs %as long as havve epochs left to do
                    
                    %____________start new epoch
                    
   
                    switch epoch_order(run,epoch_num)
                        case 1
                            tex=faces_tex;
                        case 2
                            tex=houses_tex;
                        case 3
                            tex=scram1_tex;
                        case 4
                            tex=genobjects_tex;
                      
                    end
                    %variables reset for each block
                    check=0;
                    flipcounter=epoch_start_time;
                    t=0;
                    
                    if epoch_order(run,epoch_num)==fixation
                        firsttime=1;
                        while GetSecs - epoch_start_time < block_length & ~quit
                            if firsttime
                                Screen('FillRect',w,Bkcolor);
                                Screen('FillRect',w,fcolor,[fsq(1), fsq(2), fsq(3), fsq(4)]);
                                Screen('Flip', w,0,1);
                                firsttime=0;
                            end
                            
                            %_______________check for responses
                            
                            
                            [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                            if keycode(quit_key)
                                quit=1;
                            end
                            
                            %__________________end response check
                        end
                    else %every object presentation condition
                        while GetSecs - epoch_start_time < block_length & ~quit
                            temp=randperm(4);
                            while t==temp(1)
                                temp=randperm(4);
                            end
                            t=temp(1);
                            
                            switch t
                                case 1%offset to right
                                    offsets=[offset, 0, offset, 0];
                                case 2%left
                                    offsets=[-offset, 0, -offset, 0];
                                case 3 %up
                                    offsets=[0,-offset, 0, -offset];
                                case 4%down
                                    offsets=[0,offset, 0 offset];
                            end
                            
                            
                            Screen('DrawTexture', w, tex(stim(stimcount)), [], [(area(1)+offsets(1)), (area(2)+offsets(2)), (area(3)+offsets(3)), (area(4)+offsets(4))]);
                            Screen('FillRect',w,fcolor,[fsq(1), fsq(2), fsq(3), fsq(4)]);
                            Screen('Flip', w,0,1);
                            
                            if stimcount > 1 & stim(stimcount)==stim(stimcount-1)%if repeating image
                                check=1; %check for responses
                                start_resp_time=GetSecs;
                                %new trial so increment and save data
                                trial_num=trial_num+1;
                                data=[data; trial_num, epoch_order(run, epoch_num),run, 0, 0];
                                %data: trial_num,condition, run num, corr or miss,respose time
                            end
                            stimcount=stimcount+1;
                            
                            while GetSecs - flipcounter < stim_dur & ~quit
                                %_______________check for responses
                                if check %get response for last targ letter
                                    
                                    if GetSecs > start_resp_time + resp_interval
                                        check=0;
                                    end
                                    
                                    [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                                    if keycode(quit_key)
                                        quit=1;
                                    end
                                    if find(keycode(resp_key))%if hit response key
                                        data(trial_num, 4)=1;%1 for hit response
                                        data(trial_num, 5)=GetSecs-start_resp_time;
                                        check=0;
                                    end
                                end
                                %__________________end response check
                            end
                            Screen('FillRect',w,Bkcolor);
                            Screen('FillRect',w,fcolor,[fsq(1), fsq(2), fsq(3), fsq(4)]);
                            Screen('Flip', w,0,1);
                            
                            while GetSecs - flipcounter < stim_dur + blank_dur & ~quit
                                %_______________check for responses
                                if check %get response for last targ letter
                                    
                                    if GetSecs > start_resp_time + resp_interval
                                        check=0;
                                    end
                                    
                                    [keyIsDown, secs, keycode] = KbCheck(dev_id); %check response
                                    if keycode(quit_key)
                                        quit=1;
                                    end
                                    if find(keycode(resp_key))%if hit response key
                                        data(trial_num, 4)=1;%1 for hit response
                                        data(trial_num, 5)=GetSecs-start_resp_time;
                                        check=0;
                                    end
                                end
                                %__________________end response check
                            end
                            flipcounter=flipcounter+stim_dur+blank_dur;
                        end%while time left in presentation block
                    end%if which kind of epoch
                end%if epochs left to do
            end%time left in runs
            actual_run_time=GetSecs-start_time
            
            save(fullfile('./data/',filename)); % workspace

            
            %calculate accuracy to display on screen
            temp=[];
            temp=data(:, 3)==run; %all rows of run just completed
            temp=find(temp);%get indexes of rows
            Ltemp=length(temp); %how many trials of type 1
            rt=0;
            num_corr=0;
            if Ltemp>0 % if have some trials
                for k=1:Ltemp
                    if data(temp(k), 4)==1 %got correct
                        rt=rt+ data(temp(k), 5);
                        num_corr=num_corr+1;
                    end
                end
                run_behav=(num_corr/Ltemp)*100;
                run_rt=rt/Ltemp;
            end
        end%if not quiting
        
    end%for loop for runs
    
    avg_data=[];    %data: trial_num, condition,run num, corr or miss, rt
    total_num_corr=0;
    for cond=1:8
        temp=[];
        temp=data(:, 2)==cond; %all rows of 1st cond
        temp=find(temp);%get indexes of rows
        Ltemp=length(temp); %how many trials of type 1
        rt=0;
        num_corr=0;
        if Ltemp>0 % if have some trials
            for k=1:Ltemp
                if data(temp(k), 4)==1 %got correct
                    rt=rt+ data(temp(k), 5);
                    num_corr=num_corr+1;
                    total_num_corr=total_num_corr+1;
                end
            end
            
            
            if num_corr >0
                avg_data=[avg_data; cond, num_corr/Ltemp, rt/Ltemp];
            else
                avg_data=[avg_data; cond, 0,0];
            end
        end
    end
    
    avg_data
    
    [Ttemp, y]=size(data);
    avg_total=total_num_corr/Ttemp

    save(fullfile('./data/',filename)); % workspace
    
    ListenChar(0);
    ShowCursor;
    
    Screen('CloseAll');
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.

    % it is necessary to save the error message for this try...catch in
    % case some of the clean-up code includes a try...catch that might
    % overwrite the error message (i.e., ListenChar(0); seems to do this)
    thiserror = lasterror();

    % save workspace!  otherwise we'll potentially loose all access to the data
    save('./crash_dump.mat');
    fprintf('\n#----------------------------------#\nsaving workspace in ./crash_dump.mat\n#----------------------------------#\n')
    
    % Clean up in case of an error in main code
    Screen('CloseAll'); % close any open screens/windows
    ShowCursor;         % restore cursor visibility
    ListenChar(1);      % keystrokes make it to command line/editor (Ctrl-c)
    Priority(0);        % restore normal processing priority

    % rethrow error, but print STACK first!
    %    "I'm the boss, need the info..." -Dr Evil
    display(sprintf('\n'));
    for i = 1:length(thiserror.stack)
        display(thiserror.stack(i));
    end
    rethrow(thiserror); % display error message that caused crash
end %try..catch..