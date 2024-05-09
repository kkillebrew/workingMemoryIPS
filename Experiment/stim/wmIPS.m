% script-based expeirmental code for a working memory task based off of
% Sheremata et al. 2010 [JNeuro 30(38):12581-12588]

%  Changes KWK made:
%  -  changed params.setup to include screen_num for multiple displays
%  -  updated keys for instructions to include all keys instead of just
%  response keys to match instructions
%  -  updated instuctions to say "Lef/Right Shift" and "Peripheral cue =
%  memory task on one side or either"
%  -  added the params.kb_dev_id to all the KbWaits etc. to ensure correct
%  keyboard responses
%  -  line 435 temp_feat_idx changed to randomly select orientation so as
%  to not have all same oreintation
%  -  changed spacing between stimuli for letters and font size of letters

% Questions/concerns KWK:
%  -  what does the program do if no key is pressed?
%  -  what is the alloted time for key press?
%  -  stimulus looks really small on the stim comp...does it scale with
%  monitor size/resolution?
%  -  need to make mem stim not all vertical/horizontal
%  -  how is numtrials calculated?
%  -  how is the number of blocks determined?  -- determined by memory reps
%  -  where do you keep track of whether the trial was "new" or "old"
%  -  should the trials be random? and should there be a pre stim cue
%  before each trial?
%  -  where is the variable that keeps track of what type of trial is

%  -  passive vs. memory is in blocks.task
%  -  trial type (aka left right both) trials.mem_side
%  -  new/old trials -- trials.targ_change

% line

% to do:

% open questions / parameters
% - *** should the distractors ever change features from sample to probe?
% - *** should reported performance exclude passive-viewing blocks?  else
%       performance is artificially enhanced (since passive-viewing always
%       correct)
% - *** need to adjust monitor params for behavioral (and then fMRI) setup
%         marian will check with gideon (57 cm viewing distnace, 15" macbook pro 2012 model)
% - *** is X vs. O taping into verbal WM, or visual?
% - *** what is the scanner trigger key?
% - *** do you want the ability to counterbalance targ/dist color, and/or
%       response mapping across subjects? or even within subjects (split by
%       runs?)
% - *** pre and postfix time need to be adjusted for the scanner (at least
%       12-15 s)
% - *** should we explicitly avoid trial with only one target feature for
%       all targets?

% - timing of task seems like the pressure is put on encoding, not memory
%         this is standard, but behavioral data will tell us whether to
%         extend 
% - size and extent of stimulus placement
% - general timing and size paramters


% perhaps a parametric modulation from visual-spatial sketchpad to
% phonological loop using different sets of letters?  X/O -> F/R/K
% could combine expected strategy or look at self-reported strategey 


clear params blocks trials timing data this tmp stim cue fix; % clear some variables/structures that are important
AssertOpenGL; % make sure Psychtoolbox-3 is installed properly on this system
rand('twister',sum(100.*clock)); % shuffle the random number generator

addpath('./shared/'); % for shared helper functions

%% parameters
% all necessary parameters should be stored in the params structure for
% organization and easy params.datafile saving/loading.
params.experiment = 'wm'; % brief string to describe experiment.  will be part of datafile.  ideally short and only characters (no spaces, underscores, dashes, etc.). 

params.setup = 'kyle_lab_comp';
switch params.setup
    case 'ryan_macair'
        % monitor dimensions
        params.mon_width_cm = 28.5; % height 18 cm
        params.mon_dist_cm = 66;
        params.mon_width_deg = 2*(180/pi)*atan((params.mon_width_cm/2)/params.mon_dist_cm);
        
        % when using multiple monitors set the monitor number
        params.screen_num = 0;
        
        % window size (leave as empty matrix [] for full screen)
        params.win_size = []; %[0 0 960 600];[0 0 480 300];

        % what keyboard to listen to
        [kb_id, kb_names] = GetKeyboardIndices;
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device
        
        % Size of letters needs to be set here due to effects of screen
        % size
        params.stim.letters.letter_size = 70; % size of letters through TextSize

        debugging = 1;
    case 'kyle_lab_comp'
        % *** this is where we could add details for a different setup
        % monitor dimensions
        params.mon_width_cm = 40; % height X cm
        params.mon_dist_cm = 73;
        params.mon_width_deg = 2*(180/pi)*atan((params.mon_width_cm/2)/params.mon_dist_cm);
        
        % when using multiple monitors set the monitor number
        params.screen_num = 0;
        
        % window size (leave as empty matrix [] for full screen)
        params.win_size = [];

        % what keyboard to listen to
        [kb_id, kb_names] = GetKeyboardIndices;
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Wireless Keyboard')); % use -1 for all devices, or specify a specific device
        
        % Size of letters needs to be set here due to effects of screen
        % size
        params.stim.letters.letter_size = 70; % size of letters through TextSize

        debugging = 0;
    case 'stim_comp'
        % monitor dimensions
        params.mon_width_cm = 40; % height 18 cm
        params.mon_dist_cm = 73;
        params.mon_width_deg = 2*(180/pi)*atan((params.mon_width_cm/2)/params.mon_dist_cm);
        
        % when using multiple monitors set the monitor number
        params.screen_num = 1;
        
        % window size (leave as empty matrix [] for full screen)
        params.win_size = []; %[0 0 960 600];[0 0 480 300];

        % what keyboard to listen to
        [kb_id, kb_names] = GetKeyboardIndices;
        params.kb_dev_id = kb_id(strcmp(kb_names,'Apple Keyboard')); % use -1 for all devices, or specify a specific device
        
        % Size of letters needs to be set here due to effects of screen
        % size
        params.stim.letters.letter_size = 45; % size of letters through TextSize
        
        debugging = 0;
    otherwise
        error('unrecognized setup %s',params.setup)
end

% params that determine number of trials
% num trials per block = set_size*trials_per_block + npassive
% 7 blocks per run 6 trials per block (block cons
params.blocks.npassive = 1; % number of passive viewing blocks per run
params.blocks.memory_reps = 2; % number of reps for each unique memory block (determined by a cross of other parameters, like set size and task-type)

params.blocks.set_size = 3; % list of all possible set sizes (equally drawn from)
params.blocks.mem_side = [-1 1 2]; % which side could the target stimuli appear? -1, 1, 2 {left, right, either} [DO NOT INCLUDE PASSIVE HERE SINCE IT WILL HAVE A DIFFERENT NUMBER OF REPS]

params.trials_per_block = 6; % how many trials make up a single block?


% fixation spot
params.fix.innerwidth_deg = 0.1; % width of inner portion of fixation spot in degrees
params.fix.inner_color = [0 0 0]; % RGB value (0-255) of inner fix spot
params.fix.outerwidth_deg = 0.25; % width of outer portion of fixation spot in degrees
params.fix.outer_color = [255 255 255]; % RGB value (0-255) of outer fix spot

% block cue
params.cue.width_deg = 0.5; % width of block cue
params.cue.color = [255 255 255]; % RGB value (0-255) of block cue
params.cue.eccen = 2; % eccentricity in deg of block cue

% misc
params.instruction_text_size = 24; % size of instructions on screen
params.text_color = [0 0 0]; % RGB value (0-255) of instructional text
params.bg_color   = [128 128 128]; % RGB value (0-255) of background


% response keys
KbName('UnifyKeyNames'); % cross-platform defined key names
params.abort_key  = KbName('Escape'); % pressing escape kills the code
params.return_key = KbName('Return'); % continuing after a break, or starting experiment, requires pressing the Return key
params.trigger_key = KbName('5!'); % key that serves as scanner trigger
params.resp_keys  = KbName({'LeftShift' 'RightShift'}); % subject response key(s) index for KbCheck, etc.


% stimulus timing
params.prefix_dur  = 3; % duration of prefixation time (before trial 1) in seconds
params.postfix_dur = 3; % duration of final fixation period (after last trial) in seconds
params.block_cue_dur       = 2; % duration of block-cue (cue on ONLY) in seconds
params.block_cue_blank_dur = 2; % duration of block-cue blank (AFTER cue) in seconds
params.sample_dur = 0.100; % sample duration in seconds
params.delay_dur  = 0.900; % delay duration in seconds
params.probe_dur  = 2; % probe duration in seconds
params.iti_dur    = [3]; % possible inter-trial-interval durations in seconds (randomly, but evenly, sampled)

% stimulus
params.prop_targ_change = 0.50; % probability that any given trial will contain a target change
params.stim_per_hemifield = 8; % how many stimuli are presented in each hemifield on each trial?
params.hemifield_width_deg = 3; % width of each stimulus field (defining possible stim positions) in deg
params.hemifield_height_deg = 4; % heigh of each stimulus field (defining possible stim positions) in deg
params.hemifield_eccen_deg = 5; % how many degrees should each "field" of stimuli be offset to each side?
params.stimtype = 'letters_orientation'; % type of stimulus (currently, for entire run)
% below, you MUST define params.stim.(params.stimtype).feature_values as a
% row-wise list of possible feature values to select from.  also define
% any other params that are necessary for drawing the stimulus.
switch params.stimtype
    case 'orientation'
        params.stim.orientation.feature_values = [0 90]'; % orientation of bars (0-360, 0 up- vertical, + clockwise)
        params.stim.orientation.targ_color = [255 0 0]; % RGB value (0-255) of target stim
        params.stim.orientation.dist_color = [0 0 255]; % RGB value (0-255) of distractor stim
        params.stim.orientation.bar_length_deg = 0.7; % deg
        params.stim.orientation.bar_width_deg = 0.23; % deg
        
        params.stim_min_dist = 1.1 * params.stim.orientation.bar_length_deg; % min distance between two stimuli within the same hemifield
    
    case 'color'
        params.stim.color.feature_values = [255 0 0; 0 255 0; 0 0 255]; % orientation of bars (0-360, 0 up- vertical, + clockwise)
        params.stim.color.targ_shape = 'circle'; % 'circle' or 'square'
        params.stim.color.dist_shape = 'square'; % 'circle' or 'square'
        params.stim.color.width_deg = 0.7; % width/height of each element
        
        params.stim_min_dist = 1.1 * params.stim.color.width_deg; % min distance between two stimuli within the same hemifield
        
    case 'squares'
        params.stim.squares.feature_values = [0 45]'; % orientation of squares (0-square, 45 diamond)
        params.stim.squares.targ_color = [255 0 0]; % RGB value (0-255) of target stim
        params.stim.squares.dist_color = [0 0 255]; % RGB value (0-255) of distractor stim
        params.stim.squares.width_deg = 0.7; % deg
        
        params.stim_min_dist = 1.1 * sqrt(2*params.stim.squares.width_deg^2); % min distance between two stimuli within the same hemifield
    
    case 'letters'
        tmp_letters = 'RK';%'XO';
        params.stim.letters.feature_values = transpose(tmp_letters); % orientation of bars (0-360, 0 up- vertical, + clockwise)
        params.stim.letters.targ_color = [255 0 0]; % RGB value (0-255) of target stim
        params.stim.letters.dist_color = [0 0 255]; % RGB value (0-255) of distractor stim
%         params.stim.letters.letter_size = 70; % size of letters through TextSize
        params.stim.letters.letter_font = 'Arial Bold';

        params.stim_min_dist = 2; % min distance between two stimuli within the same hemifield
    case 'letters_orientation'
        % for orientation
        params.stim.orientation.feature_values = [0 90]'; % orientation of bars (0-360, 0 up- vertical, + clockwise)
        params.stim.orientation.targ_color = [255 0 0]; % RGB value (0-255) of target stim
        params.stim.orientation.dist_color = [0 0 255]; % RGB value (0-255) of distractor stim
        params.stim.orientation.bar_length_deg = 0.7; % deg
        params.stim.orientation.bar_width_deg = 0.23; % deg
        
        % may need to make stim specific values
        params.stim_min_dist = 1.1 * params.stim.orientation.bar_length_deg; % min distance between two stimuli within the same hemifield
        
        % for letters
        tmp_letters = 'RK';%'XO';
        params.stim.letters.feature_values = transpose(tmp_letters); % orientation of bars (0-360, 0 up- vertical, + clockwise)
        params.stim.letters.targ_color = [255 0 0]; % RGB value (0-255) of target stim
        params.stim.letters.dist_color = [0 0 255]; % RGB value (0-255) of distractor stim
%         params.stim.letters.letter_size = 70; % size of letters through TextSize
        params.stim.letters.letter_font = 'Arial Bold';

        params.stim_min_dist = 2; % min distance between two stimuli within the same hemifield
    otherwise
        error('unrecognized stimtype %s',params.stimtype)
end




%% datafile
% datfiles will be <subjid>_<experimentID>_<datecode>_<fileNumber>
c = clock;
params.time_stamp = sprintf('%02d/%02d/%04d %02d:%02d:%02.0f',c(2),c(3),c(1),c(4),c(5),c(6)); % month/day/year hour:min:sec
params.datecode = datestr(now,'mmddyy');

% get input
params.subjid = input('Enter Subject Code:','s');
params.runid  = input('Enter Run:');

% determine datafile
params.datafile = sprintf('%s_%s_%s_%03d',params.subjid,params.experiment,params.datecode,params.runid);
params.datadir = '../data/';

% check to see if this file exists
if exist(fullfile(params.datadir,[params.datafile '.mat']),'file')
    tmpfile = input('File exists.  Overwrite? y/n:','s');
    while ~ismember(tmpfile,{'n' 'y'})
        tmpfile = input('Invalid choice. File exists.  Overwrite? y/n:','s');
    end
    if strcmp(tmpfile,'n')
        display('Bye-bye...');
        return; % will need to start over for new input
    end
end


%% prepare psychtoolbox
% Opens psychtoolbox window
[w, rect] = Screen('OpenWindow',params.screen_num,params.bg_color,params.win_size);

% turned off for debugging
if ~debugging
    HideCursor; % Hide mouse cursor from screen
    ListenChar(2); % stops keypresses during experiment
    topPriorityLevel = MaxPriority(w); % Maximum priority level for the current operating system
    Priority(topPriorityLevel);
end

% set center of screen
params.x_center = rect(3)/2;
xc = params.x_center; % for easy reference
params.y_center = rect(4)/2;
yc = params.y_center; % for easy reference

% establish size and ppd, need to change for crt
params.pix_per_deg = (rect(3)/params.mon_width_deg);
ppd = params.pix_per_deg; % for easy reference

% draw to screen so we know things are working...
AddText(w, 'loading, please wait...', [xc yc], params.text_color, 0, params.instruction_text_size);
Screen(w,'Flip');

% convert sizes in deg to pixels
% and define rects here so we don't have to keep replicating these long bits of code
fix.innerwidth_pix     = params.fix.innerwidth_deg * ppd;
fix.innerhalfwidth_pix = fix.innerwidth_pix / 2;
fix.inner_rect = [xc-fix.innerhalfwidth_pix, yc-fix.innerhalfwidth_pix, xc+fix.innerhalfwidth_pix, yc+fix.innerhalfwidth_pix]; % rect for drawing inner portion of fixation spot
fix.outerwidth_pix = params.fix.outerwidth_deg * ppd;
fix.outerhalfwidth_pix = fix.outerwidth_pix / 2;
fix.outer_rect = [xc-fix.outerhalfwidth_pix, yc-fix.outerhalfwidth_pix, xc+fix.outerhalfwidth_pix, yc+fix.outerhalfwidth_pix]; % rect for drawing outer portion of fixation spot

switch params.stimtype
    case 'orientation'
        stim.orientation.length_pix = params.stim.orientation.bar_length_deg * ppd;
        stim.orientation.halflength_pix = stim.orientation.length_pix / 2;
        stim.orientation.width_pix = params.stim.orientation.bar_width_deg * ppd;
        stim.orientation.halfwidth_pix = stim.orientation.width_pix / 2;
        % and associated rects...
        stim.orientation.rect = [xc-stim.orientation.halfwidth_pix yc-stim.orientation.halflength_pix xc+stim.orientation.halfwidth_pix yc+stim.orientation.halflength_pix]; % template rect for a vertical (0 deg orientation) bar at center of screen
    case 'color'
        stim.color.width_pix = params.stim.color.width_deg * ppd;
        stim.color.halfwidth_pix = stim.color.width_pix / 2;
        % and associated rects...
        stim.color.rect = [xc-stim.color.halfwidth_pix yc-stim.color.halfwidth_pix xc+stim.color.halfwidth_pix yc+stim.color.halfwidth_pix]; % template rect for all stimuli (circle or square)
    case 'squares'
        stim.squares.width_pix = params.stim.squares.width_deg * ppd;
        stim.squares.halfwidth_pix = stim.squares.width_pix / 2;
        % and associated rects...
        stim.squares.points = [...
            xc-stim.squares.halfwidth_pix yc-stim.squares.halfwidth_pix; ...
            xc+stim.squares.halfwidth_pix yc-stim.squares.halfwidth_pix; ...
            xc+stim.squares.halfwidth_pix yc+stim.squares.halfwidth_pix; ...
            xc-stim.squares.halfwidth_pix yc+stim.squares.halfwidth_pix]; % template rect for all stimuli (squares or diamonds)
    case 'letters'
        % nothing to do
    case 'letters_orientation'
        % set for orientation
        stim.orientation.length_pix = params.stim.orientation.bar_length_deg * ppd;
        stim.orientation.halflength_pix = stim.orientation.length_pix / 2;
        stim.orientation.width_pix = params.stim.orientation.bar_width_deg * ppd;
        stim.orientation.halfwidth_pix = stim.orientation.width_pix / 2;
        % and associated rects...
        stim.orientation.rect = [xc-stim.orientation.halfwidth_pix yc-stim.orientation.halflength_pix xc+stim.orientation.halfwidth_pix yc+stim.orientation.halflength_pix]; % template rect for a vertical (0 deg orientation) bar at center of screen
        % nothing for letters
    otherwise
        error('unrecognized stimtype %s',params.stimtype)
end

cue.width_pix = params.cue.width_deg * ppd;
cue.halfwidth_pix = cue.width_pix / 2;
cue.eccen_pix = params.cue.eccen * ppd;
cue.rect(1,:) = [xc-cue.halfwidth_pix-cue.eccen_pix yc-cue.halfwidth_pix xc+cue.halfwidth_pix-cue.eccen_pix yc+cue.halfwidth_pix]; % rect for a LEFT block cue
cue.rect(2,:) = [xc-cue.halfwidth_pix+cue.eccen_pix yc-cue.halfwidth_pix xc+cue.halfwidth_pix+cue.eccen_pix yc+cue.halfwidth_pix]; % rect for a RIGHT block cue
cue.rect(3,:) = [xc-cue.halfwidth_pix yc-cue.halfwidth_pix xc+cue.halfwidth_pix yc+cue.halfwidth_pix]; % rect for a PASSIVE (CENTERED) block cue

% measure the frame rate
params.frame_rate = Screen('FrameRate',w); % in seconds.  this is what the operating system is set to, does not work on some systems (e.g., lcd or laptops)
params.flip_interval = Screen('GetFlipInterval',w); % in seconds.  this is an actual measurement, should work on all systems
params.flip_interval_correction = .80 * params.flip_interval; % this should work even on laptops that don't return a FrameRate value

% Enable alpha blending for contrast manipulations (like images with alpha channels)
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


%% blocks structure
% block structure will help us create the trials structure correctly
% each block is defined by:
%    0. task: passive vs. memory
%   (all following should ideally be crossed)
%    1a. stimtype
%    1. set size
%    2. memory side (i.e., location of targets in sample and probe)
%   and possibly, in the future...
%    (3. stimtype: visual vs. verbal)
%    (4. task: maintain vs. manipulate)
% The block structure will also be important for the block-cue period,
% which will indicate set_size and mem_side to the subject for the
% subsequent 6 trials.
blocks.n = []; % (placeholder) how many total blocks for this run?
blocks.order = []; % (placeholer) actual order of blocks
[tmp1, tmp2] = meshgrid(params.blocks.set_size,params.blocks.mem_side);
unique_exp_blocks = [tmp1(:) tmp2(:)]; % unpack meshgrid output
all_exp_blocks = repmat(unique_exp_blocks,params.blocks.memory_reps,1); % account for expeirmental block repetitions
blocks.stimtype = repmat({params.stimtype},size(all_exp_blocks,1),1);
blocks.set_size = all_exp_blocks(:,1);
blocks.mem_side = all_exp_blocks(:,2);
blocks.task_label = {'passive' 'memory'};
blocks.task = 0*blocks.set_size+2;
% add passive-viewing blocks
if params.blocks.npassive
    blocks.stimtype(end+1,1) = repmat({params.stimtype},params.blocks.npassive,1);
    blocks.task(end+1,1) = ones(params.blocks.npassive,1);
    blocks.set_size(end+1,1) = repmat(Sample(params.blocks.set_size),params.blocks.npassive,1); % *** randomly selected, not sure if this is ideal
    blocks.mem_side(end+1,1) = repmat(Sample(params.blocks.mem_side),params.blocks.npassive,1); % *** randomly selected, not sure if this is ideal
end

% randomize order of blocks
[~, idx_order] = Shuffle(blocks.set_size);
for fncell = fieldnames(blocks)';
    fn = fncell{1}; % convert to string
    if size(blocks.(fn),1) == length(idx_order)
        blocks.(fn) = blocks.(fn)(idx_order,:);
    end
end

% update n and order
blocks.n = length(blocks.set_size);
blocks.order = 1:blocks.n;


%% trials structure
% initialize trial structure, which will hold all of the trial-specific information in the order that it will be presented
trials.n = NaN; % total number of trials
trials.order = []; % an ordered list of trial id for the actual presentation order
trials.block = []; % block index of this set of trials.
trials.stimtype = []; % the type of stimulus, inherently defining the features by which change can occur
trials.set_size = []; % how many target stimuli on this trial? (must be <= params.stim_per_hemifield
trials.mem_side = []; % which side are the target stimuli presented? -1, 0, 1, 2 {left, passive, right, either}
trials.task     = []; % task index, see blocks.task_label
trials.targ_change = []; % boolean, does this trial contain a target-orientation change?
trials.targ_change_idx = []; % index (into targ_ lists) of the item that changed.  zero for no targ_change.  
trials.targ_side = []; % which side are the target presented ON THIS TRIAL? -1 1 [left right]
trials.targ_feature_sample = {}; % 'feature' of each target (to-be-remembered) element during SAMPLE presentation.  what this refers to depends on stimtype
trials.targ_feature_probe  = {}; % 'feature' of each target (to-be-remembered) element during PROBE presentation.  what this refers to depends on stimtype
trials.targ_pos_deg = {}; % xy position of each target element in degrees
trials.targ_pos_pix = {}; % xy position of each target element in pixels
trials.dist_feature     = {}; % orientation of each distractor (not-to-be-remembered) element (split by side, or just concatenated??)
trials.dist_pos_deg = {}; % xy position of each distractor element in degrees
trials.dist_pos_pix = {}; % xy position of each distractor element in pixels

% loop over each block
for ib = 1:blocks.n
    % label this block (in order of creation
    trials.block = cat(1,trials.block,repmat(ib,params.trials_per_block,1));
    
    % pull directly from blocks structure
    for fncell = {'stimtype' 'set_size' 'mem_side' 'task'}
        fn = fncell{1}; % convert to string
        trials.(fn) = cat(1,trials.(fn),repmat(blocks.(fn)(ib),params.trials_per_block,1));
    end
    
    switch blocks.mem_side(ib)
        case 2 % EITHER left OR right (but not a mix)
            trials.targ_side = cat(1,trials.targ_side,RandSample([-1 1],[params.trials_per_block 1])); % randomly sample from left (-1) or right (1)
        case {-1 1} % ONLY left OR right
            trials.targ_side = cat(1,trials.targ_side,repmat(blocks.mem_side(ib),params.trials_per_block,1));
        otherwise
            error('not sure how to handle a block mem_side of %d',blocks.mem_side(ib))
    end
end

% update n and order
trials.n = length(trials.block); % extract final number of trials
trials.order = (1:trials.n)';

% randomly fill in target/distractor orientation/position, regardless of
% block or order
trials.targ_change = rand(trials.n,1)<params.prop_targ_change;
trials.targ_change_idx = []; % we will fill this in later, because we want to be careful that we don't have a situation where all the targets have the same feature value in either the sample or the probe period
for it = 1:trials.n
    % remember to account for both hemifields
    tmp_feat_idx = 1; % just something to get us into the while loop...
    while length(unique(tmp_feat_idx))==1
        tmp_feat_idx = RandSample(1:size(params.stim.(trials.stimtype{it,1}).feature_values,1),[trials.set_size(it) 1]); % random row indices into params.stim.(trials.stimtype).feature_values
    end
    trials.targ_feature_sample(it,:) = {params.stim.(trials.stimtype{it,1}).feature_values(tmp_feat_idx,:)};
    trials.targ_feature_probe(it,:) = trials.targ_feature_sample(it,:); % start by assuming no targ change, so just copy sample to probe
    if trials.targ_change(it,1)
        % then we need to account for a change when defining probe targets
        % but don't allow a change that causes all targets to have the same
        % feature value
        tmp_probe_feat = 1;  % just something to get us into the while loop...
        first_loop = 1; % force our way into the while lo
        %%iii = 0; % debugging
        while length(unique(tmp_probe_feat))==1 || first_loop
            %%iii = iii+1; % debugging
            %%fprintf('%d.',iii) % debugging
            first_loop = 0; % don't force ourselves back into this loop
            tmp_targ_idx = Randi(trials.set_size(it,1)); % random index of target that will change
            tmp_samp_feat = trials.targ_feature_probe{it,:}(tmp_targ_idx,:); % sample value for this target
            [~, tmp_samp_feat_idx] = ismember(tmp_samp_feat,params.stim.(trials.stimtype{it,1}).feature_values,'rows'); % which row index matches the current samp_feat
            tmp_valid_probe_feat_idx = [1:tmp_samp_feat_idx-1 tmp_samp_feat_idx+1:size(params.stim.(trials.stimtype{it,1}).feature_values,1)]; % all possible row indices to choose probe feature from (not including the one that matches the sample feature)
            tmp_probe_feat_idx = Sample(tmp_valid_probe_feat_idx); % the row index of the new probe feature
            tmp_probe_feat = trials.targ_feature_sample{it,:}; % start by assuming no targ change, so copy sample to probe
            tmp_probe_feat(tmp_targ_idx,:) = params.stim.(trials.stimtype{it,1}).feature_values(tmp_probe_feat_idx,:);
        end
        %%fprintf('\n') % debugging
        trials.targ_feature_probe{it,:} = tmp_probe_feat;
    end
    % verify that none of the sample/probe contain all targets of the same feature value
    if sum(cellfun(@(x) length(unique(x)),trials.targ_feature_sample)==1) || sum(cellfun(@(x) length(unique(x)),trials.targ_feature_probe)==1)
        sca; % close screen
        error('oops, looks like there is a probe or sample with all targets having the same feature value.  that shouldn''t happen');
    end
    % fill in distractor features completely randomly
    tmp_feat_idx = RandSample(1:size(params.stim.(trials.stimtype{it,1}).feature_values,1),[2*params.stim_per_hemifield-trials.set_size(it) 1]); % random row indices into params.stim.(trials.stimtype).feature_values
    trials.dist_feature(it,:) = {params.stim.(trials.stimtype{it,1}).feature_values(tmp_feat_idx,:)};
 
    % get targ and dist positions together
    left_td_pos_deg = PickPoints(params.stim_per_hemifield,[],params.stim_min_dist,[params.hemifield_width_deg params.hemifield_height_deg]);
    left_td_pos_deg(:,1) = left_td_pos_deg(:,1) - params.hemifield_eccen_deg; % offset xs
    right_td_pos_deg = PickPoints(params.stim_per_hemifield,[],params.stim_min_dist,[params.hemifield_width_deg params.hemifield_height_deg]);
    right_td_pos_deg(:,1) = right_td_pos_deg(:,1) + params.hemifield_eccen_deg; % offset xs
        
    % pack these lists, together, into cell arrays
    switch trials.targ_side(it)
        case -1 % targets on the left
            trials.targ_pos_deg{it,1} = left_td_pos_deg(1:trials.set_size(it),:);
            trials.dist_pos_deg{it,1} = cat(1,left_td_pos_deg(trials.set_size(it)+1:end,:),right_td_pos_deg);
            
        case 1 % targets on the right
            trials.targ_pos_deg{it,1} = right_td_pos_deg(1:trials.set_size(it),:);
            trials.dist_pos_deg{it,1} = cat(1,right_td_pos_deg(trials.set_size(it)+1:end,:),left_td_pos_deg);
            
        otherwise
            error('not sure how to handle a trial targ_side of %d',traisl.targ_side(it))
    end
 
    % convert from deg -> pix
    trials.targ_pos_pix{it,1} = Degrees2PixelsPTB(trials.targ_pos_deg{it,1},[xc yc], ppd);
    trials.dist_pos_pix{it,1} = Degrees2PixelsPTB(trials.dist_pos_deg{it,1},[xc yc], ppd);
end

%% data structure
% initialize and preallocate memory for data structure, which will hold
% trials-specific data like subjects resp, rts and status (correct/incorrect)
nanfill = NaN(trials.n,1);
data.resp     = nanfill; % keycode of response button pressed
data.respidx  = nanfill; % index into params.resp_keys of (first) button pressed
data.rt       = nanfill; % response time relative to probe onset
data.resptime = nanfill; % response time relative to trials start
data.status   = nanfill; % correct (1) or incorrect (0)
data.hit      = nanfill; % actual target change, resp target change
data.miss     = nanfill; % actual target change, resp no target change
data.false_alarm = nanfill;% no actual target change, resp target change
data.correct_reject = nanfill;% no actual target change, resp no target change

%% data structure
% initialize and preallocate memory for timing structure, which will hold
% timing information for all phases of the run
for scell = {'run' 'prefix' 'postfix'}
    s = scell{1}; % convert to string
    timing.([s '_start'])  = NaN;
    timing.([s '_stop'])   = NaN;
    timing.([s '_dur'])    = NaN;
end
block_nanfill = NaN(blocks.n,1);
for scell = {'block' 'block_cue' 'block_cue_blank'}
    s = scell{1}; % convert to string
    timing.([s '_start'])  = block_nanfill;
    timing.([s '_stop'])   = block_nanfill;
    timing.([s '_dur'])    = block_nanfill;
end
for scell = {'trial' 'sample' 'delay' 'probe' 'iti' 'respcheck'}
    s = scell{1}; % convert to string
    timing.([s '_start'])  = nanfill;
    timing.([s '_stop'])   = nanfill;
    timing.([s '_dur'])    = nanfill;
end
%% setup keypress queue
% better timing than KbCheck and doesn't require while loops (i.e., works in background)
allowed_keys = zeros(1,256); % by default, keys are not allowed to influence queue
allowed_keys([params.resp_keys params.abort_key]) = 1; % only listen for valid response keys and the abort key
KbQueueCreate(params.kb_dev_id,allowed_keys);
KbQueueStart(params.kb_dev_id); % start listening for keypresses


%% instructions
% Display instructions
AddText(w, 'Please maintain fixation throughout run', [xc yc], params.text_color, 6, params.instruction_text_size);
AddText(w, 'Peripheral cue = memory task on one side or either', [xc yc], params.text_color, 4.1, params.instruction_text_size);
AddText(w, 'Does any cued target change? Left Shift=NO, Right Shift=YES', [xc yc], params.text_color, 3, params.instruction_text_size);
AddText(w, 'Central cue = passive task (no memory)', [xc yc], params.text_color, 1.1, params.instruction_text_size);
AddText(w, 'Press any key when second stimulus appears', [xc yc], params.text_color, 0, params.instruction_text_size);
AddText(w, 'Press a key when ready to begin', [xc yc], params.text_color, -3, round(.80 * params.instruction_text_size)); % slightly smaller...
Screen(w,'Flip');

% wait for subject to press a key
[~, keycode] = KbPressWait(params.kb_dev_id);
while ~keycode
    [~, keycode] = KbPressWait(params.kb_dev_id);
end

%% waiting for scanner trigger (or self-start with "return" for behavioral)
AddText(w, 'Ready for the next run.', [xc yc], params.text_color, 3, params.instruction_text_size);
AddText(w, 'Waiting for trigger...', [xc yc], params.text_color, 0, params.instruction_text_size);
AddText(w, '(press return key if self-triggering)', [xc yc], params.text_color, -3, round(.80 * params.instruction_text_size)); % slightly smaller...
Screen(w,'Flip');

% wait for subject to press the Return key to begin experiment after .5 secs
[~, keycode] = KbPressWait(params.kb_dev_id);
while ~any(keycode([params.return_key params.trigger_key]))
    [~, keycode] = KbPressWait(params.kb_dev_id);
end


%% prefix period
% add fixpot (keep for duration of expeirment)
Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect);
timing.run_start = Screen('Flip',w);
timing.prefix_start = timing.run_start;
timing.prefix_stop  = WaitSecs(params.prefix_dur);


%% trial loop
ib = 0; % block counter
this.block = NaN; % start with a null value so that first trial is recognized as a new block
abort = 0;
for it = 1:trials.n;
    if abort; break; end; % kick out of trial loop if the abort key was pressed on the last trial
    
    
    % check to see if this is the start of a new block, and if so, show block cue
    if trials.block(it,1) ~= this.block
        % ...then this trial is of a different block than the last trial
                
        % mark last block as stopped (as long as it is not the first block)
        tmp_block_time = GetSecs; % use GetSecs since we can get here coming from multiple previous places
        if ib
            timing.block_stop(ib,1) = tmp_block_time;
        end

        % update block counters and record block start
        ib = ib+1;
        timing.block_start(ib,1) = tmp_block_time;
        this.block = trials.block(it,1); % update
        
        % show cue for this block
        % *** CUES MAY NEED TO BE STIMTYPE SPECIFIC, IF STIMTYPE INTERLEAVED
        if blocks.task(ib,1) == 1
            Screen('FrameOval',w,params.cue.color,cue.rect(3,:),1,1);
        else
            switch blocks.mem_side(ib,1)
                case -1 % targets on the left
                    Screen('FillOval',w,params.cue.color,cue.rect(1,:));
                case 1 % targets on the right
                    Screen('FillOval',w,params.cue.color,cue.rect(2,:));
                case 2 % targets on either side
                    Screen('FillOval',w,params.cue.color,cue.rect(1,:));
                    Screen('FillOval',w,params.cue.color,cue.rect(2,:));
                otherwise
                    error('not sure how to handle a trial targ_side of %d',traisl.targ_side(it))
            end
        end
        Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect);
        timing.block_start(ib,1) = Screen('Flip',w);
        timing.block_cue_start(ib,1) = timing.block_start(ib,1);
        WaitSecs(params.block_cue_dur - params.flip_interval_correction); % correct for the fact that we won't actually flip the screen for another frame
        
        % clear to fixation spot
        Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect);
        timing.block_cue_blank_start(ib,1) = Screen('Flip',w);
        timing.block_cue_stop(ib,1) = timing.block_cue_blank_start(ib,1);
        timing.block_cue_blank_stop(ib,1) = WaitSecs(params.block_cue_blank_dur);
    end
    
    
    % mark trial start time
    timing.trial_start(it,1) = GetSecs; % use GetSecs since we can get here coming from multiple previous places
    
    % flush the KbQueue - note that this will mean responses from previous trials will not spill over
    KbQueueFlush(params.kb_dev_id);
    
    
    % extract any trial-specific params into "this" structure
    for fncell = {'targ_feature_sample' 'targ_feature_probe' 'targ_pos_pix' 'dist_feature' 'dist_pos_pix'}
        fn = fncell{1}; % convert to string
        this.(fn) = trials.(fn){it,:}; % unpack from cell array
    end
    if trials.task(it,1) == 1
        % passive-viewing, any response is valid
        this.correct_resp = 1:length(params.resp_keys);
    else
        % working memory task, correct response depends on if there was a
        % change from sample to probe
        this.correct_resp = trials.targ_change(it,1) + 1; % targ_change 0/1 no change/change -> response 1/2 left/right
    end
    
    
    % sample period
    this.stim_phase = 'sample'; % for draw_stim
    draw_stim
    % draw fixspot, flip, wait
    Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect); % draw double fix spot
    timing.sample_start(it,1) = Screen(w,'Flip');
    WaitSecs(params.sample_dur); % this is typically really short, so don't use correction since it could be fairly substantial
    
    
    % delay period
    Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect); % draw double fix spot
    timing.delay_start(it,1) = Screen(w,'Flip');
    timing.sample_stop(it,1) = timing.delay_start(it,1);
    WaitSecs(params.delay_dur - params.flip_interval_correction); % correct for the fact that we won't actually flip the screen for another frame
    
    
    % probe period
    KbQueueFlush(params.kb_dev_id); % only consider responses that occur after the probe is displayed
    this.stim_phase = 'probe'; % for draw_stim
    draw_stim
    % draw fixspot, flip, wait
    Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect); % draw double fix spot
    timing.probe_start(it,1) = Screen(w,'Flip');
    timing.delay_stop(it,1) = timing.probe_start(it,1);
    WaitSecs(params.probe_dur - params.flip_interval_correction); % correct for the fact that we won't actually flip the screen for another frame
    
    
    % iti period
    Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect); % draw double fix spot
    timing.iti_start(it,1) = Screen(w,'Flip');
    timing.probe_stop(it,1) = timing.iti_start(it,1);
    % to avoid letting small timing errors build up, time the iti relative
    % to the run time based on the expected elapsed time up to this point.
    this.exp_trial_stop = timing.run_start + params.prefix_dur+ ib*(params.block_cue_dur+params.block_cue_blank_dur) + it*(params.sample_dur+params.delay_dur+params.probe_dur+params.iti_dur);
    timing.respcheck_start(it,1) = WaitSecs('UntilTime',this.exp_trial_stop);
    
    % check to see if we got a response and mark if it was 'correct' or 'incorrect'
    [press firstPress] = KbQueueCheck(params.kb_dev_id); % check response queue
    if press
        % should we abort?
        if firstPress(params.abort_key)
            abort = 1;
            continue
        end
        
        % first detected response key - record and validate depending on current task and trial params
        resp_keys_time = firstPress(params.resp_keys); % time (GetSecs clock) of the first time each possible resp_key was pressed since last queue flush
        
        % it is possible that more than one response key was pressed during the trial (if more than one was allowed)
        % in this case, we'll only record the first key that was pressed.
        resp_keys_time(resp_keys_time==0) = NaN; % convert 0s to NaNs, so the 0s are not identified by min in the next line of code:
        first_resp_idx = find(resp_keys_time == min(resp_keys_time));
      
        % N.B. In the insanely-virtually-impossible case that two keys
        %      were pressed at EXACTLY the same time, we'll arbitrarily
        %      select one as the recorded response.  Not sure if this is the
        %      best thing to do, but (1) it should NEVER happen and (2) it
        %      would screw up the data structure if we try to record multiple
        %      responses.
        if length(first_resp_idx) > 1
            % then multiple keys have been pressed.  arbitrarily select one as the recorded response.
            first_resp_idx = Sample(first_resp_idx);
        end
        
        % record response as the keycode of the button pressed
        data.resp(it,1) = params.resp_keys(first_resp_idx);
        data.respidx(it,1) = first_resp_idx;
        
        % record status (correct/incorrect) and response time (relative to probe onset and trial start)
        % N.B. rt is relative to stim_on; resptime is relative to trial_start, but
        %      we can recalculate offline to other times
        data.rt(it,1)       = resp_keys_time(first_resp_idx) - timing.probe_start(it,1);
        data.resptime(it,1) = resp_keys_time(first_resp_idx) - timing.trial_start(it,1);
        
        % mark status
        data.status(it,1) = ismember(first_resp_idx,this.correct_resp); % does actual response match correct response (this.correct_rep and first_resp_idx have matching indices that fit with params.resp_keys)

        % categorize as hit, miss, false_alaram, correct_reject
        % N.B. the '1*' is there to make sure these are stored as integers,
        % not logicals, else if there is a no-response trial the NaN cannot
        % be added to a logical list
        data.hit(it,1)      = 1*and(trials.targ_change(it,1),data.respidx(it,1)==2); % actual target change, resp target change
        data.miss(it,1)     = 1*and(trials.targ_change(it,1),data.respidx(it,1)==1); % actual target change, resp no target change
        data.false_alarm(it,1) = 1*and(~trials.targ_change(it,1),data.respidx(it,1)==2);% no actual target change, resp target change
        data.correct_reject(it,1) = 1*and(~trials.targ_change(it,1),data.respidx(it,1)==1);% no actual target change, resp no target change

    else
        % no button pressed
        data.resp(it,1)     = 0; % this in an index into params.resp_keys, so zero means no response
        data.rt(it,1)       = NaN; % no rt without a response
        data.resptime(it,1) = NaN; % no rt without a response
        data.status(it,1)   = 0; % non-response trials are considered incorrect
        data.hit(it,1)      = NaN; % without a response, this is ill-defined
        data.miss(it,1)     = NaN; % without a response, this is ill-defined
        data.false_alarm(it,1) = NaN; % without a response, this is ill-defined
        data.correct_reject(it,1) = NaN; % without a response, this is ill-defined

    end
    
    % save data file after every trial in case of crash
    % N.B. this will likely lead to small timing errors, especially on the
    % very first trial, because the save will take some time (~5 ms on my
    % laptop).  but this will be accounted for by adjustments to following
    % ITIs since they are measured relative to the expected time so that
    % the fMRI experiment doesn't get far off of the TR.
    save(fullfile(params.datadir,params.datafile));
    
    timing.iti_stop(it,1)   = GetSecs;
    timing.respcheck_stop(it,1) = timing.iti_stop(it,1);
    timing.trial_stop(it,1) = timing.iti_stop(it,1);
end


%% postfix period
% add fixpot (keep for duration of expeirment)
Screen('FillRect',w,params.fix.outer_color,fix.outer_rect); Screen('FillRect',w,params.fix.inner_color,fix.inner_rect);
timing.postfix_start = Screen('Flip',w);
% mark stop time of final block
timing.block_stop(ib,1) = timing.postfix_start;
% wait postfix time
timing.postfix_stop  = WaitSecs(params.postfix_dur);


%% notify subject that run is complete
perf.overall = 100 * nanmean(data.status);
perf.correct_detect = 100 * nanmean(data.status(trials.targ_change));
AddText(w, 'Run Complete!', [xc yc], params.text_color, 1, params.instruction_text_size);
AddText(w, sprintf('Performance = %.02f%%',perf.overall), [xc yc], params.text_color, -1, params.instruction_text_size);
timing.run_stop = Screen('Flip',w);

%% post-run calculations and cleanup
% calculate timing durations
for scell = {'run' 'prefix' 'postfix' 'block' 'block_cue' 'block_cue_blank' 'trial' 'sample' 'delay' 'probe' 'iti' 'respcheck'}
    s = scell{1}; % convert to string
    timing.([s '_dur']) = timing.([s '_stop']) - timing.([s '_start']);
end

% save data file
save(fullfile(params.datadir,params.datafile));

% Close window and restoﬂre normal functioning
WaitSecs(3); % make sure subject sees the final message
KbQueueRelease; % release keypress queue or process may continue to run in background
ShowCursor;
ListenChar(0);
Priority(0);
Screen('CloseAll')


%% feedback on console
% ...now that Screen is closed (since PTB writes lots of things to console): some info about last file/run for easy tracking
display(sprintf('\n\n\n*****************************************'));
display(sprintf('   last file: %s', params.datafile));
display(sprintf('   last run : %d', params.runid));
display(sprintf('   performance = %.2f%%',perf.overall));
display(sprintf('*****************************************\n'));