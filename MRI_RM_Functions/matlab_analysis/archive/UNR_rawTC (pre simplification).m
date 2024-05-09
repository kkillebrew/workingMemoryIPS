function [TC subj] = rawTC(subj_initials,runs,rois,hemisphere,vox_select,varargin)
% rawMVPA
%   [TC SUBJ] = rawTC(subj_initials,rois,hemisphere,vox_select,varargin)
%
% Currently assumed to be run from .../CLab/[Expeirment]/matlab_analysis/
% This function assumes that the PU MVPA toolbox is installed.
%
% Generalized script for ploting BOLD Time Course.  Data is first extracted
% using the Princeton Univeristy MVPA toolbox, so it is very comperable to
% the extraction of data for rawMVPA.  Note, however, that normalized
% functional data are extracted by default.
%
% Arguments
%   subj_initials = subject initials, which tells us directory and filenames for this subject.
%   runs        = list of runs to analyze.  this is different than older 'task' definition
%                 which required a subject database, which I dropped at UNR for simplicity.
%   rois        = a cell array of roi names to test, or 'visual', 'parietal', or 'all'
%                to test all visual, parietal or both.
%   hemisphere  = 'left' 'right or 'both' (defualt).  Which hemisphere to consider.
%   vox_select = defines how voxels are selected within each ROI (voxel selection).
%            {'fullFstat' thresh} = This is the (abs) of the Full_Fstat statistic...(like tstat option)
%            {'tstat' thresh} = This is the (abs) t-statistic that must be
%                   exceeded in the EXP vs FIX comparison to be included as a valid
%                   ROI voxel. (t-value of 1.98 corresponds roughly to p<0.05,
%                   depending on df).
%            {'nvox' thresh}  = This sets the max number of voxels to take from each
%                   ROI.  Voxels are selected by their t-value (positive only).  So the 'thresh'
%                   most-active voxels are selected.  Rois with less than thresh voxels
%                   will use all available voxels. See RSA.nvox output.
%            {'nvoxabs' thresh}  = This sets the max number of voxels to take from each
%                   ROI.  Voxels are selected by their absolute-value t-value (positive and negative.  
%                   So the 'thresh' most-active voxels are selected.  Rois with less than thresh voxels
%                   will use all available voxels. See RSA.nvox output.
%            {'anova' alpha} = Use a non-peaking anova just prior to cross_validation and
%                   accept voxels that are significantly modulated across condition
%                   at an alpha level of 'alpha'.
%            {'minbeta' thresh} = The beta weight must exceed 'thresh' in the EXP vs FIX comparison to
%                   be included as a valid ROI voxel.  Only positive beta weights are considred (EXP>FIX)
%            'none' = (default) no voxel selection is used. Note that this is not 
%                   necessarilly equivilent to setting a tstat threshold to 0.%   optional args [default]:
%
%   optional args [default]:
%       preTRs  [2] = how many TRs to show prior to onset of each block/event.  Must be >=0. Use zero to
%                     start timecourse aligned with event/block starts as defined by stimtimes.
%       postTRs [16] = how many TRs to show after the onset of each block/event.  Must be >0. Should span the event/block
%                     length and any time points afterwards for the HRF to come back to baseline.
%       baseline_idx [-preTRs:0] = index, relative to block/event onset (idx 0) overwhich average
%                                     signal is extracted, defined as baseline, and subtracted from timecouse.
%                                     If empty ([]), baseline is set to zero and timecourse is not shifted.
%
%       norm_type ['each'] = 'each' to normalize each event/block individually by mean activity in baseline_idx
%                            'pre' to load pre-normalized data (by adding '_norm' to prefix_base)
%                            'grand' to normalize to mean activity in baseline_idx across entire run
%
%       smoothing [0] = what smoothing level to use.
%       prefix_base ['tscvr'] = define prefix for bucket/datafiles, which defines preprocessing steps
%                 (not including smoothing or normalization).  Default value assumes undistortion was NOT performed.
%       exp_bucket_prefix [''] = prefix (no subject, but including smoothing) for bucket used by voxel selection
%                         (nvox or tstat).  if empty, then prefix supplied by prefix_base is used (with smoothing).
%                         if empty, then prefix supplied by prefix_base is used (with smoothing). Additionally,
%                         '_norm' will be added unless 'norm_type' is set to 'pre' (in which case it should be included in 'prefix_base').
%                         exp_bucket_suffix suffix is always added on the end.
%       exp_bucket_suffix ['_bucket'] = string appended to end of exp_bucket_prefix.
% %       exp_bucket_task [''] = task used in prefix for exp_bucket_prefix, if different than task argument.  If empty, 
% %                              task argument is used.
%
%       exp_reg ['EXP'] = what contrast is used for the "EXP"-based voxel selection?  By default, it is "EXP", which typically
%                         represents Experimental Condition vs. Baseline.
%
% %       reg_task [''] = the task to supply to ExtractStimTimes.  By default, and almost always, this should match the input argument 'task'.
% %                       however, there are very specific circumstances where you may wish to create an alternate set of regressor stimtimes
% %                       for the same set of runs.  For example, I added this code to account for an analysis where, by default, the "adapt"
% %                       task (which was used as input and defined the runs loaded) had regressors (and stimtimes) defined for 3 stimulus classes,
% %                       but I wanted to create timecourses collapsed across those stimulus classes.  One alternative is to create the timecourses
% %                       for each regressor separately and then average them together, but this assumes an equal number of stimtimes for each condition.
% %                       The more conservative approach was to create a new set of stimtimes that was a combination of those from the three stimulus classes
% %                       and rerun the rawTC analysis.  It would be inefficient to create another set of all the associated data needed for this analysis,
% %                       so this optional argument allows you to load in that special set of stimtimes for the same data.
% %                       N.B. this defines the output filename to ensure unique files when supplying the same task and two different reg_task's
%       reg_include ['all']  = 'all' or cell array list of regressors to include in analysis.
%       reg_exclude ['none'] = 'none' or cell array list of regressors to include in analysis.
%       block_filter ['none'] = what type of filter should we apply to remove outlier blocks.  'none' for no filtering, 'sse' to remove blocks that whos
%                               sum-squared error (over entire timecourse) exceeds 2 standard deviations of the mean SSE.  filter is applied separately
%                               for each ROI-condition.
%       output_tail [''] = string to add to the end of the output file.  An '_' will be used to separate this
%                          string from the rest of the filename.
%
%       USAGE EXAMPLE: rawTC(... ,'prefix_base','mytimecourse','reg_include',{'cond1' 'cond2'})
%
% Return
%   TC      = timecourse structure
%   SUBJ    = final subject structure

% REHMB 01.09 - created, based on rawMVPA.  Plan to setup like classification,
%               then pull out data to plot time courses.
%       09.09 - default vox_select -> 'none' and now an empty vox_select is converted to 'none'.  makes for more explicit options storage and easier to remember usage.
%       09.09 - separate bucket loading for data and EXP contrast with option for loading different buckets.
%       03.11 - major update based on previous modifications to glmFIR
%       10.11 - added reg_task option and changed outputfile name to use reg_task instead of task (note: reg_task defaults to task unless specifically changed)
%       04.14 - trying to adapt for UNR setup, which lacks subject-database and the idea of "tasks"

% TODO


%% set parameters
% subj_initials = 'KK'; % subject initials define where data is located and filenames
% runs = 1:7; % which functional runs to use?
% rois = {'V1' 'V2'}; % which ROIs to use?
% hemisphere = 'lh'; % 'lh', 'rh', or 'both'
% vox_select = {'tstat' 1.98}; % threhsold for voxel selection



%% validate parameters
% expand rois argument
rois = UNR_ROIArgs(rois);

% % hemisphere
% if nargin < 4
%     hemisphere = 'both';
% end
% if ~ismember(hemisphere,{'both' 'left' 'right'})
%     error('invalid hemisphere (%s).  Must be ''left'', ''right'', or ''both''', hemisphere)
% end
% 
% % voxel selection
% if nargin < 5 || isempty(vox_select)
%     vox_select = {'none'};
% end
% if ~iscell(vox_select)
%     vox_select = {vox_select};
% end


% parse option strings
defaults.preTRs              = 2;
defaults.postTRs             = 16;
defaults.baseline_idx        = NaN; % place-holder.  if undefined, will use -preTRs:0
defaults.norm_type           = 'each';
defaults.smoothing           = 0;
defaults.prefix_base         = 'tscvr';
defaults.exp_bucket_prefix   = '';
defaults.exp_bucket_suffix   = '_bucket';
% % defaults.exp_bucket_task     = '';
defaults.exp_reg             = 'EXP';
% % defaults.reg_conds           = 'all';
% % defaults.reg_task            = task; % by default, this matches the required input argument
defaults.reg_include         = 'all';
defaults.reg_exclude         = 'none';
defaults.block_filter        = 'none';
defaults.output_tail         = '';
defaults.verbose             = 1; % NOT FULLY IMPLEMENTED YET
options = propval(varargin,defaults);

% propgate some defaults based on potential input args
if isnan(options.baseline_idx)
    % use default value.  propagate from preTRs
    options.baseline_idx = -options.preTRs:0;
end

% validate some optional arguments
if options.preTRs < 0
    error('preTRs cannot be negative')
end
if options.postTRs < 1
    error('postTRs must be at least 1, else there would be nothing of interest to plot')
end
if ~strcmp(options.norm_type,'pre') && isempty(options.baseline_idx)
    error('must define a valid baseline_idx when not using pre-normalized data')
end

% bookkeeping
TC.block_filter.type = options.block_filter;
TC.block_filter.notes = 'block filtering performed on TC.each.raw_norm data ONLY, by filling filtered rows with NaNs';


%% save args to TC
inargs = {'subj_initials' 'runs' 'rois' 'hemisphere' 'vox_select'};
for i = 1:length(inargs)
    eval(sprintf('TC.%s = %s;',inargs{i},inargs{i}));
end
TC.options = options;


%% initialize and parameters
% INITIALIZING THE SUBJ STRUCTURE
% % % sdb = SubjectDatabase(subj_initials); % this should be defined locally for each experiment
% % TC.subj = sdb.subj;
% % TC.date = sdb.date;
% % TC.patid = sdb.patid;
% % TC.sdb = sdb; % store full subject database

% start by creating an empty subj structure using PU-MVPA toolbox
subj = init_subj('rawtc',subj_initials);

% where's the data?
filedir = fullfile('../',subj_initials,'analysis');


% update prefix for smoothing
prefix = options.prefix_base; % start with base
if options.smoothing
    prefix = [prefix 'sm' num2str(options.smoothing)];
end
% % prefix = [prefix 'sm' num2str(options.smoothing)]; % alwyas expand, even for sm0

% update prefix for pre-normalized data
if strcmp(options.norm_type,'pre')
    prefix = [prefix '_norm'];
end



%% get ROIs and index definitions and create grand_mask using all requested ROIs
gotv = 0;
% % gotp = 0;
% % gotf = 0;
goto = 0;
% % gotc = 0;
% % gotn = 0;
try
    custom_rois = GetCustomROIs(-1,'names');
catch
    custom_rois = {'NONE_DEFINED'};
end
ROI.grand_mask = []; % init to null
for ridx = 1:length(rois)
    roi = rois{ridx};
    switch roi
        case UNR_GetRetROIs(-1,'names')
            if ~gotv % only load ROIs if we need them
                ROI.visual = UNR_GetRetROIs(subj_initials,'matrix');
                gotv = 1;
            end
            rois_brik{ridx} = 'visual';
% %         case UNR_GetParietalROIs(-1,'names')
% %             if ~gotp % only load ROIs if we need them
% %                 ROI.parietal = UNR_GetParietalROIs(subj_initials,'matrix');
% %                 gotp = 1;
% %             end
% %             rois_brik{ridx} = 'parietal';
% %         case UNR_GetFrontalROIs(-1,'names')
% %             if ~gotf % only load ROIs if we need them
% %                 ROI.frontal = UNR_GetFrontalROIs(subj_initials,'matrix');
% %                 gotf = 1;
% %             end
% %             rois_brik{ridx} = 'frontal';
        case UNR_GetObjROIs(-1,'names');
            if ~goto % only load ROIs if we need them
                ROI.object = UNR_GetObjROIs(subj_initials,'matrix');
                goto = 1;
            end
            rois_brik{ridx} = 'object';
% %         case custom_rois
% %             if ~gotc % only load ROIs if we need them
% %                 ROI.custom = UNR_GetCustomROIs(subj_initials,'matrix');
% %                 gotc = 1;
% %             end
% %             rois_brik{ridx} = 'custom';
% %         case GetNoiseROIs(-1,'names');
% %             if ~gotn % only load ROIs if we need them
% %                 ROI.noise = UNR_GetNoiseROIs(subj_initials,'matrix');
% %                 gotn = 1;
% %             end
% %             rois_brik{ridx} = 'noise';
        otherwise
            error('unknown ROI source for %s.', roi)
    end
    
    % add current ROI to grand_mask
    if isempty(ROI.grand_mask)
        % first one, need to initilize size
        ROI.grand_mask = zeros(size(ROI.(rois_brik{ridx}).data.lh)); % can use either hemisphere
    end
    % merge with existing grand_mask
    switch hemisphere
        case 'both'
            % collapse across left and right hemispheres
            ROI.grand_mask = or(ROI.grand_mask, ...
                or(...
                ismember(ROI.(rois_brik{ridx}).data.lh, ROI.(rois_brik{ridx}).idx.(roi)), ...
                ismember(ROI.(rois_brik{ridx}).data.rh, ROI.(rois_brik{ridx}).idx.(roi))));
        case 'left'
            % just left hemisphere
            ROI.grand_mask = or(ROI.grand_mask, ...
                ismember(ROI.(rois_brik{ridx}).data.lh, ROI.(rois_brik{ridx}).idx.(roi)));
        case 'right'
            % just right hemisphere
            ROI.grand_mask = or(ROI.grand_mask, ...
                ismember(ROI.(rois_brik{ridx}).data.rh, ROI.(rois_brik{ridx}).idx.(roi)));
    end
end


%% create grand_mask
subj = init_object(subj,'mask','grand_mask');
subj = set_mat(subj,'mask','grand_mask',ROI.grand_mask);


%% read/load EXP Bucket HEADER/DATA
% N.B. by default, the EXP bucket is used and referred to as EXP in this code, 
% but using the exp_reg optional argument, the contrast used for EXP is adjustable.
%
% N.B. This is only used for voxel selection, and may be same as above.
% start with header only so we can figure out which colums we'll need from the Bucket.
% we should at least be able to elimiate half the data ('b' vs 't'/'F')
if ismember(vox_select{1},{'fullFstat' 'tstat' 'tstatabs' 'nvox' 'nvoxabs' 'minbeta'})
    % get header info so we can limit loading to only EXP t-contrast
    if options.verbose; fprintf('Reading Bucket header for %s contrast voxel selection...',options.exp_reg); end
% %     if isempty(options.exp_bucket_task)
% %         options.exp_bucket_task = task;
% %     end

    if isempty(options.exp_bucket_prefix)
        options.exp_bucket_prefix = prefix;
        if ~strcmp(options.norm_type,'pre')
            % if time course is based on normalized data, this has already been done.  if not, we need to for the bucket.
            options.exp_bucket_prefix = [options.exp_bucket_prefix '_norm'];
        end
    end
    
    exp_bucket_file = ['_' options.exp_bucket_prefix options.exp_bucket_suffix];
    EXP_HEADER = UNR_LoadBucket(subj_initials,exp_bucket_file,[],1);
    if options.verbose; fprintf('...done\n'); end
    
    % load EXP t-contrast data
    if options.verbose; fprintf('Loading Bucket data...'); end
    bucket_options.Format = 'matrix';
    if ismember(vox_select{1},{'minbeta'})
        exp_bucket_param = 'b'; % beta weight
    else
        exp_bucket_param = 't'; % t-statistic
    end
    if strcmp(options.exp_reg,'Full_Fstat')
        exp_bucket_param = 'other'; % t-statistic
    end
    
    if ~isfield(EXP_HEADER.idx.(exp_bucket_param),options.exp_reg)
        error('couldn''t locate %s contrast in the provided EXP bucket',options.exp_reg)
    end
    bucket_options.Frames = EXP_HEADER.idx.(exp_bucket_param).(options.exp_reg); % only load GLMs that will be used + EXP contrast.  Necessary to save memory when loading very large buckets.
   
    EXP_BUCKET = UNR_LoadBucket(subj_initials,exp_bucket_file,bucket_options,0);
    if options.verbose; fprintf('done\n'); end
end
% %% load glm-bucket for thresholding?
% if ~isempty(vox_select) & ismember(vox_select{1},{'tstat' 'nvox'})
%     fprintf('Loading Bucket (EXP contrast only) for voxel selection using thresholding...');
%     HEADER = UNR_LoadBucket(subj_initials,['_' task '_' prefix '_bucket'],[],1);
%     bucket_options.Format = 'matrix';
%     bucket_options.Frames = HEADER.idx.t.([task '_EXP']); % for thresholding purposes, we only need EXP contrast's t value
%     BUCKET = UNR_LoadBucket(subj_initials,['_' task '_' prefix '_bucket'],bucket_options);
%     fprintf('\ndone\n')
% end


%% read in the actual data using the 'grand mask' as a starting point.
% only load runs that are part of this task
for i = runs
    fns{i} = fullfile(filedir,sprintf('%s_r%02d_%s+orig',subj_initials,i,prefix));
end

% now, read and set up the actual data. load_AFNI_pattern reads in the
% EPI data from a BRIK file, keeping only the voxels active in the
% grand mask (see above)
subj = load_afni_pattern(subj,'grand_pattern','grand_mask',fns);
curpat = 'grand_pattern';


%% regressors
% % % TODO: simplify ExtractStimTimes to only provide st?  i think we convert st to reg here.
% % [st reg con] = UNR_ExtractStimTimes(subj_initials,0,options.reg_task); % will provide reg structure (NB, reg_task might be different than task)

% the following assumes all stim_times files are located in
% $SUBJ/scripts/stim_times/ and are called $SUBJ_$REGNAME.1D
regnames = {'let_left' 'let_right' 'ori_left' 'ori_right'};
for rncell = regnames
    rn = rncell{1}; % convert to string
    this_file = ['../' subj_initials '/scripts/stim_times/' subj_initials '_' rn '.1D'];
    st.(rn) = num2cell(textread(this_file,'%n')); % textread will not work in future, and currently this only works because there is a single block/run.  but it works for the current expeirment    
end
    
% convert stimtimes into regressors based on arguments 
% (for pre/post TRs relative to nearest TR).
% Throw warning if TRs don't match up that well with regressors

% which regs do we care about?
conds = fieldnames(st);

% limit ourselves to only those regressors requested
if ~strcmp(options.reg_include,'all')
    % make sure all requested conds are available
    if any(~(ismember(options.reg_include,conds)))
        error('It appears that you have requested a regressor condition that is not defined in this script')
    end
    keep = ismember(conds,options.reg_include);
else
    keep = true(size(conds));
end
if ~strcmp(options.reg_exclude,'none')
    keep = and(keep,~ismember(conds,options.reg_exclude));
end
% keep was based on conds, since that how we input reg_include/exclude, 
% but we need to index stimtimes fields with task in them
conds = conds(keep);

% store a copy of the actual conds used (more specific and used for plotting)
% N.B. that we don't change options.reg_conds, since this would alter our filename (see below)
TC.reg_labels = conds'; % row vector

% TODO - move to top as params
TRs = 2; % TR in seconds
TRn = 176; % number of TRs per run
runs_TR  = repmat(TRs,length(runs),1); % TR (s) for each loaded run
runs_nTR = repmat(TRn,length(runs),1); % n TRs for each loaded run
total_TRs = sum(runs_nTR); % total number of TRs loaded for this expeirment
runs_dur = runs_nTR .* runs_TR; % duration of each run


regs   = zeros(length(conds),total_TRs); % init to all zeros
regs_n = zeros(length(conds),1); % number of repetitions (events/blocks) of each regressor.  will use later to initialize some lists
for c = 1:length(conds)
    this_cond = conds{c};
    
    for r = runs
        % get event/block start times for this run (local times) from stimtimes structure
        local_times = st.(this_cond){r};
        
        % verify that this run included this condition.  if not, local_times will be a string: '*' or '* *'
        if ischar(local_times)
            continue
        end
        
        % adjust based on run order (add sum of all TRs for previous runs)
        global_times = local_times + sum(runs_dur(1:r-1));
        
        % convert to TRs
        % N.B. when converting to TRs, TIME=0 should correspond to TR=1.  So add 1.
        global_TR = global_times / runs_TR(r) + 1;
        if any(mod(global_TR,1))
            max_st_TR_diff = max(mod(global_TR,1));
            if max_st_TR_diff > 0.5 % TODO: make this an optional argument that defaults to 0.5, but can be overridden
                % dont' continue because things seem like they don't line up well
                error('stimtimes onsets do not exactly align with TRs.  max error = %.02f\nmax error exceeds threshold for continuing (need to create optional argument for controlling this)',max_st_TR_diff)
            else
                % just throw warning
                warning('stimtimes onsets do not exactly align with TRs.  max error = %.02f',max_st_TR_diff)
            end
        end
        global_TR_closest = round(global_TR); % find nearest TR (either direction) to use as start of event/block
        
        % update the current regressor to reflect stimtimes (starting TR of each event/block)
        this_regn = length(global_TR_closest); % number of events for this regressor for this run
        regs_n(c,1) = regs_n(c,1) + this_regn; % add to total across all runs
        regs(c,global_TR_closest) = 1; % add to current conditions regressor
    end
end

curreg = 'stimtimes';
subj = init_object(subj,'regressors',curreg);
subj = set_mat(subj,'regressors',curreg,regs);
subj = set_objfield(subj,'regressors',curreg,'condnames',conds);


%% some ROI-independent variables
TC.TRs       = -options.preTRs:options.postTRs-1; % this can also act as an index, relative to block/event onset
uTR = unique(runs_TR);
if length(uTR) > 1
    error('not sure how to set TR for runs with different TR length.  aborting.  check code and decide what to do if this ever happens')
end
TC.time      = TC.TRs * uTR;
TC.nTRs      = length(TC.TRs); % number of TRs per event/block
baseline_ref = ismember(TC.TRs,options.baseline_idx);
if ~all(ismember(options.baseline_idx,TC.TRs))
    error('requested a baseline_idx that includes TRs not spanned by preTRs and postTRs')
end



%% selectors
run_id = [];
for r = 1:length(runs)
    run_id = cat(2,run_id,repmat(r,1,runs_nTR(r)));
end
subj = init_object(subj,'selector','runs');
subj = set_mat(subj,'selector','runs',run_id);
cursel = 'runs';


%% remove "rest"
% to avoid a blocklabels warning, we can create a no_rest selector
subj = create_norest_sel(subj,curreg);
curactsel = [curreg '_norest'];


%% create block labels
% we'll need this to get block/event onsets and offsets
subj = create_blocklabels_custom(subj,curreg,cursel,'actives_selname',curactsel,'keep_contiguous',1,'fixed_contiguous_length',1); % stimtimes reqressor simply denotes event/block start times, so always a length of 1 TR
blocklabels = get_mat(subj,'selector','blocklabels');



%% Save original regressors
orig_curreg = curreg;
TC.regressor_used = orig_curreg;


%% Regressor-Averaged TimeCouse Extraction - loop through all ROIs creating temporary patterns along the way
% setup for voxel selection based on GLM.  the following do not change across ROI, so we should only do them once.
if ~isempty(vox_select)
    switch vox_select{1}
        
        % FEATURE SELECTION BASED ON t-stat threshold
        % N.B. we DON'T need to account for a hemisphere-specific mask for tstat, since each voxel is independently included/excluded
        case 'fullFstat'
            fullFstat_bucket_ref = gt(EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.other.(options.exp_reg)),vox_select{2});
        case 'tstat'
            % positive and negative t-values
            tstat_bucket_ref = gt(EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.t.(options.exp_reg)),vox_select{2});
        case 'tstatabs'
            % positive t-values only
            tstat_bucket_ref = gt(abs(EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.t.(options.exp_reg))),vox_select{2});

        % FEATURE SELECTION BASED ON n highest t-stat
        % N.B. we NEED to account for a hemisphere-specific mask for tstat, since we'll be using rank as our selection criteria
        case 'nvox'
            % positive t-values only
            nvox_t =  EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.t.(options.exp_reg));
        case 'nvoxabs'
            % positive and negative t-values
            nvox_t = abs(EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.t.(options.exp_reg)));
            
        % FEATURE SELECTION BASED ON beta weight threshold
        case 'minbeta'
            % postive beta values only
            beta_bucket_ref = gt(EXP_BUCKET.data(:,:,:,EXP_BUCKET.idx.b.(options.exp_reg)),vox_select{2});

    end
end


% need to initilize these objects so remove_object works for them on the first iteration
for r = 1:length(rois)
    this.roi_name = rois{r};
    if options.verbose; fprintf('\n***************************\nStarting %s\n', this.roi_name); end
    
    
    % create an "cortical area" mask (without regard to hemisphere)
    this.roi_mask = or(ismember(ROI.(rois_brik{r}).data.lh, ROI.(rois_brik{r}).idx.(this.roi_name)), ...
        ismember(ROI.(rois_brik{r}).data.rh, ROI.(rois_brik{r}).idx.(this.roi_name)));
    
    % merge "cortical area" mask with grand_mask, which will account for hemisphere
    this.roi_mask = and(this.roi_mask, ROI.grand_mask);
    
    % check for anatomically-valid voxels
    if ~sum(this.roi_mask(:))
        TC.nvox.(this.roi_name) = 0;
        fprintf('\tCurrent ROI (%s) mask does''t contain any voxels! Skipping',this.roi_name)
        continue
    end
    
    
    % is this voxel suffifiently activated by the EXP contrast (see options.exp_reg)?
    use_anova = 0; % default
    if ~isempty(vox_select)
        switch vox_select{1}
            case {'fullFstat'}
                % FEATURE SELECTION BASED ON Full F-stat threshold
                this.roi_mask = and(this.roi_mask, fullFstat_bucket_ref);
            case {'tstat' 'tstatabs'}
                % FEATURE SELECTION BASED ON t-stat threshold
                this.roi_mask = and(this.roi_mask, tstat_bucket_ref);
            case {'nvox' 'nvoxabs'}
                % FEATURE SELECTION BASED ON n highest t-stat
                
                % TODO: use 'this' structure
                nvox_t_roi = nvox_t(this.roi_mask); % only consider current mask voxels
                nvox_t_roi_sort = sort(nvox_t_roi(:));
                
                % TODO: not sure this works as IPS5 had 101 voxels with 100 vox threshold for e1011.  t-value exactly the same?
                cutoff_idx = min(vox_select{2},length(nvox_t_roi_sort)); % can't take more than all the voxels
                cutoff = nvox_t_roi_sort(end-cutoff_idx+1);
                this.roi_mask = and(this.roi_mask, ge(nvox_t,cutoff));
            case 'minbeta'
                % FEATURE SELECTION BASED ON beta weight threshold
                this.roi_mask = and(this.roi_mask, beta_bucket_ref);
                
            case 'anova'
                % nothing to do at the moment.  this is taken care of below.
                use_anova = 1;
            case 'none'
                % nothing to do, now or ever
            otherwise
                error('first element of vox_select (%s) must be ''fullFstat'', ''tstat'', ''tstatabs'', ''nvox'', ''nvoxabs'', ''anova'' or ''none''',vox_select{1})
        end
    end
    
    % keep track of number of voxels in ROI
    TC.nvox.(this.roi_name) = sum(this.roi_mask(:));
    
    % do we still have voxelss after voxel selection?
    if ~sum(this.roi_mask(:))
        fprintf('\tCurrent ROI (%s) mask does''t contain any voxels after thresholding! Skipping',this.roi_name)
        continue
    end
    
    % create temporary mask in subj structure
    if ~exist_object(subj,'mask','this_mask')
        subj = init_object(subj,'mask','this_mask'); % create it
    end
    subj = set_mat(subj,'mask','this_mask',this.roi_mask); % update
    this_curmask = 'this_mask';
    
    % load in the current ROI into a temporary pattern
    if exist_object(subj,'pattern','this_pattern')
        subj = remove_object(subj,'pattern','this_pattern'); % remove it
    end
    subj = create_pattern_from_mask(subj,curpat,this_curmask,'this_pattern'); % create it
    this_curpat = 'this_pattern';
    
   
    % no-peaking anova
    if use_anova
        alpha = vox_select{2};
        % remove old anova pattern
        if exist_group(subj,'pattern',[this_curpat '_anova'])
            subj = remove_group(subj,'pattern',[this_curpat '_anova']);
        end
        % remove old anova mask
        if exist_group(subj,'mask',[this_curpat '_thresh' num2str(alpha)])
            subj = remove_group(subj,'mask',[this_curpat '_thresh' num2str(alpha)]);
        end
        
        % voxel selection will create a new mask that we will use.  it also creates
        % a new pattern with the ANOVA results, but we won't use that for cross_validation,
        % so we won't update this_curpat.
        statmap_arg.use_mvpa_ver = 1;
        subj = feature_select(subj,this_curpat,curreg,cursel,'thresh',alpha,'statmap_arg',statmap_arg);
        this_curmask = [this_curpat '_thresh' num2str(alpha)];
    end
    
    
    
    % extract pattern/regressor matrices
    this_pat  = get_mat(subj,'pattern',this_curpat);
    this_reg  = get_mat(subj,'regressors',curreg);
    nRegs     = size(this_reg,1);
    this_runs = get_mat(subj,'selector','runs');
    uRuns     = unique(this_runs);
    nRuns     = length(uRuns);
    
    
    % setup roi-specific TC output
    % by run ("raw" timecourse)
    TC.session.runs.(this.roi_name) = runs'; % actual runs string id.  may not match up with uRuns, which is locally defined (ie, always 1:nruns)
    TC.session.nTRs.(this.roi_name) = runs_nTR; % n TRs per run
    TC.session.raw.(this.roi_name)  = repmat(NaN,nRuns,max(runs_nTR)); % averaged time course for entire session (non-normalized), split by runs
    TC.session.base.(this.roi_name) = repmat(NaN,nRuns,1); % run-specific "baseline" (mean)
    TC.session.raw_base.(this.roi_name) = repmat(NaN,nRuns,max(runs_nTR)); % raw - base
    TC.session.raw_norm.(this.roi_name) = repmat(NaN,nRuns,max(runs_nTR)); % 100 * (raw - base)/base (unless we've inputted preprocessing normalized data)
    
    % each block extracted separately
    TC.each.raw.(this.roi_name)      = repmat(NaN,[max(regs_n) TC.nTRs nRegs]); % nReps (max across all regs) X nTRs X nRegs
    TC.each.base.(this.roi_name)     = repmat(NaN,[max(regs_n) 1 nRegs]);       % nReps (max across all regs) X 1    X nRegs
    TC.each.raw_base.(this.roi_name) = repmat(NaN,[max(regs_n) TC.nTRs nRegs]); % nReps (max across all regs) X nTRs X nRegs
    TC.each.raw_norm.(this.roi_name) = repmat(NaN,[max(regs_n) TC.nTRs nRegs]); % nReps (max across all regs) X nTRs X nRegs
    
    % sorted by regressors
    TC.byreg.(this.roi_name)       = repmat(NaN,nRegs,TC.nTRs); % nRegs X nTRs
%     TC.byreg.baseline.(this.roi_name)  = repmat(NaN,nRegs,1);    % nRegs X 1
%     TC.byreg.avg_base.(this.roi_name)  = repmat(NaN,nRegs,TC.nTRs); % nRegs X nTRs
%     TC.byreg.avg_norm.(this.roi_name)  = repmat(NaN,nRegs,TC.nTRs); % nRegs X nTRs
    
    % grand average over all conditions (collapsed across regs)
    TC.grand.(this.roi_name)       = repmat(NaN,1,TC.nTRs); % 1 X nTRs
%     TC.grand.each.(this.roi_name)      = repmat(NaN,[sum(regs_n) TC.nTRs 1]); % nReps (sum across all regs) X nTRs (X 1)
%     TC.grand.baseline.(this.roi_name)  = NaN;                % scalar (1 x 1)
%     TC.grand.avg_base.(this.roi_name)  = repmat(NaN,1,TC.nTRs); % 1 X nTRs
%     TC.grand.avg_norm.(this.roi_name)  = repmat(NaN,1,TC.nTRs); % 1 X nTRs


    % create run-normalized session timecourses
    this_pat_mean = mean(this_pat,1); % mean over all voxels for entire session
    for i = uRuns % uRuns is an index
        this_TRs = 1:TC.session.nTRs.(this.roi_name)(i,1);
        this_ref = this_runs == i;
        this_raw = this_pat_mean(:,this_ref);
        this_base = mean(this_raw);
        
        TC.session.raw.(this.roi_name)(i,this_TRs) = this_raw;
        TC.session.base.(this.roi_name)(i,1) = this_base;
        TC.session.raw_base.(this.roi_name)(i,this_TRs) = this_raw - this_base;
        % normalize to percent signal change
        if strcmp(options.norm_type,'pre')
            % already have normalized data (from pre-processing), just copy
            TC.session.raw_norm.(this.roi_name)(i,this_TRs) = this_raw - this_base; % units are already percent signal change, same as basline subtracted data
        else
            % normalize each event/block by run mean
            TC.session.raw_norm.(this.roi_name)(i,this_TRs) = 100 * (this_raw - this_base) / this_base; % convert to percent signal change
        end
    end

    % loop over regressors and calculate time course
    %grand_bi = 0;
    for i = 1:nRegs
        this_blocklabels = blocklabels .* this_reg(i,:); % ignore blocks that aren't part of this regressor
        byreg_bi = 0;
        for b = unique(this_blocklabels)
            % b==0 represents fixation and/or non-active time points.  skip em'.
            if b ~= 0
                byreg_bi = byreg_bi+1;
                %grand_bi = grand_bi+1;
                idx = find(diff(this_blocklabels == b)==-1); % yeilds regressor onset index (ie, transition in this_blocklabels from 'b' to 'not b'...
                
                % Expand stimtimes to full regressor based on preTRs and postTRs
                reg_start = idx - options.preTRs;
                reg_stop  = idx + options.postTRs - 1;
                reg_idx   = reg_start:reg_stop;
                
                % validate regressors
                if min(reg_idx) < 1
                    error('trying to access TPs (idx %d) that occurred before the start of the experiment (idx 1). try increasing postTRs.',min(reg_idx));
                end
                if max(reg_idx) > size(this_pat,2)
                    error('trying to access TPs (idx %d) that occurred after the end of the experiment (idx %d).  try decreasing postTRs.',max(reg_idx),size(this_pat,2));
                end
                
                TC.each.raw.(this.roi_name)(byreg_bi,:,i) = nanmean(this_pat(:,reg_idx)); % this is or a single event/block, sorted by regressor

                % extract for baseline for this event/block
                if sum(baseline_ref)
                    TC.each.base.(this.roi_name)(byreg_bi,1,i) = mean(TC.each.raw.(this.roi_name)(byreg_bi,baseline_ref,i));
                else
                    TC.each.base.(this.roi_name)(byreg_bi,1,i) = 0;
                end
                this_base = TC.each.base.(this.roi_name)(byreg_bi,1,i);
                TC.each.raw_base.(this.roi_name)(byreg_bi,:,i) = TC.each.raw.(this.roi_name)(byreg_bi,:,i) - this_base;
                
                % normalize to percent signal change
                if strcmp(options.norm_type,'pre')
                    % already have normalized data (from pre-processing), just copy
                    TC.each.raw_norm.(this.roi_name)(byreg_bi,:,i) = TC.each.raw_base.(this.roi_name)(byreg_bi,:,i);
                else
                    % normalize each event/block by preceeding baseline indices
                    TC.each.raw_norm.(this.roi_name)(byreg_bi,:,i) = 100 * (TC.each.raw_base.(this.roi_name)(byreg_bi,:,i) / this_base);
                end
            end
        end
        
        
        % this is where we could filter blocks based on some metric of 'outlierness' (eg, deviation from sum-squared error)
        switch options.block_filter
            case 'sse'
                this_mean = nanmean(TC.each.raw_norm.(this.roi_name)(:,:,i));
                this_diff = TC.each.raw_norm.(this.roi_name)(:,:,i) - repmat(this_mean,max(regs_n),1);
                this_diff2 = this_diff.^2;
                block_sse  = sum(this_diff2,2); % sum-squared error (over entire timecourse) for each block compared to all others in this ROI-Condition
                mean_sse = mean(block_sse);
                std_sse  = std(block_sse);
                thresh   = mean_sse + 2* std_sse; % maximum allowed sum-squared error (2 standard deviations above average SSE)
                block_tofilter = gt(block_sse,thresh); % reference list for blocks (rows) that we will keep
            case 'none'
                block_tofilter = false(max(regs_n),1);
            otherwise
                error('invalid option for block_filter (%s)',block_filter)
        end
        % bookkeeping
        this_nreg = sum(~isnan(TC.each.raw_norm.(this.roi_name)(:,1,i))); % if first value is NaN, then there is no data for that row - this can occur with unbalanced designs
        TC.block_filter.nblocks_orig.(this.roi_name)(1,i) = this_nreg;
        TC.block_filter.nblocks_filtered.(this.roi_name)(1,i) = sum(block_tofilter);
        TC.block_filter.proportion_filtered.(this.roi_name)(1,i) = sum(block_tofilter)/this_nreg;
        % apply filter
        % note, we can't simply remove the rows, else the number of rows per roi/condition gets messed up.
        %       but we can fill those rows with NaNs
        TC.each.raw_norm.(this.roi_name)(block_tofilter,:,i) = NaN;
        
        
        % reg-specific average, error and baseline
        % we will average over raw data, but which one depends on what we want the final output to be (norm_type)
        switch options.norm_type
            case 'each'
                % we have already normalized each event/block individually, now we just average over those
                TC.byreg.(this.roi_name)(i,:) = nanmean(TC.each.raw_norm.(this.roi_name)(:,:,i));
                
                % save multiple error measurements for later plotting
                TC.byreg.err.stdev.(this.roi_name)(i,:) = nanstd(TC.each.raw_norm.(this.roi_name)(:,:,i));
                TC.byreg.err.sem.(this.roi_name)(i,:)   = nanste(TC.each.raw_norm.(this.roi_name)(:,:,i));
                TC.byreg.err.ci95.(this.roi_name)(i,:)  = 1.96 * nanste(TC.each.raw_norm.(this.roi_name)(:,:,i)); % assuming large enough population to use 1.96

            otherwise
                error('no other normalization type except ''each'' is setup yet')
        end
    end
    
    % for grand data, if we have already normalized to 'each' event/block, we can just average over each.raw_norm
    % otherwise, we'll have to average over each.raw, recalculate baseline (and reset each.base), re-normalize raw_norm and re-calc byreg,
    % and finally copy raw_norm.
    switch options.norm_type
        case 'each'
            % extract each block, ignoring regressors, but skipping NaN rows (in case of unequal Ns across regs)
            this_data_byreg = [];
            for ii = 1:nRegs
                this_regdata = TC.each.raw_norm.(this.roi_name)(:,:,ii);
                this_regdata = this_regdata(~isnan(this_regdata(:,1)),:); % avoid averaging in NaN rows
                this_data_byreg = cat(1,this_data_byreg,this_regdata);
            end
            
            % grand average, error and baseline
            TC.grand.(this.roi_name) = nanmean(this_data_byreg);
            
            % save multiple error measurements for later plotting
            TC.grand.err.stdev.(this.roi_name) = nanstd(this_data_byreg);
            TC.grand.err.sem.(this.roi_name)   = nanste(this_data_byreg);
            TC.grand.err.ci95.(this.roi_name)  = 1.96 * nanste(this_data_byreg); % assuming large enough population to use 1.96
            
    end

end % rois

fprintf('_________________________________\n%s Complete\n', mfilename);

%% store TC structure to harddrive
% % outfile = [num2str(subj_initials) '_' options.reg_task]; % uses 'reg_task' instead of 'task' to avoid confusion in stored files
outfile = [num2str(subj_initials)]; % uses 'reg_task' instead of 'task' to avoid confusion in stored files

if ~strcmp(hemisphere,'both')
    outfile = [outfile '_' hemisphere(1) 'h'];
end

if options.smoothing
    outfile = [outfile '_sm' num2str(options.smoothing)];
end

% normalization type
outfile = [outfile '_' options.norm_type 'norm'];

% block filtering
if ~strcmp(options.block_filter,'none')
    outfile = [outfile '_' options.block_filter 'filter'];
end

if ~isempty(vox_select) && ~strcmp(vox_select{1},'none')
    outfile = [outfile '_' vox_select{1} num2str(vox_select{2})];
end

if ~isempty(options.output_tail)
    outfile = [outfile '_' options.output_tail];
end

% remove any '.' from filename so matlab treats it as a .mat file
outfile = ['./stored/' regexprep(outfile,'\.','') '_TC'];

if options.verbose; fprintf('Output in: %s\n', outfile); end

save(outfile,'TC');

return