function bucket = UNR_LoadBucket(subj_initials, bucket_id, options, header_only, subdir)
% UNR_LoadBucket(subj_initials, bucket_id, options, header_only)
%   bucket = UNR_LoadBucket(subj_initials, bucket_id, options, header_only)
%
% Load a bucket file into a structure, including index values for all regressors/contrasts.
%
% Arguments
%   subj_initials = subject's initials (determines directory and filename)
%   bucket_id = string identifier for bucket file. (i.e. '_allruns_norm_bucket').
%                  Excludes subj, and trailing '+orig...'.
%   options = options structure as defined by BrikLoad or a string
%             'matrix' or 'vector' (default) passed into BrikLoad
%   header_only = Boolean, just load header info (no data or errMessage returned)
%                 This will superceed 'options' and it will be ignored.
%   subdir = (optional, default = '') subdirectory of analysis to load from (e.g., 'perms')
%
% Output
%   bucket = bucket structure

% REHBM 10.07.08
%       02.28.09 - optional format argument added
%       06.04.09 - added option for loading header_only
%                - format changed to options structure (see BrikLoad) or string for previous usage (format)
%       02.03.10 - optional argument for loading from a subdirectory of the analysis folder
%       04.2014  - created UNR-version

%% validate arguments
if nargin < 3 || isempty(options)
   options.Format = 'vector'; 
elseif ischar(options)
    f = options;
    
    if ~ismember(f,{'vector','matrix'})
        error('UNR_LoadBucket - for string usage, options (%s) must be ''vector'' or ''matrix''',f);
    end
    
    clear options;
    options.Format = f;    
end

if nargin < 4 | isempty(header_only)
    header_only = 0;
end

if nargin < 5 | isempty(subdir)
    subdir = '';
end


%% load bucket
datadir = fullfile('../',subj_initials,'/analysis/',subdir);
bucketfile = fullfile(datadir,[subj_initials bucket_id '+orig']);
bucket.file = bucketfile; % for easy acess outside of this function

if ~exist([bucketfile '.HEAD'],'file')
    error('bucket file (%s) doesn''t exist',bucketfile)
end

if header_only
    % load bucket from disk
    [bucket.err bucket.info] = BrikInfo(bucketfile);
    if bucket.err
        bucket.errMessage = 'Error reading header info';
    else
        bucket.errMessage = '';
    end
else
    % load bucket from disk
    [bucket.err bucket.data bucket.info bucket.errMessage] = BrikLoad(bucketfile,options);
end

% check for errors
if bucket.err
    error(bucket.errMessage)
end

%define which indices of BRIK files correspond to conditions and ROIs in
%experiments. Take AFNI index+1 (no 0 in MATLAB).  We can read this directly from the bucket info.
brik_labels = textscan(bucket.info.BRICK_LABS,'%s','Delimiter','~'); % a nicer format
brik_labels = brik_labels{1};

%% parse Brik Labels from bucket.info for easy indexing.
% Sort by Coef, Tstats and Fstats.
% init bucket.idx structure
for fncell = {'b' 't' 'f' 'other'}
    bucket.idx.(fncell{1}) = [];
end
bucket_idx = 0; % we need to use a counter in case we also included options.Frames.  BRIK_LABS doesn't account for that filtering.
for i = 1:length(brik_labels)
    % was this part of the requested list of Frames?
    if isfield(options,'Frames') & ~ismember(options.Frames,i)
        continue; % user doesn't want to load this
    end
    bucket_idx = bucket_idx+1;
    
    % determine how to file the index in the bucket.idx structure
    if regexp(brik_labels{i},'Coef$')
        % this is a beta coefficient label
        %%%lab = brik_labels{i}(1:(end-7)); % FAILS WHEN NUMBER OF COEF IN FIR MODEL EXCEEDS 9
        lab = regexprep(brik_labels{i},'#.+_Coef$','');
        fn  = 'b';
    elseif regexp(brik_labels{i},'Tstat$')
        % this is a t-statistic label
        %%%lab = brik_labels{i}(1:(end-8));
        lab = regexprep(brik_labels{i},'#.+_Tstat$','');
        fn  = 't';
    elseif ~isempty(regexp(brik_labels{i},'Fstat$')) && ~ismember(brik_labels{i},{'Full_Fstat'})
        % this is a f-statistic label (newer)
        %%%lab = brik_labels{i}(1:(end-8));
        lab = regexprep(brik_labels{i},'#.+_Fstat$','');
        fn  = 'f';
    elseif ~isempty(regexp(brik_labels{i},'F-stat$')) && ~ismember(brik_labels{i},{'Full F-stat'})
        % this is a f-statistic label (newer)
        %%%lab = brik_labels{i}(1:(end-8));
        lab = regexprep(brik_labels{i},' F-stat$','');
        fn  = 'f';
    else
        % unknown type, file as 'other'
        lab = brik_labels{i};
        fn  = 'other';
    end

    % some cosmetic formatting of label
    % trim _GLT
    if regexp(lab,'_GLT$')
        lab = lab(1:(end-4));
    end
    
    % remove '-', '#' and ' ' since they are invalid Matlab fieldname characters
    lab = regexprep(lab,'[\[\] #-]','_');
    lab = strtrim(lab); % remove leading and trailing whitespace

    % file the index appropriately
    if ~isfield(bucket.idx.(fn),lab)
        % init, just add it
        bucket.idx.(fn).(lab) = bucket_idx;
    else
        % we already have one of these, so it is probably an FIR model with multiple
        % coefficients for each regressor.  append this one to the end of current list.
        bucket.idx.(fn).(lab)(end+1) = bucket_idx;
    end
      
        
end
bucket.nGLMs = bucket_idx; % total number of LOADED GLM entries in bucket (regressors+contrasts, betas and stats)

