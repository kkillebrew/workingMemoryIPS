function [TC] = UNR_average_rawTC(all_subj_initials,hemisphere,vox_select,varargin)
% average_rawTC
%   [TC] = average_rawTC(expid,hemisphere,vox_select,OPTIONAL_ARGS)
%
% Load in rawTC analysis results (TC strctures) from multiple all_subj_initials, concat them together and
% re-run analysis.
%
% Arguments
%   See UNR_rawTC for argument explaination.
%
% Return
%   TC    = average representation similarity structure across all_subj_initials

% Hidden options
% suffix [''] = string id of perm bucket for filename.  Used internally, not meant for direct manipulation.
% nosave [0] = boolean, should we suppress saving the file to disk?  unnecessary for permutations

% REHMB 03.11 - created as modification of average_glmFIR for rawTC

% TODO


%% validate arguments

% hemisphere
if nargin < 2
    hemisphere = 'both';
end
if ~ismember(hemisphere,{'both' 'left' 'right'})
    error('invalid hemisphere (%s).  Must be ''left'', ''right'', or ''both''', hemisphere)
end

% voxel selection
if nargin < 3 || isempty(vox_select)
    vox_select = {'none'};
end
if ~iscell(vox_select)
    vox_select = {vox_select};
end

defaults.preTRs              = 2;
defaults.postTRs             = 16;
defaults.baseline_idx        = NaN; % place-holder.  if undefined, will use -preTRs:0
defaults.norm_type           = 'each';
defaults.datafile_prefix     = 'tscvrdtm';
defaults.bucket_prefix       = 'tscvrsm4_norm_bucket';
defaults.exp_reg             = 'EXP';
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


%% determine filename for these options
%   (taken from final output of rawTC)
file_suffix = '';

if strcmp(hemisphere,'both')
    file_suffix = [file_suffix '_both'];
else
    file_suffix = [file_suffix '_' hemisphere(1) 'h'];
end

% normalization type
file_suffix = [file_suffix '_' options.norm_type 'norm'];

% block filtering
if ~strcmp(options.block_filter,'none')
    file_suffix = [file_suffix '_' options.block_filter 'filter'];
end

if ~isempty(vox_select) && ~strcmp(vox_select{1},'none')
    file_suffix = [file_suffix '_' vox_select{1} num2str(vox_select{2})];
end

if ~isempty(options.output_tail)
    file_suffix = [file_suffix '_' options.output_tail];
end

% remove any '.' from filename so matlab treats it as a .mat file
file_suffix = [regexprep(file_suffix,'\.','') '_TC'];



%% load in subj files with requested options
for i = 1:length(all_subj_initials)
    thisfile = fullfile('./stored',[all_subj_initials{i} file_suffix]);
    tmp = load(thisfile); % brings in TC structure
        
    % check options to make sure they match arguments
    if ~isequal(options,tmp.TC.options)
        disp(options)
        disp(tmp.TC.options)
        error('options are not consistent between arguments and expid %d.  maybe you aren''t including an optional argument in the average_rawTC call that was included in the individual datasets?',all_subj_initials(i))
    end
    
    % storage
    if i==1
        % start with this one, much of this will be overwritten
        TC = tmp.TC;
        fns = {'subj_initials' 'runs'};
        for k = 1:length(fns)
            TC.(fns{k}) = {};
        end
        TC.missing = {};
    end
    
    
    % pack some data into cell form
    TC.subj_initials{i} = tmp.TC.subj_initials;
    fns = {'subj_initials' 'runs'};
    for k = 1:length(fns)
        TC.(fns{k}){i} = tmp.TC.(fns{k});
    end

    
    % average timecourse with previously loaded one
    for sortingcell = {'byreg' 'grand'}
        this_sorting = sortingcell{1};
        
        rois = TC.rois;
        for ri = 1:length(rois)
            this_roi = rois{ri};
            if isfield(tmp.TC.(this_sorting),this_roi)
                % timcourse
                TC.(this_sorting).(this_roi)(:,:,i) = tmp.TC.(this_sorting).(this_roi);
                
                %                 % perms
                %                 if TC.options.nperms
                %                     TC.(this_sorting).perm.dist_square.(this_roi)(:,:,:,i) = tmp.TC.(this_sorting).perm.dist_square.(this_roi);
                %                 end
                
                % nvox
                TC.nvox.(this_roi)(:,i) = tmp.TC.nvox.(this_roi);
                TC.nvox_anat.(this_roi)(:,i) = tmp.TC.nvox_anat.(this_roi);
                TC.nvox_prop.(this_roi)(:,i) = tmp.TC.nvox_prop.(this_roi);
            else
                % not to screen and in TC struct
                fprintf('\n************************************************\n')
                warning('ROI %s not found for subj %s.  Skipping...',this_roi,all_subj_initials{i})
                TC.missing{end+1} = [num2str(all_subj_initials{i}) '_' this_roi];
            end
            
        end
    end
end



%% Average time courses and update variance estimates
for sortingcell = {'byreg' 'grand'}
    this_sorting = sortingcell{1};
    
    %     rois = fieldnames(TC.(this_sorting).dist_square);
    rois = TC.rois;

    for ri = 1:length(rois)
        this_roi = rois{ri};

        if isfield(TC.(this_sorting).err.stdev,this_roi)
            % udpate associated err measurements first, before we overwrite with just the mean
            TC.(this_sorting).err.stdev.(this_roi) = std(TC.(this_sorting).(this_roi),[],3);
            TC.(this_sorting).err.sem.(this_roi)   = ste(TC.(this_sorting).(this_roi),3);
            TC.(this_sorting).err.ci95.(this_roi)  = 1.96 * ste(TC.(this_sorting).(this_roi),3); % assuming large enough population to use 1.96
            %TC.(this_sorting).err.ci95.(this_roi)  = ci(TC.(this_sorting).(this_roi),95,3); % assuming large enough population to use 1.96
                
            % time course
            TC.(this_sorting).(this_roi) = mean(TC.(this_sorting).(this_roi),3);
            
            
            %         % perm
            %         if TC.options.nperms
            %             TC.(this_sorting).perm.dist_square.(this_roi) = mean(TC.(this_sorting).perm.dist_square.(this_roi),4);
            %         end
            % nvox
            if strcmp(this_sorting,'byreg') % make sure we only do this once
                TC.nvox.(this_roi)(:,end+1) = mean(TC.nvox.(this_roi));
                TC.nvox_anat.(this_roi)(:,end+1) = mean(TC.nvox_anat.(this_roi));
                TC.nvox_prop.(this_roi)(:,end+1) = mean(TC.nvox_prop.(this_roi));
            end
        end
    end
    % udpate NOTES on err bars
    TC.(this_sorting).NOTES = 'err measurements are over all_subj_initials';
end


if options.verbose; fprintf('_________________________________\n%s Complete\n', mfilename); end


TC.subj_initials = cell2mat(all_subj_initials); % hack so that UNR_plotTC will work with averaged data

%% store TC structure to harddrive
expid_tag = cell2mat(all_subj_initials); %mat2str(all_subj_initials);
expid_tag = regexprep(expid_tag,'[','('); % [ -> (
expid_tag = regexprep(expid_tag,']',')'); % ] -> )
expid_tag = regexprep(expid_tag,' ','_'); % <space> -> _
outfile = fullfile('./stored',[expid_tag file_suffix]);

save(outfile,'TC');
