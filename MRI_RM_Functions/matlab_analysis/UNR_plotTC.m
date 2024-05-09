function [dataout] = plotTC(TC,rois,sorting,plottype,reg_include,errbars,figID)
% plotTC
%   [dataout] = plotTC(TC,[rois],[sorting],[plottype],[reg_include],[errbars],[figID])
%
% Quick plots for TimeCourse analysis (rawTC)
%
% Arguments
%   TC = TimeCourse structure returned by rawTC.m
%   rois = 'visual', 'parietal',  'all', or 'avail' (default), or cell array with list of specific rois
%   sorting = 'byreg' or 'grand' (default) to collapse across regressors
%              N.B. 'grand' averaging is done in this code so we can apply reg_include.  not just read from TC structure.
%   plottype = {plottype plotargs}
%              'timecourse' or 'tc' for timecourse plots for each roi
%              'nvox', 'nvox_anat' or 'nvox_prop' for number of voxels
%              nvox into proportion of voxels of total ROI)
%              {'mean' [TRs]} for calculating means over peak timepoints
%              {'index' [TRs]} for calculating standard index (A-B)/(A+B)
%              {'dprime' [TRs]} for calculating dprime index (A-B)/sqrt((var(A)+var(B))/2)
%              {'diff' [TRs]} for calculating (A-B)
%                   for 'mean', 'index', 'dprime' and 'diff' plotargs defines what TRs (see TC.TRs) to use extract timecourse data from
%   reg_include = cell array of regressor labels (strings) or regressor pairs (cell array of strings) to include in plot.
%                 default is 'all' for all regressors individually.  Pairs are required for plottype{1} == 'index', 'dprime' or 'diff'
%   errbars = 'stdev' 'sem' 'ci95' (default) or 'none' .  automatically set to 'none' if 'reg_include' is not 'all' since we can't do voxel-wise error bars at this stage
%   figID = if not empty, supplies a figure handle to draw the resulting plots to


% REHBM 03.05.09 - created, following plotTC.m
%       03.21.11 - re-worked using plotFIR as template.  updated for new plotTC.
%       10.21.11 - changed from reference- to index-based regressor extraction (for reg_include) to maintain order of reg_include input for plotting with predictable colors
%       11.08.11 - added dprime and diffplottype (integrated with index), altered the way regressors are input for index and dprime plottype (now input as cell array pairs)

%% validate arguments
if nargin < 2 || isempty(rois)
    rois = 'avail';
end

if nargin < 3 || isempty(sorting)
    sorting = 'grand';
end

if nargin < 4 || isempty(plottype)
    plottype = 'tc';
end
if ~iscell(plottype)
    plottype = {plottype};
end
if length(plottype) > 1
    plotargs = plottype{2};
end
plottype = plottype{1}; % back to string

% check to see if we requested 'session' summary plots, in which case, plottype AND sorting should be 'session'
if strcmp(plottype,'session') || strcmp(sorting,'session')
    sorting  = 'session';
    plottype = 'session';
end

if nargin < 5 || isempty(reg_include)
    reg_include = 'all';
end

if nargin < 6 || isempty(errbars)
    errbars = 'ci95';
end

if nargin < 7 || isempty(figID)
    % no figure ID supplied, open a new one
    %figID = figure('MenuBar','none','OuterPosition',[440 358 448 394]); % defualt OuterPosition = [440   358   560   493]
    figID = figure('MenuBar','none','Color','w','PaperOrientation','landscape');
end
figure(figID); % focus figure to draw

% expand rois argument
if strcmp(rois,'avail')
    rois = TC.rois;
else
    rois = UNR_ROIArgs(rois);
end


% create reg selector list - and validate regressor structure while we're at it!
if strcmp(reg_include,'all')
    regidx = 1:length(TC.reg_labels);
elseif ischar(reg_include)
    % only one regressor defined
    regidx = 1;
elseif ischar(reg_include{1})
    % cell array of strings (list of individual ROIs)
    %errbars = 'none'; % because error bars have been calculated over all regs in rawTC
    regidx = repmat(NaN,size(reg_include));
    for i = 1:length(reg_include)
        regidx(i) = find(ismember(TC.reg_labels,reg_include{i}));
    end
else
    % we probably have a cell array of cell arrays of strings (list of ROI pairs)
    % N.B. we expect PAIRs of ROIs, so length of each sub-cell should be EXACTLY 2
    if ~all(cellfun(@length,reg_include)==2)
        error('when specifying ROI pairs in reg_include, each pair (cell array) must contain exactly 2 ROIs')
    end
    
    % pairs are only valid for certain plottypes
    if ~ismember(plottype,{'index' 'dprime' 'diff'})
        error('reg_include ROI pairs are not valid for plottype of %s',plottype)
    end
    
    regidx = repmat({[NaN NaN]},size(reg_include));
    for i = 1:length(reg_include)
        regidx{i}(1) = find(ismember(TC.reg_labels,reg_include{i}(1))); % first of pair
        regidx{i}(2) = find(ismember(TC.reg_labels,reg_include{i}(2))); % second of pair
    end
end


% update vox_select in TC structure for marking plots
if isempty(TC.vox_select) || strcmp(TC.vox_select{1},'none')
    TC.vox_select = {'' ''};
end


% init output, just in case
dataout.subj_initials = TC.subj_initials;
dataout.rois = rois;
dataout.sorting = sorting;
dataout.plottype = plottype;
% data.reg_labels will account for reg_include, AND reorders them to match reg_include input (or defaults to that of data structure for 'all')
if iscell(regidx)
    for i = 1:length(regidx)
        % convert to a single string for each pair
        this_pair = TC.reg_labels(regidx{i});
        dataout.reg_labels{i} = [this_pair{1} '-' this_pair{2}];
    end
else
    dataout.reg_labels = TC.reg_labels(regidx);
end
dataout.errbars = errbars;
dataout.TPs = TC.TRs; % get directly from TC structure (unlike FIR)
dataout.time = TC.time;
 

%% plot it

nROIs = length(rois);

% setup paneling in case we need it
nrow = ceil(sqrt(nROIs));
ncol = floor(sqrt(nROIs));
if nrow*ncol < nROIs
    nrow = nrow + 1;
end


switch plottype
    case {'tc' 'timecourse'}
        % for axis matching
        ymin = Inf;
        ymax = -Inf;

        for r = 1:nROIs
            roi = rois{r};
            if strcmp(roi,'null')
                % skip, no labeling
                continue
            end
            
            subplot(nrow,ncol,r)
            hold on;
            
            plotsymb = 'o-';%plotsymb = 'r--';
        
            
            if isfield(TC.(sorting),roi)
                % mark baseline measurement
                if ~isempty(TC.options.baseline_idx)
                    baseline_ref = ismember(TC.TRs,TC.options.baseline_idx);
                    baseline_xtime = TC.time(baseline_ref);
                    if length(baseline_xtime) == 1
                        % padd it a little so it is more obvious
                        pad = 0.25 * (TC.time(2)-TC.time(1));
                        baseline_xtime = [baseline_xtime baseline_xtime] + [-pad pad];
                    end
                    plot(baseline_xtime,zeros(length(baseline_xtime),1),'Color',[.75 .75 .75],'LineWidth',2); %plot(TC.TRs(baseline_ref),zeros(sum(baseline_ref),1),'Color',[.75 .75 .75],'LineWidth',2);
                end
            
                % extract the data we care about
                this_data = TC.byreg.(roi)(regidx,:);
                if strcmp(sorting,'grand')
                    % we'll need to recalculate grand average and variance estimates after filtering out the regressors requested in the arguments
                    this_data_byreg = [];
                    for ii = regidx
                        this_regdata = TC.each.raw_norm.(roi)(:,:,ii);
                        this_regdata = this_regdata(~isnan(this_regdata(:,1)),:); % avoid averaging in NaN rows
                        this_data_byreg = cat(1,this_data_byreg,this_regdata);
                    end
                    this_data_byreg = this_data; % for errbars
                    this_data = mean(this_data,1);
                end
                
                % the main data
                plot(TC.time,this_data',plotsymb); %plot(TC.TRs,this_data',plotsymb);
                
                % error bars
                if ~strcmp(errbars,'none')
                    if ismember(errbars,fieldnames(TC.(sorting).err))
                        if strcmp(sorting,'grand')
                            switch errbars
                                case 'stdev'
                                    this_err = std(this_data_byreg,[],1);
                                case 'sem'
                                    this_err = ste(this_data_byreg,1);
                                case 'ci95'
                                    this_err = ci(this_data_byreg);%1.96 * ste(this_data_byreg,1);
                                otherwise
                                    error('unrecognized errbars value (%s)',errbars)
                            end
                        else
                            this_err = TC.byreg.err.(errbars).(roi)(regidx,:);
                        end
                        if numel(TC.TRs) == 1
                            % then we need to fake the multiple regressors by doubling our datapoints
                            % or all datapoints are treated as a single regressor
                            errorbar(repmat(TC.time,size(this_data,1),2)',repmat(this_data,1,2)',repmat(this_err,1,2)',plotsymb); %errorbar(repmat(TC.TRs,size(this_data,1),2)',repmat(this_data,1,2)',repmat(this_err,1,2)',plotsymb)
                        else
                            % plot it at once
                            errorbar(repmat(TC.time,size(this_data,1),1)',this_data',this_err',plotsymb); %errorbar(repmat(TC.TRs,size(this_data,1),1)',this_data',this_err',plotsymb)
                        end
                    else
                        error('TC.(%s).err does not define requested errbars type (%s)',sorting,errbars)
                    end
                end
                
                % figure out a good axis size
                minx = min(TC.time)-1; %minx = min(TC.TRs)-1;
                maxx = max(TC.time)+1; %maxx = max(TC.TRs)+1;
                xlim([minx maxx]);
                
                
                if ~strcmp(errbars,'none')
                    miny = min(min(this_data-this_err));
                    maxy = max(max(this_data+this_err));
                else
                    miny = min(min(this_data));
                    maxy = max(max(this_data));
                end
                ypadpc = 0.01* (maxy-miny);
                miny = miny-ypadpc;
                maxy = maxy+ypadpc;
                
                ymin = min(miny,ymin);
                ymax = max(maxy,ymax);
                
                set(gca,'YTick',[-10:.5:10 90:0.5:110]);
            end
            
            % legend on last plot
            if r == nROIs && strcmp(sorting,'byreg')
                this_leg = regexprep(TC.reg_labels(regidx),'_','\\_');
                if ~isempty(TC.options.baseline_idx)
                    this_leg = [{'baseline'} this_leg];
                end
                legend(this_leg,'Location','NorthWest');
            end
            
            % labels and title
            title(roi)
            xlabel('time (s) relative to block/event onset')% xlabel(sprintf('TRs (%fs) from Block/Event onset',sdb.TR/1000))
            ylabel('% Signal Change (-baseline)');
            vline(0,'k--'); % event/block onset
            
            hold off;
        end
        
        % set axes to match
        for r = 1:nROIs
            subplot(nrow,ncol,r);
            if ymin~=ymax && ~isinf(ymin) && ~isinf(ymax) && ~isnan(ymin) && ~isnan(ymax) % else there is likely no data to scale
                ylim([ymin ymax]);
            end
        end
        
        % output the timecourses
        dataout.(plottype).tc = this_data;
        
    case {'session'}
        
        close(figID); % we'll open our own for session plots
        
        for r = 1:nROIs
            roi = rois{r};
            if strcmp(roi,'null')
                % skip, no labeling
                continue
            end
            
            
            figure('MenuBar','none','Color','w','PaperOrientation','landscape');
            % 1. full session TC (with vertical lines denoting run separations)
            subplot(2,1,1)
            hold on;
            plotsymb = '-';%plotsymb = 'r--';

            if isfield(TC.session.raw_norm,roi)
                % extract the data we care about
                this_full = TC.session.raw_norm.(roi)';
                this_full = this_full(:)';
                this_full = this_full(~isnan(this_full)); % NaNs occur "between" runs when the number of TRs is not the same across every run (ie, NaNs are end-padding in the matrix format)

                % the main data
                plot(this_full,plotsymb);
                
                % mark the start of each run
                runends = cumsum(TC.session.nTRs.(roi)); % last TR of each run
                runends = runends(1:end-1); % no need to mark last run's endtime
                runends = runends + 0.5; % shift vertical lines by half a TR so the division occurs *between* runs
                vline(runends,'k:'); % mark with dotted black vertical lines
                
                title('run-averaged timecourse');
                xlabel('TR')
                ylabel('percent signal change');
                xlim([-1 length(this_full)+1]);
            end
            
            % 2. run-averaged TC (warning: may not make sense for some designs, like event-related or highly randomized designs)
            subplot(2,1,2)
            hold on;
            plotsymb = '-';
            
            if isfield(TC.session.raw_norm,roi)
                % extract the data we care about
                this_runmean = TC.session.raw_norm.(roi);
                this_runmean = nanmean(this_runmean); % average over runs (ignore NaN end-padding)
                
                % the main data
                plot(this_runmean,plotsymb);
                
                % check for consistency of nTRs over runs
                if length(unique(TC.session.nTRs.(roi))) > 1
                    warning('number of TRs is not consistent across runs.  interpret run-averaged timecourse with caution')
                    title('run-averaged timecourse (caution: nTRs not consistent across runs!)');
                else
                    title('run-averaged timecourse');
                end
                xlabel('TR')
                ylabel('percent signal change');
                xlim([-1 length(this_runmean)+1]);
            end
            
            % mark plot
            MarkPlot(sprintf('plotTC +/- %s  -  %s  %s  %s  %s  %s  %s  %s', errbars, mat2str(TC.subj_initials), TC.task, TC.hemisphere,[TC.vox_select{1} num2str(TC.vox_select{2})],TC.options.output_tail,roi));

            
            % 3. run-separated TC plots (in new figure)
            figure('MenuBar','none','Color','w','PaperOrientation','landscape');
            
            % setup paneling in case we need it
            nRuns = size(TC.session.raw_norm.(roi),1);
            nrow = ceil(sqrt(nRuns));
            ncol = floor(sqrt(nRuns));
            if nrow*ncol < nRuns
                nrow = nrow + 1;
            end
            
            if isfield(TC.session.raw_norm,roi)
                % extract the data we care about
                this_byruns = TC.session.raw_norm.(roi);

                % for axis matching
                ymin = Inf;
                ymax = -Inf;
                
                % loop over runs
                for i = 1:nRuns
                    subplot(nrow,ncol,i)
                    
                    this_run = TC.session.runs.(roi){i}; % run string, i is the matching index
                    this_rundata = this_byruns(i,~isnan(this_byruns(i,:))); % ignore NaN end-padding
                    
                    % the main data
                    plot(this_rundata,plotsymb);
                
                    title(this_run);
                    %xlabel('TR')
                    %ylabel('percent signal change');
                    
                    
                    % figure out a good axis size
                    xlim([-1 length(this_rundata)+1]);
                    
                    miny = min(this_rundata);
                    maxy = max(this_rundata);
                    ypadpc = 0.01 * (maxy-miny);
                    miny = miny-ypadpc;
                    maxy = maxy+ypadpc;
                    
                    ymax = max([ymax abs(miny) abs(maxy)]);        %max(maxy,ymax);
                    ymin = -max([ymax abs(miny) abs(maxy)]); %min(miny,ymin);
                end
                
                
                % set axes to match
                for i = 1:nRuns
                    subplot(nrow,ncol,i);
                    if ymin~=ymax && ~isinf(ymin) && ~isinf(ymax) && ~isnan(ymin) && ~isnan(ymax) % else there is likely no data to scale
                        ylim([ymin ymax]);
                    end
                end
            end
            
            % mark plot
            MarkPlot(sprintf('plotTC +/- %s  -  %s  %s  %s  %s  %s  %s  %s', errbars, mat2str(TC.subj_initials), TC.task, TC.hemisphere,[TC.vox_select{1} num2str(TC.vox_select{2})],TC.options.output_tail,roi));

            
            
            hold off;
        end
        
        % output the timecourses
        dataout.(plottype).full    = this_full;
        dataout.(plottype).runmean = this_runmean;
        dataout.(plottype).run     = this_byruns;
        
    case 'mean'
        % extract mean activation over some timepoints (TRidx) defined by plotargs
        
        % extract plotargs
        peakTRs = plotargs;
        
        % validate arguments
        if ~all(ismember(peakTRs,TC.TRs))
            error('invalid TRs supplied for peak mean extraction.  Must be from set defined in TC.TRs')
        end
        TRidx = find(ismember(TC.TRs,peakTRs)); % convert to idx
        
        m = repmat(NaN,nROIs,length(regidx));
        sem  = repmat(NaN,nROIs,length(regidx)); % SEM over extracted timepoints
        ci95 = repmat(NaN,nROIs,length(regidx)); % 95% confidence interval over extracted timepoints
        %p  = repmat(NaN,nROIs,length(regidx));
        %sig = zeros(nROIs,length(regidx)); % sig by two tailed permutation test?
        for r = 1:nROIs
            roi = rois{r};
            
            if isfield(TC.byreg,roi)
                % extract the data we care about
                this_data = TC.byreg.(roi)(regidx,:); % full time course for both regs
                this_peak_data = this_data(:,TRidx); % extract peak TC.TRs
                
                m(r,:)   = mean(this_peak_data,2);% average over extracted timepoints
                sem(r,:) = ste(this_peak_data,2);% standard error over extracted timepoints
                ci95(r,:)  = ci(this_peak_data,[],2); %1.96 * sem(r,:); % 95% confidence interval (estimated - assumes large N)
            end
        end
        
        % ignore error bars for now since they would be over timepoints.  a separate script should loop over subjects and re-calculate error bars across subjects.
        e_bars = ci95; % copy an error matrix, or []
        %         if size(m,1) == 1
        %             % then we are only plotting one ROI, so we need to pad m to make sure things are grouped appropriately
        %             barerr2(rois', [m;repmat(NaN,size(m))], [e_bars;repmat(NaN,size(m))], 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'});
        %         else
        barerr2(rois', m, e_bars, 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'});
        %         end
        
        % labels and title
        title(mat2str(TC.subj_initials));
        ylabel(sprintf('mean activation'));
        this_leg = regexprep(dataout.reg_labels,'_','\\_');
        legend(this_leg)
        
        % output
        dataout.(plottype).mean = m;
        %         dataout.(plottype).e = e;
        
    case {'index' 'dprime' 'diff'}
        % calculate standard index or dprime index over some timepoints (TRidx) defined by plotargs
        %    standard index = (A+B) / (A-B)
        %    dprime index   = (A-B) / sqrt((var(A)+var(b))/2)
        %       refs: Pinsk et al, 2009; Afraz et al, 2006; Grill-Spector et al, 2006erratum
        %    diff           = A-B
        
        % extract plotargs
        peakTRs = plotargs;
        
        % validate arguments
        if ~all(ismember(peakTRs,TC.TRs))
            error('invalid TRs supplied for %s calculation.  Must be from set defined in TC.TRs',plottype)
        end
        if strcmp(reg_include,'all')
            error('%s calculation requires reg_include to be explicitly defined',plottype)
        end
        if ~strcmp(sorting,'byreg') || ~iscell(regidx) % already verified that a sub-cells are of length 2 during argument validation
            error('for %s calculation, you must use ''byreg'' sorting and include a list of regressor pairs',plottype)
        end
        
        TRidx = find(ismember(TC.TRs,peakTRs)); % convert to idx
        
        metric = repmat(NaN,nROIs,length(regidx));
        %e  = repmat(NaN,nROIs,1);
        %ci95 = repmat(NaN,nROIs,2); % 95% confidence interval, pulled from permutation test
        %p  = repmat(NaN,nROIs,1);
        %sig = zeros(nROIs,1); % sig by two tailed permutation test?
        for r = 1:nROIs
            roi = rois{r};
            
            if isfield(TC.byreg,roi)
                % loop over all requested ROI-pairs
                for k = 1:length(regidx)
                    this_regidx = regidx{k};
                    
                    % extract the data we care about
                    this_data = TC.byreg.(roi)(this_regidx,:); % full time course for both regs
                    peak_data = this_data(:,TRidx); % extract peak TC.TRs
                    peak_means = mean(peak_data,2); % mean over peak TC.TRs
                                        
                    switch plottype
                        case 'index'
                            if min(peak_means) < 0
                                warning('\n\tsome ''peak'' values are negative for %s. index will be invalid',roi)
                            end
                            % by using regidx, peak_means will inherently be in the order requested by reg_labels
                            % STANDARD INDEX = A-B/A+B
                            metric(r,k) = diff(peak_means(end:-1:1)) / sum(peak_means); 
                        case 'dprime'
                            % by using regidx, peak_means will inherently be in the order requested by reg_labels
                            
                            % for dprime, we need to calculate variance (of peak_mean) over blocks
                            % N.B. variance [== sqrt(standard deviation)] of mean over blocks
                            block_data = TC.each.raw_norm.(roi)(:,TRidx,this_regidx); % extract relevant timepoints for each block and condition
                            peak_vars  = squeeze(var(mean(block_data,2))); % extract mean of each block, calculate variance across blocks, collapse to 2 values
                            
                            % DPRIME = A-B/sqrt([var(A)-var(B)]/2)
                            metric(r,k) = diff(peak_means(end:-1:1)) / sqrt(sum(peak_vars)/2);
                        case 'diff'
                            % simple difference
                            % DIFF = A-B
                            metric(r,k) = diff(peak_means(end:-1:1));
                    end
                end
            end
        end
        
        barerr2(rois', metric, [], 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'});
        
        % labels and title
        title(mat2str(TC.subj_initials));
        switch plottype
            case 'index'
                ylabel('standard index [ (A-B) / (A+B) ]');
            case 'dprime'
                ylabel('d'' index [ (A-B) / sqrt((var(A)+var(B)/2) ]');
            case 'diff'
                ylabel('difference [ A-B ]');
        end
        
        this_leg = regexprep(dataout.reg_labels,'_','\\_');
        legend(this_leg)
        
        %%%miny = -1;
        %%%maxy = 1;
        %ylim([miny maxy]);
        
        % output
        dataout.(plottype).metric = metric;
        %         dataout.(plottype).e = e;
              
    case {'nvox' 'nvox_anat' 'nvox_prop'}
        n  = repmat(NaN,nROIs,1);
        for r = 1:nROIs
            roi = rois{r};
            
            if isfield(TC.(plottype),roi)
                % actual data
                n(r,1) = TC.(plottype).(roi)(end); % the last value in nvox* is the average of multiple datasets if this TC struct is an average of multiple datasets
            end
        end
        
        barerr2(rois', n, [], 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'});
        
        % labels and title
        title(mat2str(TC.subj_initials));
        ylabel(sprintf('%s',plottype));
        %ylim([miny maxy]);
        
        % output
        dataout.(plottype).n = n;

    otherwise
        close(figID); % closes last figure opened
        error('invalid plottype (%s). not defined in plotTC.m',plottype)
end


%% final figure stuff
if ~strcmp(sorting,'session')
    MarkPlot(sprintf('plotTC +/- %s  -  %s  %s  %s  %s  %s', errbars, TC.subj_initials, TC.hemisphere,[TC.vox_select{1} num2str(TC.vox_select{2})],TC.options.output_tail));
end
