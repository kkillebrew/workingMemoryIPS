%%%function [GROUP] = group_TC(subj_initials,suffix,rois,sorting,reg_include,plottype,varargin)


% group_TC_script
%   To use this script, first you need to set up an options structure that defines the arguments listed below.  Then
%   you simply call group_TC_script.  This will allow you to maintain all of the internal variables in the current workspace,
%   but also keeps the main body of this script in a central location for use with multiple experiemnts.
%
% Required Fields of G (for "group") strucure (e.g., G.subj_initials = {'4101' '4102'};)
%     subj_initials      = cell array of experimental ids. First part of stored rawTC filenames.
%     suffix      = suffix of the stored rawTC filenames.
%     rois        = cell array of rois to include, or 'all' for all rois that are part of the loaded files.
%     sorting     = 'byreg' or 'grand' (passed into UNR_plotTC).
%     reg_include = cell array of regressors (or regressor pairs) to include in the analysis (passed into UNR_plotTC).
%     plottype    = passed into UNR_plotTC, see UNR_plotTC for options.  For 'tc' and 'session' there
%                   is no actual group analysis to perform.  This script is just a shortcut for
%                   batch plotting and printing.
% Optional Fields of 'options' structure [default]
%     stattype ['rm_anova'] = 'rm_anova', 'friedman' or 'none' to skip statistical tests
%     posthooc_test ['hsd'] = passed to multcompare 'ctype'.  but note that multcompare does not take
%                             into account random-effects, so those provided here are for fixed-effects
%                             and aren't optimal for within-subjects designs.
%     errbars ['none']      = passed into UNR_plotTC.
%     plot_group [1]        = should we plot the group data (with cross-subject error bars)?
%     plot_each [0]         = should we leave each of the individual plots open?  (they are necessarily created by calling UNR_plotTC for each dataset)
%     print_each [0]        = should we automatically print each of the individual plots?  Requires plot_each==1
%     dock_figs [0]         = should we put all the individual dataset plots in the same window?
%     save_figs [0]         = should we save a pdf version of each figure in a ./figures/ subdirectory?
%
% Main output stored in G (for "group") structure

% REHMB 01.09 - created, from load group_TC
%       05.14 - modified for UNR environment


%% validate input stucture
if ~exist('G','var')
    error('G structure does not exist')
end
req_fields = {'subj_initials' 'suffix' 'rois' 'sorting' 'reg_include' 'plottype'};
if any(~ismember(req_fields,fieldnames(G)))
    missing_fields = req_fields(~ismember(req_fields,fieldnames(G)));
    m = repmat({', '},1,2*length(missing_fields)-1);
    m(1:2:length(m)) = missing_fields;
    error('the following required fields are not defined in G structure: %s',cell2mat(m))
end
if ~isfield(G,'options')
    G.options = struct();
end
    
% parse options
defaults.stattype            = 'rm_anova';
defaults.posthoc_test        = 'hsd';
defaults.errbars             = 'none';
defaults.plot_group          = 1;
defaults.plot_each           = 0;
defaults.print_each          = 0;
defaults.dock_figs           = 0;
defaults.save_figs           = 1;
G.options = propval(G.options,defaults);


%% parse additional variables from options
switch G.plottype{1}
    case 'session'
        warning('for plottype ''session'', this script will plot each subject''s data separetely.  no such thing as group analysis for this data');
        G.outfield = 'not applicable'; % overwrite 
        G.stattype = 'none'; % overwrite 
    case 'tc'
        warning('for plottype ''tc'', this script will plot each subject''s data separately.  for group analysis of timecourses, use average_glmTC');
        G.outfield = 'not applicable'; % overwrite 
        G.stattype = 'none'; % overwrite 
    case 'mean'
        G.outfield = G.plottype{1};
    case {'index' 'dprime' 'diff'}
        G.outfield = 'metric';
    case {'nvox' 'nvox_prop' 'nvox_anat'}
        G.outfield = 'n';
    otherwise
        error('outfield for group analysis not defined for plottype %s',G.plottype{1})
end


if G.options.print_each
    a = input('*** are you sure you want to print each plot? (y for ''yes''): ','s');
    if ~strcmp(a,'y')
        fprintf('\tPhew, glad I checked!\n')
        G.options.print_each = 0;
    else
        fprintf('\tAwesome, let''s do it!\n')
    end
end




%% load the data
fprintf('\n')
for i = 1:length(G.subj_initials);
    fprintf('processing %s...',G.subj_initials{i})
    this_file = fullfile('./stored/',[G.subj_initials{i} G.suffix '.mat']);
    if exist(this_file,'file')
        load(this_file);
    else
        error('can''t locate %s',this_file)
    end

    % set up a figure for storing all the resulting plots
    if G.options.dock_figs
        this_fig = figure('Name',G.subj_initials{i},'WindowStyle','docked');
    else
        this_fig = figure('Name',G.subj_initials{i});
    end
    
    % call UNR_plotTC
    o = UNR_plotTC(TC,G.rois,G.sorting,G.plottype,G.reg_include,G.options.errbars,this_fig);
    
    % print and/or close
    if ~G.options.plot_each
        close
    elseif G.options.print_each
        orient landscape
        print -dpsc
        close
        if ismember(G.plottype{1},{'session'})
            % then UNR_plotTC created two figures, print the second one too
            orient landscape; print -dpsc; close            
        end
    end
    %%G.reg_labels{i} = TC.reg_labels;
        
    
    % for 'tc' or 'session', there is no actual group analysis.  we can continue on.
    if ismember(G.plottype{1},{'tc' 'timecourse' 'session'})
        % no group statistics to run
        fprintf('complete\n',G.subj_initials{i})
        continue
    end
    
    
    % extract relevant info
    % unpack matrix output (for multiple ROIs and multple conditions) into some useful vector format
    this_m = o.(G.plottype{1}).(G.outfield); % matrix-formatted output of UNR_plotTC (nROIs x nREGS)
    this_v = this_m'; this_v = this_v(:); % vector format
    this_roi = repmat(o.rois,length(o.reg_labels),1); this_roi = this_roi(:)';% roi labels
    this_reg = repmat(o.reg_labels,length(o.rois),1)'; this_reg = this_reg(:)';% reg labels
    if i > 1
        % verify the roi and reg labels are consistent with previous dataset
        if any(~(strcmp(this_roi,G.data_roi))) || any(~(strcmp(this_reg,G.data_reg)))
            error('inconsistent return matrices across datasets.  roi and/or reg labels don''t appear to match.  check this_roi/reg and G.data_roi/reg')
        end
    else
        % we'll need these to (1) compare with other datasets and (2) for grouping variables for anovan
        G.data_roi = this_roi;
        G.data_reg = this_reg;
    end
    G.data(i,:) = this_v';

    fprintf('complete\n',G.subj_initials{i})
end

if ismember(G.plottype{1},{'tc' 'timecourse' 'session'})
    % no group statistics to run
    return
end

if all(isnan(G.data(:)))
    % then there was no valid data to extract from any ROI in any session
    warning('no valid data extracted from any dataset for ROI=%s',mat2str(cell2mat(G.rois)))
    return
end


% convert list of subjects to one big string for marking plots and filenames
subj_initials_str = ['{' G.subj_initials{1}];
for i = 2:length(G.subj_initials)
    subj_initials_str = [subj_initials_str ' ' G.subj_initials{i}];
end
subj_initials_str = [subj_initials_str '}'];


%% generate stats and associated plots
switch G.options.stattype
    case 'rm_anova'
        % repeated-measures anova treating subject as a random effect
        
        % reformat data for anovan
        nsubj = length(G.subj_initials);
        this_subj = repmat(1:nsubj,1,size(G.data,2))'; % arbitrary subject id label (numerals)
        this_data = G.data(:); % dependent variable - vector
        this_roi  = repmat(G.data_roi,nsubj,1); this_roi = this_roi(:); % roi grouping variable
        this_reg  = repmat(G.data_reg,nsubj,1); this_reg = this_reg(:); % reg (condition) grouping variable
        
        anova_factors = 'byroi'; % 'omnibus' for anova over rois and regressors; 'byroi' for separate anovas for each ROI;
        switch anova_factors
            case 'omnibus'
                % setup some variables for the anovan and multcompare
                if length(G.rois) > 1 && length(G.reg_include) > 1
                    mcomp_dims = [1 2];
                elseif length(G.rois) > 1
                    mcomp_dims = [2];
                elseif length(G.reg_include) > 1
                    mcomp_dims = [1];
                else
                    error('hmmmm, there doesn''t appear to be any factors (roi or regs) to analyze')
                end

                % THE ANOVA
                [G.anovan_p G.anovan_tab G.anovan_stats] = anovan(this_data,{this_reg this_roi this_subj},'varnames',{'reg' 'roi' 'subj'},'random',3,'display','off');
                f = figure('MenuBar','none','Color','w');
                [G.multcompare.comp G.multcompare.means G.multcompare.h G.multcompare.gnames] = multcompare(G.anovan_stats,'ctype',G.options.posthoc_test,'dimension',mcomp_dims);
                              
                % mark info on figures
                title(G.rois)
                MarkPlot(sprintf('%s   %s   %s   %s',subj_initials_str,G.suffix,G.sorting,G.plottype{1}))

                % --- embed anova table in multcompare figure for easy print/viewing
                % extract the table info we want to print
                % convert entries to strings
                table = G.anovan_tab(:,[1:3 5:8]);
                for k = 1:size(table,2) % cols
                    maxsize = -Inf;
                    for i = 1:size(table,1) % rows
                        if ~ischar(table{i,k})
                            table{i,k} = num2str(table{i,k});
                        end
                        maxsize = max(maxsize,length(table{i,k}));
                    end
                    % all entries must be strings of the same length.  first/last characters are the separator and provide a buffer for the display
                    for i = 1:size(table,1) % rows
                        table{i,k} = [' ' table{i,k} repmat(' ',1,maxsize-length(table{i,k})) ' '];
                    end
                end
                % plot it above multcompare stuff
                xs = xlim;
                ys = ylim;
                newys = ys +[0 0.30*ys(2)]; % increase y axis by 25%
                ylim(newys)
                text(mean(xs),newys(2),cell2mat(table),'FontName','Courier','HorizontalAlignment','center','VerticalAlignment','top');
                
                if G.options.save_figs
                    % as pdf and .fig (for multcomparisons)
                    fpre = sprintf('%s%s__multcompare__%s-%s',subj_initials_str,G.suffix,G.options.posthoc_test,cell2mat(G.rois));
                    orient landscape
                    print('-dpdf',sprintf('./figures/%s.pdf',fpre));
                    saveas(gcf,sprintf('./figures/%s',fpre),'fig');
                end
                
            case 'byroi'
                % TO DO: integrate 'omnibus' and 'byroi' code with some if switches to minimize code replication (for marking plots and saving figures, etc)
                
                % loop over ROIs and run a separate ANOVA for each one
                for i_roi = 1:length(G.rois)
                    r = G.rois{i_roi};
                    
                    % select out the appropriate data
                    ref = strcmp(this_roi,r);
                    roi_subj = this_subj(ref);
                    roi_data = this_data(ref);
                    if all(isnan(roi_data))
                        % nothing to analyze for this ROI
                        continue
                    end
                    roi_reg  = this_reg(ref);
                    mcomp_dims = 1; % only 1 factor (besides subject) to consider

                    % THE ANOVA
                    [G.anovan_p.(r) G.anovan_tab.(r) G.anovan_stats.(r)] = anovan(roi_data,{roi_reg roi_subj},'varnames',{'reg' 'subj'},'random',2,'display','off');
                    f = figure('MenuBar','none','Color','w');
                    [G.multcompare.comp.(r) G.multcompare.means.(r) G.multcompare.h.(r) G.multcompare.gnames.(r)] = multcompare(G.anovan_stats.(r),'ctype',G.options.posthoc_test,'dimension',mcomp_dims);
                    
                    % mark info on figures
                    title(r); % mark which ROI this is
                    MarkPlot(sprintf('%s   %s   %s   %s',subj_initials_str,G.suffix,G.sorting,G.plottype{1}))
                    
                    % --- embed anova table in multcompare figure for easy print/viewing
                    % extract the table info we want to print
                    % convert entries to strings
                    table = G.anovan_tab.(r)(:,[1:3 5:8]);
                    for k = 1:size(table,2) % cols
                        maxsize = -Inf;
                        for i = 1:size(table,1) % rows
                            if ~ischar(table{i,k})
                                table{i,k} = num2str(table{i,k});
                            end
                            maxsize = max(maxsize,length(table{i,k}));
                        end
                        % all entries must be strings of the same length.  first/last characters are the separator and provide a buffer for the display
                        for i = 1:size(table,1) % rows
                            table{i,k} = [' ' table{i,k} repmat(' ',1,maxsize-length(table{i,k})) ' '];
                        end
                    end
                    % plot it above multcompare stuff
                    xs = xlim;
                    ys = ylim;
                    newys = ys +[0 0.30*ys(2)]; % increase y axis by 25%
                    ylim(newys)
                    text(mean(xs),newys(2),cell2mat(table),'FontName','Courier','HorizontalAlignment','center','VerticalAlignment','top');
                    
                    if G.options.save_figs
                        % as pdf and .fig (for multcomparisons)
                        fpre = sprintf('%s%s__multcompare__%s-%s',subj_initials_str,G.suffix,G.options.posthoc_test,r);
                        orient landscape
                        print('-dpdf',sprintf('./figures/%s.pdf',fpre));
                        saveas(gcf,sprintf('./figures/%s',fpre),'fig');
                    end
                    
                end
                
            case 'otherwise'
                error('invalid option for anova_factors (%s)',anova_factors)
        end
        
    case 'friedman'
        % Friedman's non-parametric within-subjects ANOVA
        error('needs to be setup still.  see rm_anova for code structure')
        
        % matlab's FRIEDMAN performs a non-parametric two-way ANOVA, adjusting for row effects.
        %   thus, subject will always be the row-factor and only one of roi and reg can be the column factor.
        if length(G.rois) > 1 && length(G.reg_include) > 1
            error('cannot perform friedman test with multiple factors (since subject will always be one of the factors')
% %         elseif length(G.rois) > 1
% %             % analysis over rois
% %             friedman_group = G.data_roi;
% %         elseif length(G.reg_include) > 1
% %             friedman_group = G.data_reg;
% %         else
% %             error('hmmmm, there doesn''t appear to be any factors (roi or regs) to analyze')
        end
        
        % G.data should already be in the appropriate format since each row represents change along a single dimention (roi OR regressor)
        [G.friedman_p G.friedman_tab G.friedman_stats] = friedman(G.data,1); % rows are subj, columns are rois/conditions
        f = figure('MenuBar','none','Color','w');
        [G.multcompare.comp G.multcompare.means G.multcompare.h G.multcompare.gnames] = multcompare(G.friedman_stats,'ctype',G.options.posthoc_test);
        title(G.rois)
        
     case 'ttest'
        % two-tailed t-test (ho: mean==0) at alpha = 0.05
        [G.ttest_h G.ttest_p] = ttest(G.data,0,0.05,'both');

        
        % setup the figure
        f = figure('MenuBar','none','Color','w');

        % reshape mean and error into roi x regressor (mean/error across subjects)
        nreg = length(G.reg_include);
        nroi = length(G.rois);
        d_mean = reshape(nanmean(G.data),nreg,nroi)';
        
        ttest_err = 'ci95'; % 'sem' or 'ci95' %%% Need to integrate into arguments and integrate with group plot errbars
        switch ttest_err
            case 'sem'
                d_err  = reshape(nanste(G.data),nreg,nroi)'; % SEM
            case 'ci95'
                d_err  = reshape(ci(G.data),nreg,nroi)'; % 95% CI
            otherwise
                error('invalid ttest_err value (%s)',ttest_err)
        end
        d_pval = reshape(G.ttest_p,nreg,nroi)'; % p-value from ttest

        barerr2(G.rois,  d_mean, d_err, 'p_vals', d_pval, 'p_thr', [0.05 0.01 0.001], 'p_mark', {'*' '+' '#'});

        
        % plot labels
        ylabel(sprintf('%s +/- %s',G.plottype{1},ttest_err));
        if iscell(G.reg_include{1})
            % then we need to convert regressor pairs to single string
            this_leg = {}; % init
            for i = 1:length(G.reg_include)
                this_pair = G.reg_include{i};
                this_leg{i} = [this_pair{1} ' vs. ' this_pair{2}];
            end
        else
            this_leg = G.reg_include;
        end
        
        % plot individual subject data
        % N.B. - this is tricky, because for multiple regressors, we probably want to only connect lines
        %         within ROIs, but for one regressor, probably want to connect all datapoints for each subject
        %      - also, if multiple regs, then bars are offset from XTick, and we'd need to recalc that offset here
        plot_subj_data = 0;
        if plot_subj_data
            hold on
            x = repmat(1:size(G.data,2),size(G.data,1),1);
            %plot(x(:),G.data(:),'or');
            plot(x',G.data','o--');
            
            % udpate legend to include single subject lines
            this_leg = cat(2,this_leg,G.subj_initials);
            hold off
        end
        
        % show legend
        legend(regexprep(this_leg,'_','\\_'),'Location','SouthWest');
        
        
        % mark info on figures
        MarkPlot(sprintf('%s   %s   %s   %s',subj_initials_str,G.suffix,G.sorting,G.plottype{1}))
               
        if G.options.save_figs
            % as pdf and .fig (for multcomparisons)
            fpre = sprintf('%s%s__ttestbar__%s-%s',subj_initials_str,G.suffix,G.plottype{1},G.rois{:});
            orient landscape
            print('-dpdf',sprintf('./figures/%s.pdf',fpre));
            saveas(gcf,sprintf('./figures/%s',fpre),'fig');
        end
        
    case 'signrank'
        error('this should be used for checking if each roi/condition distribution is greater than 0 (or some other value) - needs to be updated (see ttest code)')
        
        % one-tailed sign rank, median>0? at alpha = 0.05
        for i = 1:size(G.data,2)
            [G.signrank_p(1,i) G.signrank_h(1,i)] = signrank(G.data(:,i),0,0.05);
            G.signrank_p(1,i) = G.signrank_p(1,i)/2; % div by 2 for 1-sided test
        end
        
        f = figure('MenuBar','none','Color','w');
        barerr2(G.rois, median(G.data)', [], 'p_vals', G.signrank_p', 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'},'p_height',max(G.data(:))*1.05);
        %barerr2(G.rois, median(G.data)', [], 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'},'p_height',max(G.data(:))*1.05);
        hold on
        x = repmat(1:length(G.data),size(G.data,1),1);
        %plot(x(:),G.data(:),'o--r');
        plot(x',G.data','o--');
        leg = [{'median'} G.subj_initials];
        leg = regexprep(leg,'_','\\_');
        legend(leg);
        hold off
        
    case 'anova'
        error('outdated - should be using ''rm_anova'' for within-subjects statistics - code probably needs updating')
        % group stats and plot
        
        % for multiple ROIs, anova is across ROIs and requires that there is only one condition
        % for a single ROI, anova is across conditions and requires multiple conditions
        % N.B. we've already done the error checking to verify that multple ROIs and conditions are supplied
        if length(G.rois) > 1
            anova_labels = G.rois;
        elseif length(G.reg_include) > 1
            anova_labels = G.reg_include;
        else
            error('not sure how we could have even gotten here!')
        end
        [G.anova_p G.anova_tab G.anova_stats] = anova1(G.data,anova_labels,0);
        f = figure('MenuBar','none','Color','w');
        multcompare(G.anova_stats,'ctype',G.posthoc_test);
        
    case 'kruskalwallis'
        error('outdated - should be using ''friedman'' test as a within-subjects non-parametric test - code probably needs updating')
        
        % group stats and plot
        [G.kruskalwallis_p G.kruskalwallis_tab G.kruskalwallis_stats] = kruskalwallis(G.data,G.rois,0);
        f = figure('MenuBar','none','Color','w');
        multcompare(G.kruskalwallis_stats,'ctype',G.posthoc_test);
        
    case 'none'
        % nothing to do, skip statistical tests.  useful for quickly plotting the group results
        
    otherwise
        error('stattype (%s) not defined',G.stattype)
end




%% group plot
if G.options.plot_group
    % TODO: integrate these into argument structure so they can be manipulated by group_TC_wrapper
    %       should also probably integrate this with the execution of the RM_ANOVA (or Friedman) above, but that would be tricky.  could loop over ROIs and plot/save figs separately.  this will simply alter where the ROI loop comes into play.
    G.group_plot.within_subjects_variance_correction = 'byroi'; % 'grand', 'byroi', or 'none' (see below for details)
    G.group_plot.errbars = 'ci95'; % 'ci95' or 'sem' 
    
    figure('MenuBar','none','Color','w');
    nreg = length(G.reg_include);
    nroi = length(G.rois);
    d_mean = reshape(nanmean(G.data),nreg,nroi)'; % roi x regressor (mean across subjects)


    % within-subjects confidence intervals (Loftus and Masson, 1994; Cousineau, 2005; Morey, 2008)
    switch G.group_plot.within_subjects_variance_correction
        case 'none'
            % no correction for between-subjects variance
            d_sem = reshape(nanste(G.data),nreg,nroi)'; % roi x regressor (sem across subjects)
            d_ci  = reshape(ci(G.data,95),nreg,nroi)'; % roi x regressor (onfidence interval across subjects) (N.B. invalid when including NaN values - overestimates N for each comparison with a NaN)
    
        case 'grand'
            % assumes an omnibus ANOVA including ROI as a factor (ie, over ROI and Condition)
            
            % grand mean calculated over all ROIs and conditions
            norm_data = G.data - repmat(nanmean(G.data,2),1,size(G.data,2)) + nanmean(G.data(:));
            
            d_sem = reshape(nanste(norm_data),nreg,nroi)'; % roi x regressor (sem across subjects)
            d_ci  = reshape(ci(norm_data,95),nreg,nroi)'; % roi x regressor (onfidence interval across subjects) (N.B. invalid when including NaN values - overestimates N for each comparison with a NaN)
            
            % Morey, 200 notes that Cousineau's method is biased due to induced correlation in the data from the above step
            % total number of within-subjects conditions (product of number of levels across all factors)
            M = size(G.data,2);
            d_sem = d_sem * (M/(M-1)); % Morey, 2008 correction
            d_ci  = d_ci * (M/(M-1));  % Morey, 2008 correction

        case 'byroi'
            % assumes separate ANOVAs for each ROI (so only Condition as a factor)
            % in this case, data is normalized for each ROI independently (ie, not considering data from other ROIs when
            % calculating within-subject variance or adding back grand mean)
            
            % normalize data for each ROI separately
            norm_data = repmat(NaN,size(G.data));
            for i = 1:nroi % loop over ROIs
                r = G.rois{i};
                ref = strcmp(G.data_roi,r); % ref list for G.data columns corresponding to this ROI
                this_data = G.data(:,ref); % extract this ROIs data
                norm_data(:,ref) = this_data - repmat(nanmean(this_data,2),1,size(this_data,2)) + nanmean(this_data(:)); % normalize ROI data
            end
                
            d_sem = reshape(nanste(norm_data),nreg,nroi)'; % roi x regressor (sem across subjects)
            d_ci  = reshape(ci(norm_data,95),nreg,nroi)'; % roi x regressor (confidence interval across subjects)  (N.B. invalid when including NaN values - overestimates N for each comparison with a NaN)
            
            % Morey, 2008 notes that Cousineau's method is biased due to induced correlation in the data from the above step
            % total number of within-subjects conditions (nregs - assumes ROI-specific ANOVAs)
            M = nreg;
            d_sem = d_sem * (M/(M-1)); % Morey, 2008 correction
            d_ci  = d_ci * (M/(M-1));  % Morey, 2008 correction
            
        otherwise
            error('invalid option for G.group_plot.within_subjects_variance_correction')
    end
    
    switch G.group_plot.errbars
        case 'sem'
            this_err = d_sem;
        case 'ci95'
            this_err = d_ci;
        otherwise
            error('invalid option for G.group_plot.errbars')
    end
        
    barerr2(G.rois', d_mean, this_err, 'p_vals', [], 'p_thr', [0.05 0.01 0.001], 'p_mark', {'#'  '*' '+'});
    ylabel(sprintf('%s +/- %s (within-subject correction = %s)',G.plottype{1},G.group_plot.errbars,G.group_plot.within_subjects_variance_correction))
    if iscell(G.reg_include{1})
        % then we need to convert regressor pairs to single string
        this_leg = {}; % init
        for i = 1:length(G.reg_include)
            this_pair = G.reg_include{i};
            this_leg{i} = [this_pair{1} ' vs. ' this_pair{2}];
        end
    else
        this_leg = G.reg_include;
    end
    legend(this_leg);
    
    MarkPlot(sprintf('%s   %s   %s   %s',subj_initials_str,G.suffix,G.sorting,G.plottype{1}))
    
    if G.options.save_figs
        % as pdf
        fpre = sprintf('%s%s__bar__%s-%s',subj_initials_str,G.suffix,G.plottype{1},cell2mat(G.rois));
        orient landscape
        print('-dpdf',sprintf('./figures/%s.pdf',fpre));
    end
end
