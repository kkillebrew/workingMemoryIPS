% This is a wrapper script for the centralized group_TC_script.  See that script for help.
% This is a special verion o fhte wrapper wscript that will perform the
% same analysis on each hemipshere, combine the data into a single G
% structure, and perform a 3-way rm_anova (the code for the ANOVA was
% modified from group_TC_script.

clear G options

%---------------------------------------------------------------------------------------------------------
% required fields - defined in G structure
%       rois        = cell array of rois to include, or 'all' for all rois that are part of the loaded files.
%       sorting     = 'byreg' or 'grand' (passed into plotTC).
%       reg_include = cell array of regressors to include in the analysis (passed into plotTC).
%       plottype    = passed into plotTC, see plotTC for options.  For 'tc' and 'session' there
%                     is no actual group analysis to perform.  This script is just a shortcut for
%                     batch plotting and printing.
G.subj_initials = {'CB' 'GC' 'GG' 'JV' 'KK' 'KM' 'LS' 'NS'}; % currently those with rough IPS rois


% trial-wise or block-wise analysis?
analysis_type = 'blockwise'; % 'blockwise' or 'trialwise'
switch analysis_type
    case 'blockwise'
        this_suffix = '_prenorm_tstat1_blockwise_both_vs_pass_TC'; % no hemisphere
        %this_suffix = '_prenorm_blockwise_TC'; % no hemisphere
        
        peak_range = 2:19; % entire block
        %peak_range = 2:10; % first half of block
        %peak_range = 11:19; % second half of block
        %peak_range = 3:3:18; % just the trial peaks
    case 'trialwise'
        this_suffix = '_eachnorm_tstat1_trialwise_TC'; % no hemisphere
        peak_range = 2:3; % TimePoint index of peak range for mean and index/dprime
    otherwise
        error('invalid analysis_type (%s)',analysis_type)
end

% ROIs to include
G.rois = {'pIPS' 'mIPS'} ;% IPS0-2 vs. IPS3/4
%G.rois = {'pIPS' 'mIPS' 'IPS5'} ;% IPS0-2 vs. IPS3/4 vs. IPS5
G.rois = {'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' 'IPS5'} ;% all IPS ROIs


G.sorting = 'byreg';


which_analysis = 'nvox_prop';
switch which_analysis
    case {'nvox' 'nvox_prop' 'nvox_anat'}
        G.plottype = {which_analysis}; % see plotTC
        G.reg_include = {'let_left'};
        invert_lh_data = 0; % boolean, should we invert the LH data before the ANOVA?
        fix_yaxis_for_100= 0; % boolean, should we adjust visible y-axis for a plot centered on 100?
        
    case 'mean'
        G.plottype = {'mean' peak_range}; % for TC time course
        G.reg_include = {'let_left' 'let_right' 'ori_left' 'ori_right'};
        %G.reg_include = {'let_left' 'let_right'};
        %G.reg_include = {'ori_left' 'ori_right'};
        invert_lh_data = 0; % boolean, should we invert the LH data before the ANOVA?
        fix_yaxis_for_100 = 1; % boolean, should we subtract 100 from all data points before the ANOVA (and plots)?

    case {'dprime' 'index' 'diff'}
        G.plottype = {which_analysis peak_range}; % for TC time course
        G.reg_include = {{'let_left' 'let_right'} {'ori_left' 'ori_right'}};
        invert_lh_data = 1; % boolean, should we invert the LH data before the ANOVA?
        fix_yaxis_for_100 = 0; % boolean, should we subtract 100 from all data points before the ANOVA (and plots)?
        
    case 'tc'
        error('for group timecoursea analysis, just load the group-averaged dataset and plot directly')
        
    otherwise
        error('which_analysis %s not defined',which_analysis)
end



%---------------------------------------------------------------------------------------------------------
% optional fields - defined in G.options structure
G.options.stattype = 'none'; % 'rm_anova' or 'friedman' or 'ttest'
G.options.posthoc_test = 'hsd';%'dunn-sidak'; % see multcompare's 'ctype' option

G.options.errbars = 'none'; % passed into plotTC

G.options.plot_group = 1;
G.options.plot_each  = 0;
G.options.print_each = 0;
G.options.dock_figs  = 0;
G.options.save_figs  = 1;


%---------------------------------------------------------------------------------------------------------
% call the centralized script, once for each hemisphere

% RH
hemisphere = 'rh';
G.suffix = ['_' hemisphere this_suffix];
UNR_group_TC_script
if fix_yaxis_for_100
    ylim([100 100.5]);
end
Grh = G;

% LH
hemisphere = 'lh';
G.suffix = ['_' hemisphere this_suffix];
UNR_group_TC_script
if fix_yaxis_for_100
    ylim([100 100.5]);
end

% COMBINE
n_per_hemi = length(G.rois)*length(G.reg_include);
G.data_hs = repmat({'lh' 'rh'},n_per_hemi,1); G.data_hs = G.data_hs(:)';
G.data_roi = repmat(G.data_roi,1,2);
G.data_reg = repmat(G.data_reg,1,2);

if invert_lh_data
    % invert the data from the left hemisphere so all (expected) indices are positive
    % this is useful for converting left-right into contra-ipsilateral comparisons
    G.data = [-G.data Grh.data]; 
else
    G.data = [G.data Grh.data]; 
end


G.suffix = regexprep(G.suffix,'_lh_','_hs_');
G.hemis = {'lh' 'rh'};

% UPDATE
G.options.stattype = 'rm_anova'; % 'rm_anova' or 'friedman' or 'ttest'
G.options.posthoc_test = 'hsd';%'dunn-sidak'; % see multcompare's 'ctype' option





% generate stats and associated plots
switch G.options.stattype
    case 'rm_anova'
        % repeated-measures anova treating subject as a random effect
        
        % reformat data for anovan
        nsubj = size(G.data,1);%length(G.subj_initials);
        this_subj = repmat(1:nsubj,1,size(G.data,2))'; % arbitrary subject id label (numerals)
        this_data = G.data(:); % dependent variable - vector
        this_hs   = repmat(G.data_hs,nsubj,1); this_hs = this_hs(:); % hs grouping variable
        this_roi  = repmat(G.data_roi,nsubj,1); this_roi = this_roi(:); % roi grouping variable
        this_reg  = repmat(G.data_reg,nsubj,1); this_reg = this_reg(:); % reg (condition) grouping variable
        
        anova_factors = 'omnibus'; % 'omnibus' for anova over rois and regressors; 'byroi' for separate anovas for each ROI;
        switch anova_factors
            case 'omnibus'
                % setup some variables for the anovan and multcompare
                mcomp_dims = find([length(G.reg_include) length(G.rois) length(G.hemis)] > 1);
                if ~sum(mcomp_dims) 
                    error('hmmmm, there doesn''t appear to be any factors (roi or regs) to analyze')
                end

                % THE ANOVA
                [G.anovan_p G.anovan_tab G.anovan_stats] = anovan(this_data,{this_reg this_roi this_hs this_subj},'varnames',{'reg' 'roi' 'hs' 'subj'},'random',4,'display','off','model','full');
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
                    roi_hs = this_hs(ref);
                    
                    mcomp_dims = find([length(G.reg_include) length(G.hemis)] > 1);
                    if ~sum(mcomp_dims) 
                        error('hmmmm, there doesn''t appear to be any factors (roi or regs) to analyze')
                    end
                    
                    % THE ANOVA
                    [G.anovan_p.(r) G.anovan_tab.(r) G.anovan_stats.(r)] = anovan(roi_data,{roi_reg roi_hs roi_subj},'varnames',{'reg' 'hs' 'subj'},'random',3,'display','off','model','full');
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
        
    case 'none'
        % nothing to do, skip statistical tests.  useful for quickly plotting the group results
        
    otherwise
        error('stattype (%s) not defined',G.stattype)
end


% print anovan results to screen
G.anovan_tab(:,[1:7])
