% This is a wrapper script for the centralized group_TC_script.  See that script for help.

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
analysis_type = 'trialwise'; % 'blockwise' or 'trialwise'
switch analysis_type
    case 'blockwise'
        this_suffix = '_prenorm_tstat1_blockwise_TC'; % no hemisphere
        peak_range = 1:19; % TimePoint index of peak range for mean and index/dprime
    case 'trialwise'
        this_suffix = '_eachnorm_tstat1_trialwise_TC'; % no hemisphere
        peak_range = 2:3; % TimePoint index of peak range for mean and index/dprime
    otherwise
        error('invalid analysis_type (%s)',analysis_type)
end

hemisphere = 'rh'; % 'lh' or 'rh'
G.suffix = ['_' hemisphere this_suffix];



% ROIs organized by Figure

% Figure X - description
%G.rois = {'V1' 'IPS0' 'IPS1' 'IPS2'};

% small test set
G.rois = {'V1' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' 'IPS5' 'hFEF' 'pIPS' 'mIPS'};




G.sorting = 'byreg';


which_analysis = 'dprime';
switch which_analysis
    case 'nvox'
        G.plottype = {'nvox'}; % see plotTC
        G.reg_include = {'let_left'};
        
    case 'mean'
        G.plottype = {'mean' peak_range}; % for TC time course
        G.reg_include = {'let_left' 'let_right' 'ori_left' 'ori_right'};

    case {'dprime' 'index' 'diff'}
        G.plottype = {which_analysis peak_range}; % for TC time course
        G.reg_include = {{'let_left' 'let_right'} {'ori_left' 'ori_right'}};
        
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
% call the centralized script
% for i_roi = 1:length(toloop.rois)
%     G.rois = toloop.rois{i_roi};
% 
%     % now that we've updated some arguments, do
UNR_group_TC_script
% end
