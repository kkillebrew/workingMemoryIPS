% quick batch script for running scripts across multiple subjects


% adjustable params
subjs = {'CB' 'GC' 'GG' 'JV' 'KK' 'KM' 'LS' 'NS'};
hemis = {'lh' 'rh'};%{'left' 'right'}

analysis_type = 'blockwise'; % 'blockwise' or 'trialwise'
switch analysis_type
    case 'blockwise'
        % no hemisphere in file name here...
        this_suffix = '_prenorm_tstat1_blockwise_both_TC'; 
        %this_suffix = '_prenorm_blockwise_TC'; % no hemisphere

    case 'trialwise'
        this_suffix = '_eachnorm_tstat1_trialwise_TC'; % no hemisphere
    otherwise
        error('invalid analysis_type (%s)',analysis_type)
end


rois = {'V1' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' 'IPS5' 'hFEF' 'pIPS' 'mIPS'};
rois = {'pIPS' 'mIPS'};

sorting = 'byreg';  % 'byreg' or 'grand'
%sorting = 'grand';

reg_include = {'let_both' 'let_pass' 'ori_both' 'ori_pass'};
%reg_include = {'let_left' 'let_right' 'ori_left' 'ori_right'};
%reg_include = {'ori_left' 'ori_right'};
%reg_include = {'let_left' 'let_right'};

plot_individual = 0;
plot_group = 1;



% inidividual subject plots
for subjcell = subjs
    this_subj = subjcell{1};
    for hscell = hemis
        hs = hscell{1}; % convert to string
        if plot_individual
            fprintf('%s %s\n',this_subj,hs)

            this_dfile = [this_subj '_' hs this_suffix];
            load(['./stored/' this_dfile]);
            UNR_plotTC(TC,rois,sorting,'tc',reg_include,'none');
            
            subplot(2,1,1);vline(0:6:36);subplot(2,1,2);vline(0:6:36);
        end
    end
end


% group averaged plots
for hscell = hemis
    hs = hscell{1}; % convert to string
    if plot_group
        fprintf('%s %s\n','group',hs)
        this_dfile = [cell2mat(subjs) '_' hs this_suffix];
        load(['./stored/' this_dfile]);
        UNR_plotTC(TC,rois,sorting,'tc',reg_include,'none');
        
        subplot(2,1,1);vline(0:6:36);subplot(2,1,2);vline(0:6:36);
        
    end
end