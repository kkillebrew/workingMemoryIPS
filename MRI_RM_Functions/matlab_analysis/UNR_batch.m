% quick batch script for running scripts across multiple subjects


% adjustable params
subjs = {'CB' 'GC' 'GG' 'JV' 'KK' 'KM' 'LS' 'NS'};
hemis = {'left' 'right'};%{'left' 'right'}
rois = {'V1' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' 'IPS5' 'hFEF' 'pIPS' 'mIPS'};

recalc_indiv_TC = 1;
plot_individual = 0;
recalc_group_TC = 1;
plot_group = 0;


% relatively stable params
vox_select = {'tstat' 1.96};%{'tstat' 1.00};
exp_reg = 'both_vs_pass'; % GLM regressor (or contrast) for voxel selection

analysis_type = 'blockwise'; % 'blockwise' or 'trialwise'

datafile_prefix = 'tscvrdtm_norm';
norm_type = 'pre'; % 'each' or 'pre'

bucket_prefix = 'tscvrsm4_norm_bucket_REML'; %'tscvrsm4_norm_TRIALS_bucket_REML'; 

switch analysis_type
    case 'blockwise'
        baseline_idx = [];%-5:-2;
        nPostTRs = 26;
        nPreTRs = 10;
        dprime_TRs = 2:19;
    case 'trialwise'
        baseline_idx = 1;%-4:-2;
        nPostTRs = 10; %26
        nPreTRs = 10;
        dprime_TRs = 2:3;
    otherwise
        error('invalid analysis_type (%s)',analysis_type)
end

% inidividual subject plots
for hscell = hemis
    hs = hscell{1}; % convert to string
    
    if recalc_indiv_TC
        for subjcell = subjs
            subj = subjcell{1}; % convert to string
            fprintf('\n== subj %s\n',subj);
            
            TC = UNR_rawTC(...
                subj,...
                1:7,...
                rois,...
                hs,...
                vox_select,...
                'exp_reg',exp_reg,...
                'datafile_prefix',datafile_prefix,...
                'bucket_prefix',bucket_prefix,...
                'postTRs',nPostTRs,...
                'preTRs',nPreTRs,...
                'norm_type',norm_type,...
                'baseline_idx',baseline_idx,...
                'output_tail',[analysis_type '_' exp_reg]);
            
            if plot_individual
                UNR_plotTC(TC,[],'byreg','tc',[],'none');
            end
        end
    end
end

% group averaged plots
for hscell = hemis
    hs = hscell{1}; % convert to string
    
    if recalc_group_TC
        TC = UNR_average_rawTC(...
            subjs,...
            hs,...
            vox_select,...
            'exp_reg',exp_reg,...
            'datafile_prefix',datafile_prefix,...
            'bucket_prefix',bucket_prefix,...
            'postTRs',nPostTRs,...
            'preTRs',nPreTRs,...
            'norm_type',norm_type,...
            'baseline_idx',baseline_idx,...
            'output_tail',[analysis_type '_' exp_reg]);
        
        if plot_group
            UNR_plotTC(TC,[],'byreg','tc',[],'none');
            UNR_plotTC(TC,[],'byreg',{'dprime' dprime_TRs},{{'let_left' 'let_right'} {'ori_left' 'ori_right'}},'sem')
        end
    end
end