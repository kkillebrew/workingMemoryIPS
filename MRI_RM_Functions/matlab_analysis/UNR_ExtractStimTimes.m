function [stimtimes] = UNR_ExtractStimTimes(subj_initials,analysis_type)
% function to extract stimtimes from textfiles that were used as input to
% 3dDevconvolve


% the following assumes all stim_times files are located in
% $SUBJ/scripts/stim_times/ and are called $SUBJ_$REGNAME.1D
%regnames = {'let_left' 'let_right' 'ori_left' 'ori_right'};
%regnames = {'let_both' 'let_pass' 'ori_both' 'ori_pass'};
regnames = {'let_left' 'let_right' 'ori_left' 'ori_right' 'let_both' 'let_pass' 'ori_both' 'ori_pass'};
for rncell = regnames
    rn = rncell{1}; % convert to string
    this_file = ['../' subj_initials '/scripts/stim_times/' subj_initials '_' rn '.1D'];
    first_block_trial = num2cell(textread(this_file,'%n')); % textread will not work in future, and currently this only works because there is a single block/run.  but it works for the current expeirment    

    tmp = regexp(analysis_type,'^.*wise','match');
    switch tmp{1}
        case 'blockwise'
            stimtimes.(rn) = first_block_trial;
        case 'trialwise'
            stimtimes.(rn) = cellfun(@(x) x:6:x+30,first_block_trial,'UniformOutput',false);
        otherwise
            error('undefined analysis_type (%s)',analysis_type)
    end
end


