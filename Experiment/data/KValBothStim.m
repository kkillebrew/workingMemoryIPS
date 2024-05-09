clear all;
close all;

% file_list={'KWK_orient_pilot_wm_020514_001', 'CDB_orient_pilot_wm_020514_001','GG_orient_pilot_wm_020614_001','DM_orient_pilot_wm_021114_001','JEV_orient_pilot_wm_021914_001','RS_orient_pilot_wm_021914_001'};
% file_list={'KWK_letters_pilot_wm_021114_001','CDB_letters_pilot_wm_021114_001','GG_letters_pilot_wm_021114_001', 'DM_letters_pilot_wm_021114_001'};
% file_list_letters={'KWK_lettersRK_pilot_wm_021214_001','CDB_lettersKR_pilot_wm_021214_001','GG_lettersKR_pilot_wm_021214_001','DM_lettersRK_pilot_wm_021214_001','JEV_letters_pilot_wm_021914_001','RS_letters_pilot_wm_021914_001'};
% file_list = {'KK_wm_030414_001','KK_wm_030414_002','KK_wm_030414_003','KK_wm_030414_004','KK_wm_030414_005','KK_wm_030414_006','KK_wm_030414_007'};
file_list = {'GG_wm_030414_001','GG_wm_030414_002','GG_wm_030414_003','GG_wm_030414_004','GG_wm_030414_005','GG_wm_030414_006','GG_wm_030414_007'};

% calculate k val K = S*(HR-FAR)
% S=Set Size HR = Hit Rate FAR = False Alarm Rate
% HR consists of hits and correct rejects and FAR consists of FA and miss
% Depends on how many actually changed and how many didn't change (aka
% which trials were new and which were old)

% Kyle, to calculate K for a given condition, you need to know the set size (S),
% the hit rate (HR, proportion of hits), the false alarm rate (FAR, proportion of false alarms)
% and then use the formula:   K =  S * (HR - FAR).  Your hit rate depends on how many trials there
% were per condition and how many of those trials were OLD trials.  So if you had 50 trials per
% condition and half (25) were OLD trials and let's say, for example, that there were 20 "hit"
% trials, then you would have a hit rate of  HR = 20/25 = 0.80.    The idea holds for calculating
% the false alarm rate, but using the other half of the trials (NEW trials). Feel free to ask
% me more about this if you get stuck.

% Make sure to take out any trials that are NaN (no response) and to also
% calculate the precentage of response that were misses

%  -  passive vs. memory is in trials.task
%  -  trials per block -- trials.block
%  -  trial type (aka left right both) trials.mem_side
%  -  new/old trials -- trials.targ_change
%  -  is response what it should be - data.status

accuracyOrientation=[];
accuracyPassiveOrientation=[];
KLeftOrientation=[];
KRightOrientation=[];
KBothOrientation=[];

accuracyLetters=[];
accuracyPassiveLetters=[];
KLeftLetters=[];
KRightLetters=[];
KBothLetters=[];

xLabelsOrientation={'Left','Right','Either'};
xLabelsLetters={'Orientation', 'Letters'};

%% Calculate the accuracy and K scores for each participant and condition for orientation
for p=1:length(file_list)
    load(file_list{p});
    
    count1=1;
    count2=1;
    for i=1:trials.n
        if strcmp(trials.stimtype{i},'orientation')
            data.orientation.status(count1) = data.status(i);
            count1=count1+1;
        elseif strcmp(trials.stimtype{i},'letters')
            data.letters.status(count2) = data.status(i);
            count2=count2+1;
        end
    end
    
    perf.orientation.overall = 100 * nanmean(data.status(strcmp(trials.stimtype,'orientation')));
    perf.letters.overall = 100 * nanmean(data.status(strcmp(trials.stimtype,'letters')));
    
    %% Determine accuracy without passive viewing
    accuracyCountOrientation=1;
    accuracyCountLetters=1;
    countO=1;
    countL=1;
    if strcmp(trials.stimtype{i},'orientation')
        for i=1:trials.n
            if trials.task(i)==2
                if data.status(i)==1
                    accuracyCountOrientation=accuracyCountOrientation+1;
                end
                countO=countO+1;
            end
        end
    elseif strcmp(trials.stimtype{i},'letters')
        for i=1:trials.n
            if trials.task(i)==2
                if data.status(i)==1
                    accuracyCountLetters=accuracyCountLetters+1;
                end
                countL=countL+1;
            end
        end
    end
    
    accuracyLetters(p)=100*(accuracyCountLetters/countL);
    accuracyPassiveLetters(p)=perf.letters.overall;
    accuracyOrientation(p)=100*(accuracyCountOrientation/countO);
    accuracyPassiveOrientation(p)=perf.orientation.overall; 
    
    %% Determine the K score for each condition (left, right, both)
    newTotalLeftOrientation = 0;
    newCorrectLeftOrientation = 0;
    oldCorrectLeftOrientation = 0;
    newTotalRightOrientation = 0;
    newCorrectRightOrientation = 0;
    oldCorrectRightOrientation = 0;
    newTotalBothOrientation = 0;
    newCorrectBothOrientation = 0;
    oldCorrectBothOrientation = 0;
    oldTotalLeftOrientation = 0;
    oldTotalRightOrientation = 0;
    oldTotalBothOrientation = 0;
    
    newTotalLeftLetters = 0;
    newCorrectLeftLetters = 0;
    oldCorrectLeftLetters = 0;
    newTotalRightLetters = 0;
    newCorrectRightLetters = 0;
    oldCorrectRightLetters = 0;
    newTotalBothLetters = 0;
    newCorrectBothLetters = 0;
    oldCorrectBothLetters = 0;
    oldTotalLeftLetters = 0;
    oldTotalRightLetters = 0;
    oldTotalBothLetters = 0;
    
    % Determine for both old and new trial types accuracy rate
    for i=1:trials.n
        if strcmp(trials.stimtype{i},'orientation')
            % Make sure passive viewing is not included
            if trials.task(i)==2
                % Determines which condition you are in
                if trials.mem_side(i)==-1
                    % Was it a new or an old trial
                    if trials.targ_change(i)==1
                        newTotalLeftOrientation = newTotalLeftOrientation+1;
                        % Did you respond correctly
                        if data.status(i)==1
                            newCorrectLeftOrientation = newCorrectLeftOrientation+1;
                        end
                    else
                        oldTotalLeftOrientation = oldTotalLeftOrientation+1;
                        if data.status(i)==0
                            oldCorrectLeftOrientation = oldCorrectLeftOrientation+1;
                        end
                    end
                elseif trials.mem_side(i)==1
                    if trials.targ_change(i)==1
                        newTotalRightOrientation = newTotalRightOrientation+1;
                        if data.status(i)==1
                            newCorrectRightOrientation = newCorrectRightOrientation+1;
                        end
                    else
                        oldTotalRightOrientation = oldTotalRightOrientation+1;
                        if data.status(i)==0
                            oldCorrectRightOrientation = oldCorrectRightOrientation+1;
                        end
                    end
                elseif trials.mem_side(i)==2
                    if trials.targ_change(i)==1
                        newTotalBothOrientation = newTotalBothOrientation+1;
                        if data.status(i)==1
                            newCorrectBothOrientation = newCorrectBothOrientation+1;
                        end
                    else
                        oldTotalBothOrientation = oldTotalBothOrientation+1;
                        if data.status(i)==0
                            oldCorrectBothOrientation = oldCorrectBothOrientation+1;
                        end
                    end
                end
            end
        end
    end
    
    % Determine for both old and new trial types accuracy rate
    for i=1:trials.n
        if strcmp(trials.stimtype{i},'letters')
            % Make sure passive viewing is not included
            if trials.task(i)==2
                % Determines which condition you are in
                if trials.mem_side(i)==-1
                    % Was it a new or an old trial
                    if trials.targ_change(i)==1
                        newTotalLeftLetters = newTotalLeftLetters+1;
                        % Did you respond correctly
                        if data.status(i)==1 
                            newCorrectLeftLetters = newCorrectLeftLetters+1;
                        end
                    else
                        oldTotalLeftLetters = oldTotalLeftLetters+1;
                        if data.status(i)==0
                            oldCorrectLeftLetters = oldCorrectLeftLetters+1;
                        end
                    end
                elseif trials.mem_side(i)==1
                    if trials.targ_change(i)==1
                        newTotalRightLetters = newTotalRightLetters+1;
                        if data.status(i)==1
                            newCorrectRightLetters = newCorrectRightLetters+1;
                        end
                    else
                        oldTotalRightLetters = oldTotalRightLetters+1;
                        if data.status(i)==0
                            oldCorrectRightLetters = oldCorrectRightLetters+1;
                        end
                    end
                elseif trials.mem_side(i)==2
                    if trials.targ_change(i)==1
                        newTotalBothLetters = newTotalBothLetters+1;
                        if data.status(i)==1
                            newCorrectBothLetters = newCorrectBothLetters+1;
                        end
                    else
                        oldTotalBothLetters = oldTotalBothLetters+1;
                        if data.status(i)==0
                            oldCorrectBothLetters = oldCorrectBothLetters+1;
                        end
                    end
                end
            end
        end
    end
    
%     if newTotalLeftLetters == 0
%         HRLeftLetters = 1;
%     else
        HRLeftLetters = newCorrectLeftLetters/newTotalLeftLetters;
%     end
    FARLeftLetters = oldCorrectLeftLetters/oldTotalLeftLetters;
    
    HRRightLetters = newCorrectRightLetters/newTotalRightLetters;
    FARRightLetters = oldCorrectRightLetters/oldTotalRightLetters;
    
    HRBothLetters = newCorrectBothLetters/newTotalBothLetters;
    FARBothLetters = oldCorrectBothLetters/oldTotalBothLetters;
    
    KLeftLetters(p)=params.blocks.set_size*(HRLeftLetters-FARLeftLetters);
    KRightLetters(p)=params.blocks.set_size*(HRRightLetters-FARRightLetters);
    KBothLetters(p)=params.blocks.set_size*(HRBothLetters-FARBothLetters);
    
    
    
    
%     if newTotalLeftOrientation == 0
%         HRLeftOrientation = 1;
%     else
        HRLeftOrientation = newCorrectLeftOrientation/newTotalLeftOrientation;
%     end
    FARLeftOrientation = oldCorrectLeftOrientation/oldTotalLeftOrientation;
    
    HRRightOrientation = newCorrectRightOrientation/newTotalRightOrientation;
    FARRightOrientation = oldCorrectRightOrientation/oldTotalRightOrientation;
    
    HRBothOrientation = newCorrectBothOrientation/newTotalBothOrientation;
    FARBothOrientation = oldCorrectBothOrientation/oldTotalBothOrientation;
    
    KLeftOrientation(p)=params.blocks.set_size*(HRLeftOrientation-FARLeftOrientation);
    KRightOrientation(p)=params.blocks.set_size*(HRRightOrientation-FARRightOrientation);
    KBothOrientation(p)=params.blocks.set_size*(HRBothOrientation-FARBothOrientation);
    
    
    
    
    disp(sprintf('\n%s%.2f%s%s','Accuracy without passive trials for orientation: ', accuracyOrientation(p), '% for ', file_list{p}));
    disp(sprintf('%s%.2f%s%s\n','Accuracy with passive trials (inflated) for orientation: ', accuracyPassiveOrientation(p), '% for participant ', file_list{p}));
    
    disp(sprintf('%s%.1f%s','K Score for left trials for orientation: ', KLeftOrientation(p), file_list{p}));
    disp(sprintf('%s%.1f%s','K Score for right trials for orientation: ', KRightOrientation(p), file_list{p}));
    disp(sprintf('%s%.1f%s%s\n','K Score for both trials for orientation: ', KBothOrientation(p), file_list{p}));
    
    disp(sprintf('\n%s','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    
    disp(sprintf('\n%s%.2f%s%s','Accuracy without passive trials for letters: ', accuracyLetters(p), '% for ', file_list{p}));
    disp(sprintf('%s%.2f%s%s\n','Accuracy with passive trials (inflated) trials for letters: ', accuracyPassiveLetters(p), '% for ', file_list{p}));
    
    disp(sprintf('%s%.1f%s','K Score for left trials for letters: ', KLeftLetters(p), file_list{p}));
    disp(sprintf('%s%.1f%s','K Score for right trials for letters: ', KRightLetters(p), file_list{p}));
    disp(sprintf('%s%.1f%s%s\n','K Score for both trials for letters: ', KBothLetters(p), file_list{p}));
    
    disp(sprintf('\n%s','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
end

disp(sprintf('\n%s','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
%% Calculate, plot, and display averages
meanLeftOrientation = mean(KLeftOrientation);
meanRightOrientation = mean(KRightOrientation);
meanBothOrientation = mean(KBothOrientation);
meanLeftLetters = mean(KLeftLetters);
meanRightLetters = mean(KRightLetters);
meanBothLetters = mean(KBothLetters);

meanAccuracyOrientation = mean(accuracyOrientation);
meanAccuracyPassiveOrientation = mean(accuracyPassiveOrientation);
meanAccuracyLetters = mean(accuracyLetters);
meanAccuracyPassiveLetters = mean(accuracyPassiveLetters);

disp(sprintf('\n%s%.2f%s','Average accuracy without passive trials for orientation: ', meanAccuracyOrientation, '%'));
disp(sprintf('%s%.2f%s\n','Average accuracy with passive trials (inflated) for orientation: ', meanAccuracyPassiveOrientation,'%'));
disp(sprintf('\n%s%.2f%s','Average accuracy without passive trials for Letters: ', meanAccuracyLetters, '%'));
disp(sprintf('%s%.2f%s\n','Average accuracy with passive trials (inflated) for Letters: ', meanAccuracyPassiveLetters,'%'));

disp(sprintf('%s%.1f','K Score for left trials for orientation: ', meanLeftOrientation));
disp(sprintf('%s%.1f','K Score for right trials for orientation: ', meanRightOrientation));
disp(sprintf('%s%.1f%s\n','K Score for both trials for orientation: ', meanBothOrientation));
disp(sprintf('%s%.1f','K Score for left trials for Letters: ', meanLeftLetters));
disp(sprintf('%s%.1f','K Score for right trials for Letters: ', meanRightLetters));
disp(sprintf('%s%.1f%s\n','K Score for both trials for Letters: ', meanBothLetters));

figure()
subplot(1,2,1)
c = [meanLeftOrientation meanLeftLetters; meanRightOrientation meanRightLetters; meanBothOrientation  meanBothLetters];
y=c(1:3,:);
bar(y);
str = {'','Average K Score Across Participants',''}; % cell-array method
xlabel('Hemifield','FontSize',15);
ylabel('K Score','FontSize',15);
set(gca, 'XTickLabel',xLabelsOrientation, 'XTick',1:numel(xLabelsOrientation))
title(str,'FontSize',15,'FontWeight','bold');
legend({'Orientation','Letters'});

subplot(1,2,2)
c = [meanAccuracyOrientation meanAccuracyLetters];
bar(c);
str = {'','Average Accuracy Across Participants',''}; % cell-array method
xlabel('Condition','FontSize',15);
ylabel('Accuracy (%)','FontSize',15);
set(gca, 'XTickLabel',xLabelsLetters, 'XTick',1:numel(xLabelsLetters));
title(str,'FontSize',15,'FontWeight','bold');
set(gca,'ylim',[0,100]);



