clear all 
close all

% RM Email
% Looking at the neural data there appears to be an indication that early
% trials in a block show larger effects (for the spatial bias metric, at
% least) than later trials in a block.  I am wondering whether there is a
% corresponding change in behavioral performance across the 6 trials in a
% block, which could indicate some fatigue or something going on.  Can you
% create plots of the behavioral performance by condition, as a function of
% trial order within a block?  I'm thinking of a line plot with a separate
% line for each condition (left-letters, right-letters, left-ori, right-ori)
% Cowan's K on the y-axis, with the x-axis being Trial Order in a Block (1-6).


clear all;
close all;

file_list{1} = {'AW_wm_072214_001','AW_wm_072214_002','AW_wm_072214_003','AW_wm_072214_004','AW_wm_072214_005','AW_wm_072214_006','AW_wm_072214_007'};
file_list{2} = {'CDB_wm_031414_001','CDB_wm_031414_002','CDB_wm_031414_003','CDB_wm_031414_004','CDB_wm_031414_005','CDB_wm_031414_006','CDB_wm_031414_007'};
file_list{3} = {'GG_wm_030414_001','GG_wm_030414_002','GG_wm_030414_003','GG_wm_030414_004','GG_wm_030414_005','GG_wm_030414_006','GG_wm_030414_007'};
file_list{4} = {'GPC_wm_031814_001','GPC_wm_031814_002','GPC_wm_031814_003','GPC_wm_031814_004','GPC_wm_031814_005','GPC_wm_031814_006','GPC_wm_031814_007',};
file_list{5} = {'JEV_wm_031414_001','JEV_wm_031414_002','JEV_wm_031414_003','JEV_wm_031414_004','JEV_wm_031414_005','JEV_wm_031414_006','JEV_wm_031414_007',};
file_list{6} = {'KK_wm_030414_001','KK_wm_030414_002','KK_wm_030414_003','KK_wm_030414_004','KK_wm_030414_005','KK_wm_030414_006','KK_wm_030414_007'};
file_list{7} = {'KM_wm_041114_001','KM_wm_041114_002','KM_wm_041114_003','KM_wm_041114_004','KM_wm_041114_005','KM_wm_041114_006','KM_wm_041114_007'};
file_list{8} = {'LS_wm_041114_001','LS_wm_041114_002','LS_wm_041114_003','LS_wm_041114_004','LS_wm_041114_005','LS_wm_041114_006','LS_wm_041114_007'};
file_list{9} = {'MG2_wm_040214_001','MG2_wm_040214_002','MG2_wm_040214_003','MG2_wm_040214_004','MG2_wm_040214_005','MG2_wm_040214_006','MG2_wm_040214_007'};
file_list{10} = {'NS_wm_040414_001','NS_wm_040414_002','NS_wm_040414_003','NS_wm_040414_004','NS_wm_040414_005','NS_wm_040414_006','NS_wm_040414_007',};

% Calculate the accuracy accros all conditions
% Parses the data by trial number
for a=1:length(file_list)
    for p=1:length(file_list{a})
        load(file_list{a}{p});
        
        i=0;
        for MS = [-1 1 2] % left, right, either
            
            i=i+1;       % increment the dimension corr to task side
            
            k=0;
            for TK = {'letters' 'orientation'}
                k=k+1;
                
                % Holder variables for the 6 trials per condition
                trial = data.status(and(strcmp(trials.stimtype,TK),    and(trials.task==2,trials.mem_side==MS)));
                
                % p=run a=participant i=task side k=task
                trial1(p,a,i,k) = trial(1);
                trial2(p,a,i,k) = trial(2);
                trial3(p,a,i,k) = trial(3);
                trial4(p,a,i,k) = trial(4);
                trial5(p,a,i,k) = trial(5);
                trial6(p,a,i,k) = trial(6);
                
            end
        end
    end
end

% Average across runs then subject which is precent correct/accuracy
accTrial1 = mean(mean(trial1));
accTrial2 = mean(mean(trial2));
accTrial3 = mean(mean(trial3));
accTrial4 = mean(mean(trial4));
accTrial5 = mean(mean(trial5));
accTrial6 = mean(mean(trial6));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To calculate the K-Score first find the HR and FAR using the trials
% in the block for every run

% make an arrayfor every "trial 1,2,3,etc" for each run and participant and wheather or not it was a hit or fa
counter=0;
for a=1:length(file_list)
    for p=1:length(file_list{a})
        load(file_list{a}{p});
        i=0;
        for MS = [-1 1 2] % left, right, either
            
            i=i+1;       % increment the dimension corr to task side
            
            k=0;
            for TK = {'letters' 'orientation'}
                k=k+1;
                
                % determine the correct block and store it in an array, for
                % if they anwered correct and whether it was a hit or fa
                answerArray = data.status(and(strcmp(trials.stimtype,TK),and(trials.task==2,trials.mem_side==MS)));
                actualArray = trials.targ_change(and(strcmp(trials.stimtype,TK),and(trials.task==2,trials.mem_side==MS)));
                
                for m=1:6
                    
                    % have a counter of the new array
                    counter = counter+1;
                    
                    % col 1=corect/incorrect 2=present/not present
                    % 3=trial number 4=task 5=task side
                    trialArray(counter,1) = answerArray(m);
                    trialArray(counter,2) = actualArray(m);
                    trialArray(counter,3) = m;
                    % k=1=letters    k=2=orientation
                    trialArray(counter,4) = k;
                    % -1=left 1=right 2=either
                    trialArray(counter,5) = MS;
                    trialArray(counter,6) = a;
                    trialArray(counter,7) = p;
                    
                end
            end
        end
    end
end

for a=1:length(file_list)
    i=0;
    for MS = [-1 1 2] % left, right, either
        
        i=i+1;       % increment the dimension corr to task side
        
        k=0;
        for TK = {'letters' 'orientation'}
            k=k+1;
            
            for m=1:6
                partHR(m,i,k,a) = mean(trialArray(and((trialArray(:,2)==1),(and((trialArray(:,3)==m),(and((trialArray(:,4)==k),(and((trialArray(:,5)==MS),...
                    (trialArray(:,6)==a)))))))),1));
                partFAR(m,i,k,a) = mean(~trialArray(and((trialArray(:,2)==0),(and((trialArray(:,3)==m),(and((trialArray(:,4)==k),(and((trialArray(:,5)==MS),...
                    (trialArray(:,6)==a)))))))),1));
                
                partKVal(m,i,k,a) = 3*(partHR(m,i,k,a)-partFAR(m,i,k,a));
                
                if isnan(partHR(m,i,k,a))
                    partHR(m,i,k,a) = 0;
                end
                if isnan(partFAR(m,i,k,a))
                    partFAR(m,i,k,a) = 0;
                end
                if isnan(partKVal(m,i,k,a))
                    partKVal(m,i,k,a) = 0;
                end
            end
        end
    end
end


% Collapse over participants
meanKVal = mean(partKVal,4);
steKVal = ste(partKVal,4);


% Now determine the HR and FAR for each trial
% Calculate K score based on runs. Each participant should have 6 K scores
% per condition. Then average accross participants.

% for i=1:10
%    participantHR =
%    participantFAR =
% end
%
% % i=trial j=task MS=task side
% for i=1:6
%     % 1=letters 2=orientation
%     for j=1:2
%         % -1=left 1=right 2=either
%         k=0;
%         for MS = [-1 1 2]
%
%             k=k+1;
%             HR(i,j,k) = trialArray(and((trialArray(:,2)==1),(and((trialArray(:,3)==i),(and((trialArray(:,4)==j),((trialArray(:,5)==MS))))))));
%             FAR(i,j,k) = ~trialArray(and((trialArray(:,2)==0),(and((trialArray(:,3)==i),(and((trialArray(:,4)==j),((trialArray(:,5)==MS))))))));
%
%
%             % Finally calculate the K-val for each trial and condition
%             % SS*(HR-FAR)
%             %             KVal(i,j,k) = 3*(HR(i,j,k)-FAR(i,j,k));
%         end
%     end
% end


x = 1:6;
meany1 = meanKVal(:,1,1);
meany2 = meanKVal(:,2,1);
meany3 = meanKVal(:,3,1);

stey1 = steKVal(:,1,1);
stey2 = steKVal(:,2,1);
stey3 = steKVal(:,3,1);

figure()
subplot(1,2,1)
hold on
plot(x,meany1,x,meany2,x,meany3);
errorbar([meany1, meany2, meany3], [stey1, stey2, stey2])
str = {'','Cowan''s K Average Across Trials for Letters',''}; % cell-array method
xlabel('Trial','FontSize',15);
ylabel('Cowan''s K','FontSize',15);
set(gca,'XTick',0:7,'XLim',[0,7],'YTick',0:.5:3,'YLim',[0,3])
title(str,'FontSize',15,'FontWeight','bold');
legend('right','left','both');
hold off

x = 1:6;
meany1 = meanKVal(:,1,2);
meany2 = meanKVal(:,2,2);
meany3 = meanKVal(:,3,2);

stey1 = steKVal(:,1,2);
stey2 = steKVal(:,2,2);
stey3 = steKVal(:,3,2);

subplot(1,2,2)
hold on
plot(x,meany1,x,meany2,x,meany3);
errorbar([meany1, meany2, meany3], [stey1, stey2, stey2])
str = {'','Cowan''s K Average Across Trials for Orientation',''}; % cell-array method
xlabel('Trial','FontSize',15);
ylabel('Cowan''s K','FontSize',15);
set(gca,'XTick',0:7,'XLim',[0,7],'YTick',0:.5:3,'YLim',[0,3])
title(str,'FontSize',15,'FontWeight','bold');
legend('Right','Left','Both');
hold off


% collapsing across left, right, both (just for fun)
lrbmeanKVal = mean(meanKVal,2);
lrbsteKVal = ste(meanKVal,2);

x = 1:6;
meany1 = lrbmeanKVal(:,1,1);
meany2 = lrbmeanKVal(:,1,2);

stey1 = lrbsteKVal(:,1,1);
stey2 = lrbsteKVal(:,1,2);

figure()
hold on
plot(x,meany1,x,meany2);
errorbar([meany1, meany2], [stey1, stey2])
str = {'','Combined Cowan''s K Average Across Left, Right, and Both',''}; % cell-array method
xlabel('Trial','FontSize',15);
ylabel('Cowan''s K','FontSize',15);
set(gca,'XTick',0:7,'XLim',[0,7],'YTick',0:.5:3,'YLim',[0,3])
title(str,'FontSize',15,'FontWeight','bold');
legend('Letters','Orientated Bars');
hold off

counter = 0;
for a = 1:10
    for i = 1:2
        for j = 1:3
            for k = 1:6
                
                counter = counter + 1;
                rawdata(counter,1) = partKVal(k,j,i,a);
                rawdata(counter,2) = i;
                rawdata(counter,3) = j;
                rawdata(counter,4) = k;
                
            end
        end
    end
end

a = squeeze(reshape(partKVal,[6*3*2 1 1 10]));
a = a';







