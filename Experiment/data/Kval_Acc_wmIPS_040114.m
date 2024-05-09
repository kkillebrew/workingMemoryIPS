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


xLabelsK={'Left','Right','Either'};
xLabelsAcc={'Orientation', 'Letters'};

for a=1:length(file_list)
    for p=1:length(file_list{a})
        load(file_list{a}{p});
        
        % Accuracy overall for letters and orientation
        accOri(a,p,1) = 100 * mean(data.status(and(trials.task==2,strcmp(trials.stimtype,'orientation'))));
        accLet(a,p,1) = 100 * mean(data.status(and(trials.task==2,strcmp(trials.stimtype,'letters'))));
        
        MS = 0;
        for i=2:4
            if i==2
                MS = 1;
            elseif i==3
                MS = -1;
            elseif i==4
                MS = 2;
            end
            % Accuracy for left right either
            accOri(a,p,i) = 100 * mean(data.status(and(trials.mem_side==MS,and(trials.task==2,strcmp(trials.stimtype,'orientation')))));
            accLet(a,p,i) = 100 * mean(data.status(and(trials.mem_side==MS,and(trials.task==2,strcmp(trials.stimtype,'letters')))));
        end
        
        % Calc HR and FAR
        % 1=left 2=right 3=both
        MS = 0;
        for i=1:3
            if i==1
                MS = 1;
            elseif i==2
                MS = -1;
            elseif i==3
                MS = 2;
            end
            HRLet(a,p,i) = mean(data.status(and(strcmp(trials.stimtype,'letters'),and(trials.task==2,and(trials.mem_side==MS,trials.targ_change==1)))));
            FARLet(a,p,i) = mean(~data.status(and(strcmp(trials.stimtype,'letters'),and(trials.task==2,and(trials.mem_side==MS,trials.targ_change==0)))));
            HROri(a,p,i) = mean(data.status(and(strcmp(trials.stimtype,'orientation'),and(trials.task==2,and(trials.mem_side==MS,trials.targ_change==1)))));
            FAROri(a,p,i) = mean(~data.status(and(strcmp(trials.stimtype,'orientation'),and(trials.task==2,and(trials.mem_side==MS,trials.targ_change==0)))));
            
            if isnan(HRLet(a,p,i))
                HRLet(a,p,i) = 0;
            end
            if isnan(FARLet(a,p,i))
                FARLet(a,p,i) = 0;
            end
            if isnan(HROri(a,p,i))
                HROri(a,p,i) = 0;
            end
            if isnan(FAROri(a,p,i))
                FAROri(a,p,i) = 0;
            end
        end
        
        % Calculating K for (1=left 2=right 3=both)
        for i=1:3
            KLet(a,i) = params.blocks.set_size*(mean(HRLet(a,:,i)-mean(FARLet(a,:,i))));
            KOri(a,i) = params.blocks.set_size*(mean(HROri(a,:,i)-mean(FAROri(a,:,i))));
        end
        
    end
    % Plot the means for each participant
    figure()
    subplot(1,3,1)
    c = [mean(accOri(a,:,1)) mean(accLet(a,:,1))];
    bar(c);
    str = {'',sprintf('Accuracy for Participant %d',a),''}; % cell-array method
    xlabel('Stimulus Type','FontSize',15);
    ylabel('Accuracy (%)','FontSize',15);
    set(gca, 'XTickLabel',xLabelsAcc, 'XTick',1:numel(xLabelsAcc),'ylim',[0,100])
    title(str,'FontSize',15,'FontWeight','bold');
    
    subplot(1,3,2)
    c = [KOri(a,1), KLet(a,1); KOri(a,2), KLet(a,2); KOri(a,3), KLet(a,3)];
    bar(c);
    str = {'',sprintf('Cowan K for Participant %d',a),''}; % cell-array method
    xlabel('Stimulus Type','FontSize',15);
    ylabel('Cowans K','FontSize',15);
    set(gca, 'XTickLabel',xLabelsK, 'XTick',1:numel(xLabelsK),'ylim',[0,3])
    title(str,'FontSize',15,'FontWeight','bold');
    legend(xLabelsAcc);
    
    subplot(1,3,3)
    c = [mean(accOri(a,:,2)) mean(accLet(a,:,2)); mean(accOri(a,:,3)) mean(accLet(a,:,3)); mean(accOri(a,:,4)) mean(accLet(a,:,4))];
    bar(c);
    str = {'',sprintf('Accuracy for Memory Side for Participant %d',a),''}; % cell-array method
    xlabel('Stimulus Type','FontSize',15);
    ylabel('Accuracy (%)','FontSize',15);
    set(gca, 'XTickLabel',xLabelsK, 'XTick',1:numel(xLabelsK),'ylim',[0,100])
    title(str,'FontSize',15,'FontWeight','bold');
    legend(xLabelsAcc);
    
end

meanAccOri = mean(mean(accOri(:,:,1),2));
stdAccOri = std(mean(accOri(:,:,1),2))/length(mean(accOri(:,:,1),2));
meanAccLet = mean(mean(accLet(:,:,1),2));
stdAccLet = std(mean(accLet(:,:,1),2))/length(mean(accLet(:,:,1),2));

meanLeftOri = mean(mean(accOri(:,:,2),2));
stdLeftOri = std(mean(accOri(:,:,2),2))/length(mean(accOri(:,:,2),2));
meanLeftLet = mean(mean(accLet(:,:,2),2));
stdLeftLet = std(mean(accLet(:,:,2),2))/length(mean(accLet(:,:,2),2));

meanRightOri = mean(mean(accOri(:,:,3),2));
stdRightOri = std(mean(accOri(:,:,3),2))/length(mean(accOri(:,:,3),2));
meanRightLet = mean(mean(accLet(:,:,3),2));
stdRightLet = std(mean(accLet(:,:,3),2))/length(mean(accLet(:,:,3),2));

meanBothOri = mean(mean(accOri(:,:,4),2));
stdBothOri = std(mean(accOri(:,:,4),2))/length(mean(accOri(:,:,4),2));
meanBothLet = mean(mean(accLet(:,:,4),2));
stdBothLet = std(mean(accLet(:,:,4),2))/length(mean(accLet(:,:,4),2));

meanKLeftOri = mean(KOri(:,1));
stdKLeftOri = std(KOri(:,1))/length(KOri(:,1));
meanKLeftLet = mean(KLet(:,1));
stdKLeftLet = std(KLet(:,1))/length(KLet(:,1));

meanKRightOri = mean(KOri(:,2));
stdKRightOri = std(KOri(:,2))/length(KOri(:,2));
meanKRightLet = mean(KLet(:,2));
stdKRightLet = std(KLet(:,2))/length(KLet(:,2));

meanKBothOri = mean(KOri(:,3));
stdKBothOri = std(KOri(:,3))/length(KOri(:,3));
meanKBothLet = mean(KLet(:,3));
stdKBothLet = std(KLet(:,3))/length(KLet(:,3));

% Plot the accuracy means
figure()
subplot(1,3,1)
c = [meanAccOri meanAccLet];
b = [stdAccOri stdAccLet];
bar(c);
hold on
errorbar(c,b,'k.','LineWidth',2);
str = {'','Average Accuracy',''}; % cell-array method
xlabel('Stimulus Type','FontSize',15);
ylabel('Accuracy (%)','FontSize',15);
set(gca, 'XTickLabel',xLabelsAcc, 'XTick',1:numel(xLabelsAcc),'ylim',[0,100])
title(str,'FontSize',15,'FontWeight','bold');

% Plot K Val means
subplot(1,3,2)
c = [meanKLeftOri, meanKLeftLet; meanKRightOri, meanKRightLet; meanKBothOri, meanKBothLet];
b = [stdKLeftOri, stdKLeftLet; stdKRightOri, stdKRightLet; stdKBothOri, stdKBothLet];
h = bar(c);
% Code from barweb.m to align the error bars with the bar groups
set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'XTicklabel','Modelo1|Modelo2|Modelo3')
set(get(gca,'YLabel'),'String','U')
hold on;
numgroups = size(c, 1); 
numbars = size(c, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, c(:,i), b(:,i), '.k', 'linewidth', 2);
end
%%%%%%%%
str = {'','Cowan K Average',''}; % cell-array method
xlabel('Stimulus Type','FontSize',15);
ylabel('Cowans K','FontSize',15);
set(gca, 'XTickLabel',xLabelsK, 'XTick',1:numel(xLabelsK),'ylim',[0,3])
title(str,'FontSize',15,'FontWeight','bold');
legend(xLabelsAcc);

% Plot Acc for mem side
subplot(1,3,3)
c = [meanLeftOri meanLeftLet; meanRightOri meanRightLet; meanBothOri meanBothLet];
b = [stdLeftOri stdLeftLet;stdRightOri stdRightLet;stdBothOri stdBothLet];
h=bar(c);
% Code from barweb.m to align the error bars with the bar groups
set(h,'BarWidth',1);    % The bars will now touch each other
set(gca,'XTicklabel','Modelo1|Modelo2|Modelo3')
set(get(gca,'YLabel'),'String','U')
hold on;
numgroups = size(c, 1); 
numbars = size(c, 2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
      x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
      errorbar(x, c(:,i), b(:,i), '.k', 'linewidth', 2);
end
%%%%%%%%
str = {'','Average Accuracy for Memory Side',''}; % cell-array method
xlabel('Stimulus Type','FontSize',15);
ylabel('Accuracy (%)','FontSize',15);
set(gca, 'XTickLabel',xLabelsK, 'XTick',1:numel(xLabelsK),'ylim',[0,100])
title(str,'FontSize',15,'FontWeight','bold');
legend(xLabelsAcc);





