% Gets the time for the GLM for wmIPS

clear all;
close all;

subID = 'LS';
% pathOut = '/Users/C-Lab/Google Drive/Lab Projects/Marians Stuff/wmIPS/Experiment/Timing/';
pathOut = sprintf('/Volumes/imaginguser/MR_DATA/BerryhillLab/WMIPS/%s/scripts/stim_times/',subID);
pathIn = '/Users/C-Lab/Google Drive/Lab Projects/Marians Stuff/wmIPS/Experiment/data/';

% file_list = {'GG_wm_030414_001.mat','GG_wm_030414_002.mat','GG_wm_030414_003.mat','GG_wm_030414_004.mat','GG_wm_030414_005.mat','GG_wm_030414_006.mat','GG_wm_030414_007.mat'};
% file_list = {'KK_wm_030414_001.mat','KK_wm_030414_002.mat','KK_wm_030414_003.mat','KK_wm_030414_004.mat','KK_wm_030414_005.mat','KK_wm_030414_006.mat','KK_wm_030414_007.mat'};
% file_list = {'CDB_wm_031414_001.mat','CDB_wm_031414_002.mat','CDB_wm_031414_003.mat','CDB_wm_031414_004.mat','CDB_wm_031414_005.mat','CDB_wm_031414_006.mat','CDB_wm_031414_007.mat'};
% file_list = {'JEV_wm_031414_001.mat','JEV_wm_031414_002.mat','JEV_wm_031414_003.mat','JEV_wm_031414_004.mat','JEV_wm_031414_005.mat','JEV_wm_031414_006.mat','JEV_wm_031414_007.mat'};
% file_list = {'GPC_wm_031814_001.mat','GPC_wm_031814_002.mat','GPC_wm_031814_003.mat','GPC_wm_031814_004.mat','GPC_wm_031814_005.mat','GPC_wm_031814_006.mat','GPC_wm_031814_007.mat'};
% file_list = {'MG2_wm_040214_002.mat','MG2_wm_040214_003.mat','MG2_wm_040214_004.mat','MG2_wm_040214_005.mat','MG2_wm_040214_006.mat','MG2_wm_040214_007.mat',};
% file_list = {'NS_wm_040414_001','NS_wm_040414_002','NS_wm_040414_003','NS_wm_040414_004','NS_wm_040414_005','NS_wm_040414_006','NS_wm_040414_007',};
% file_list = {'KM_wm_041114_001','KM_wm_041114_002','KM_wm_041114_003','KM_wm_041114_004','KM_wm_041114_005','KM_wm_041114_006','KM_wm_041114_007'};
file_list = {'LS_wm_041114_001','LS_wm_041114_002','LS_wm_041114_003','LS_wm_041114_004','LS_wm_041114_005','LS_wm_041114_006','LS_wm_041114_007'};

for a = 1:length(file_list)
    load(sprintf('%s%s',pathIn,file_list{a}));
    % Get the times of each block in a run
    t1 = timing.block_cue_blank_stop - timing.run_start;
    t2 = timing.block_cue_start - timing.run_start;
    
    % Time of orientation left non passive
    oriLeft(a,1) = t1(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==-1));
    % Time of orientation right non passive
    oriRight(a,1) = t1(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==1));
    % Time of orientation both non passive
    oriBoth(a,1) = t1(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==2));
    
    % Time of letters left non passive
    letLeft(a,1) = t1(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==-1));
    % Time of letters right non passive
    letRight(a,1) = t1(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==1));
    % Time of letters both non passive
    letBoth(a,1) = t1(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==2));
    
    % Times of passive letters
    letPassive(a,1) = t1(and(strcmp(blocks.stimtype,'letters'),blocks.task==1));
    % Times of passive orient
    oriPassive(a,1) = t1(and(strcmp(blocks.stimtype,'orientation'),blocks.task==1));
    
    % Times of cue onset for orientation left non passive
    oriLeftCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==-1));
    % Times of cue onset for orientation right non passive
    oriRightCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==1));
    % Times of cue onset for orientation both non passive
    oriBothCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'orientation'),blocks.task==2),blocks.mem_side==2));
    
    % Times of cue onset for letters left non passive cue
    letLeftCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==-1));
    % Times of cue onset for letters right non passive cue
    letRightCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==1));
    % Times of cue onset for letters both non passive cue
    letBothCue(a,1) = t2(and(and(strcmp(blocks.stimtype,'letters'),blocks.task==2),blocks.mem_side==2));
    
    % Times of passive letters cue
    letPassiveCue(a,1) = t2(and(strcmp(blocks.stimtype,'letters'),blocks.task==1));
    % Times of passive orient cue
    oriPassiveCue(a,1) = t2(and(strcmp(blocks.stimtype,'orientation'),blocks.task==1));
end

    % Writing to the individual files for use with the glm functions
    fileID = fopen(sprintf('%s%s_ori_left.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriLeft);
    fclose(fileID);
    fileID = fopen(sprintf('%s%s_ori_right.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriRight);
    fclose(fileID);
    fileID = fopen(sprintf('%s%s_ori_both.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriBoth);
    fclose(fileID);
    
    fileID = fopen(sprintf('%s%s_let_left.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letLeft);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_let_right.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letRight);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_let_both.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letBoth);
    fopen(fileID);
    
    fileID = fopen(sprintf('%s%s_let_pass.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letPassive);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_ori_pass.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriPassive);
    fopen(fileID);
    
    fileID = fopen(sprintf('%s%s_let_left_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letLeftCue);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_let_right_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letRightCue);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_let_both_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letBothCue);
    fopen(fileID);
    
    fileID = fopen(sprintf('%s%s_ori_left_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriLeftCue);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_ori_right_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriRightCue);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_ori_both_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriBothCue);
    fopen(fileID);
    
    fileID = fopen(sprintf('%s%s_let_pass_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',letPassiveCue);
    fopen(fileID);
    fileID = fopen(sprintf('%s%s_ori_pass_cue.1D',pathOut,subID),'w');
    fprintf(fileID,'%f\n',oriPassiveCue);
    fopen(fileID);