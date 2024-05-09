clear all;
close all;

% file_list={'KWK_orient_pilot_wm_020514_001', 'CDB_orient_pilot_wm_020514_001','GG_orient_pilot_wm_020614_001','DM_orient_pilot_wm_021114_001','JEV_orient_pilot_wm_021914_001','RS_orient_pilot_wm_021914_001'};
% file_list={'KWK_letters_pilot_wm_021114_001','CDB_letters_pilot_wm_021114_001','GG_letters_pilot_wm_021114_001', 'DM_letters_pilot_wm_021114_001'};
% file_list_letters={'KWK_lettersRK_pilot_wm_021214_001','CDB_lettersKR_pilot_wm_021214_001','GG_lettersKR_pilot_wm_021214_001','DM_lettersRK_pilot_wm_021214_001','JEV_letters_pilot_wm_021914_001','RS_letters_pilot_wm_021914_001'};
file_list{1} = {'KK_wm_030414_001','KK_wm_030414_002','KK_wm_030414_003','KK_wm_030414_004','KK_wm_030414_005','KK_wm_030414_006','KK_wm_030414_007'};
file_list{2} = {'GG_wm_030414_001','GG_wm_030414_002','GG_wm_030414_003','GG_wm_030414_004','GG_wm_030414_005','GG_wm_030414_006','GG_wm_030414_007'};

% 1 = oldCorrect
% 2 = oldTotal
% 3 = newCorrect
% 4 = newTotal

% 1 = left
% 2 = right
% 3 = both

% targ_change = 1   means it was old (the target changed)
% targ_change = 2   means it was new (the target did not change)
% data_status = 0   means response was incorrect
% data_status = 1   means response was correct

countOrient = zeros(3,4,2);
countLetters = zeros(3,4,2);
a=1;

for a=1:length(file_list)
    for p=1:length(file_list)
        load(file_list{a}{p});
        % Determine for both old and new trial types accuracy rate
        for i=1:trials.n
            if strcmp(trials.stimtype{i},'letters')
                % Make sure passive viewing is not included
                if trials.task(i)==2
                    % Determines which condition you are in
                    if trials.mem_side(i)==-1
                        % Was it a new or an old trial
                        if trials.targ_change(i)==1
                            countLetters(1,4,a) = countLetters(1,4,a)+1;        % Old trial
                            % Did you respond correctly
                            if data.status(i)==1
                                countLetters(1,3,a) = countLetters(1,3,a)+1;     % Responded it old
                            end
                        else
                            countLetters(1,2,a) = countLetters(1,2,a)+1;     % New trial
                            if data.status(i)==0
                                countLetters(1,1,a) = countLetters(1,1,a)+1;  % Responded it was old
                            end
                        end
                    elseif trials.mem_side(i)==1
                        if trials.targ_change(i)==1
                            countLetters(2,4,a) = countLetters(2,4,a)+1;
                            if data.status(i)==1
                                countLetters(2,3,a) = countLetters(2,3,a)+1;
                            end
                        else
                            countLetters(2,2,a) = countLetters(2,2,a)+1;
                            if data.status(i)==0
                                countLetters(2,1,a) = countLetters(2,1,a)+1;
                            end
                        end
                    elseif trials.mem_side(i)==2
                        if trials.targ_change(i)==1
                            countLetters(3,4,a) = countLetters(3,4,a)+1;
                            if data.status(i)==1
                                countLetters(3,3,a) = countLetters(3,3,a)+1;
                            end
                        else
                            countLetters(3,2,a) = countLetters(3,2,a)+1;
                            if data.status(i)==0
                                countLetters(3,1,a) = countLetters(3,1,a)+1;
                            end
                        end
                    end
                end
            end
        end
        
        
        % Determine for both old and new trial types accuracy rate
        for i=1:trials.n
            if strcmp(trials.stimtype{i},'orientation')
                % Make sure passive viewing is not included
                if trials.task(i)==2
                    % Determines which condition you are in
                    if trials.mem_side(i)==-1
                        % Was it a new or an old trial
                        if trials.targ_change(i)==1
                            countOrient(1,4) = countOrient(1,4)+1;
                            % Did you respond correctly
                            if data.status(i)==1
                                countOrient(1,3) = countOrient(1,3)+1;
                            end
                        else
                            countOrient(1,2) = countOrient(1,2)+1;
                            if data.status(i)==0
                                countOrient(1,1) = countOrient(1,1)+1;
                            end
                        end
                    elseif trials.mem_side(i)==1
                        if trials.targ_change(i)==1
                            countOrient(2,4) = countOrient(2,4)+1;
                            if data.status(i)==1
                                countOrient(2,3) = countOrient(2,3)+1;
                            end
                        else
                            countOrient(2,2) = countOrient(2,2)+1;
                            if data.status(i)==0
                                countOrient(2,1) = countOrient(2,1)+1;
                            end
                        end
                    elseif trials.mem_side(i)==2
                        if trials.targ_change(i)==1
                            countOrient(3,4) = countOrient(3,4)+1;
                            if data.status(i)==1
                                countOrient(3,3) = countOrient(3,3)+1;
                            end
                        else
                            countOrient(3,2) = countOrient(3,2)+1;
                            if data.status(i)==0
                                countOrient(3,1) = countOrient(3,1)+1;
                            end
                        end
                    end
                end
            end
        end
    end
end