function rois = UNR_ROIArgs(rois)
% ROIArgs
%   rois = ROIArgs(rois)
%
% Extract rois from an rois argument for rawMVPA, glmMVPA, and MVPA_Summary, etc.
% Provides centralized code for multiple functions.
%
% Arguments
%   rois = 'visual', 'parietal',  'all' (default), or cell array with list of specific rois
%
% Returns
%   rois = expanded roi cell array

% REHBM 03.05.09 - pulled code from rawMVPA, glmMVPA and MVPA_Summary


%% validate arguments



%% setup rois for each set of visual areas
visual_rois   = UNR_GetRetROIs(-1,'names');
% % parietal_rois = UNR_GetParietalROIs(-1,'names');
% % frontal_rois  = UNR_GetFrontalROIs(-1,'names');
object_rois   = UNR_GetObjROIs(-1,'names');

% % % may not exist, defined locally for each experiment
% % try
% %     custom_rois = GetCustomROIs(-1,'names');
% % catch
% %     custom_rois = {};
% % end

% % % may not exist, must be hand-drawn for each experiment
% % try
% %     noise_rois = GetNoiseROIs(-1,'names');
% % catch
% %     noise_rois = {};
% % end


%% expand rois argument into cell array
if ischar(rois)
    switch rois
        case 'all'
            %rois = cat(1,visual_rois,parietal_rois,frontal_rois,object_rois,custom_rois,noise_rois);
            rois = cat(1,visual_rois,object_rois);
        case 'visual'
            rois = visual_rois;
% %         case 'parietal'
% %             rois = parietal_rois;
% %         case 'frontal'
% %             rois = frontal_rois;
        case 'object'
            rois = object_rois;
% %         case 'custom'
% %             rois = custom_rois;
% %         case 'noise'
% %             rois = noise_rois;
        otherwise
            % go with what we were given
            rois = {rois};
    end
end