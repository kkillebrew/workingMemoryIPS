function rROIs = UNR_GetRetROIs(subj_initials,format)
% UNR_GetRetROIs(subj_initials)
%   rROIs = UNR_GetRetROIs(subj_initials)
%
% Load retinotopic ROIs from disk.  Return as structure rROIs.
%
% Arguments
%   subj_initials = initials of subject
%
% Output
%   rROIs = return structure.

% REHBM 10.07.08
%       04.12.10 - added 'names' support
%       05.13.10 - added quick switch (internal) for using strict ROI definition
%       04.2014  - updated for UNR scripts


%% validate arguments
if nargin < 2
    format = 'vector';
end


%% define indicies and return names if requested
rROIs.idx.V1    = 1; %index into ROI BRIK
rROIs.idx.V2v   = 2;
rROIs.idx.V2d   = 3;
rROIs.idx.V2    = [rROIs.idx.V2v rROIs.idx.V2d];
rROIs.idx.V1_V2 = [rROIs.idx.V1 rROIs.idx.V2];
rROIs.idx.V3v   = 4;
rROIs.idx.V3d   = 5;
rROIs.idx.V3    = [rROIs.idx.V3v rROIs.idx.V3d];
rROIs.idx.EVC   = [rROIs.idx.V1 rROIs.idx.V2 rROIs.idx.V3];

rROIs.idx.V3a   = 150;
rROIs.idx.V3b   = 151;
rROIs.idx.IPS0  = 152;
rROIs.idx.IPS1  = 153;
rROIs.idx.IPS2  = 154;
rROIs.idx.IPS3  = 155;
rROIs.idx.IPS4  = 156;
rROIs.idx.IPS5  = 157;
rROIs.idx.SPL1  = 158;
rROIs.idx.pIPS  = [rROIs.idx.IPS0 rROIs.idx.IPS1 rROIs.idx.IPS2];
rROIs.idx.mIPS  = [rROIs.idx.IPS3 rROIs.idx.IPS4];

rROIs.idx.hFEF  = 200;

rROIs.idx.LO1   = 100;
rROIs.idx.LO2   = 101;
rROIs.idx.LO    = [rROIs.idx.LO1 rROIs.idx.LO2];
rROIs.idx.TO1   = 102;
rROIs.idx.TO2   = 103;
rROIs.idx.TO    = [rROIs.idx.TO1 rROIs.idx.TO2];

rROIs.idx.V4    = 50;
rROIs.idx.VO1   = 51;
rROIs.idx.VO2   = 52;
rROIs.idx.VO    = [rROIs.idx.VO1 rROIs.idx.VO2];
rROIs.idx.PHC1  = 53;
rROIs.idx.PHC2  = 54;
rROIs.idx.PHC   = [rROIs.idx.PHC1 rROIs.idx.PHC2];

% rROIs.idx.ALL_v = [rROIs.idx.V1 rROIs.idx.V2 rROIs.idx.V3 rROIs.idx.V4 rROIs.idx.V3a rROIs.idx.V7];


if strcmp(format,'names')
    % just return the ROIs defined by this function
    rROIs = fieldnames(rROIs.idx);
    rROIs = rROIs(~ismember(rROIs,{'V1v' 'V1d' 'V2v' 'V2d' 'V1_V2' 'V3v' 'V3d' 'EVC' 'LO' 'TO' 'VO' 'PHC'}));
    return
end


%% get data
datadir = ['../' subj_initials '/analysis/'];

% retinotopy ROIs
suffix = '_retROIs';

[rROIs.err.lh, rROIs.data.lh, rROIs.info.lh, rROIs.errMessage.lh] = BrikLoad([datadir subj_initials suffix '_lh+orig'],format);
if rROIs.err.lh
    error('error loading lh visual ROIs file (see above).')
end
[rROIs.err.rh, rROIs.data.rh, rROIs.info.rh, rROIs.errMessage.rh] = BrikLoad([datadir subj_initials suffix '_rh+orig'],format);
if rROIs.err.rh
    error('error loading rh visual ROIs file (see above).')
end
