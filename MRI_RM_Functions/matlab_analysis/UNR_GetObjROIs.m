function objROIs = UNR_GetObjectROIs(subj_initials,format)
% UNR_GetObjectROIs(subj_initials)
%   objROIs = UNR_GetObjectROIs(subj_initials)
%
% Load object ROIs from disk.  Return as structure objROIs.
%
% Arguments
%   subj_initials = initials of subject
%
% Output
%   objROIs = return structure.

% REHBM 10.07.08
%       04.12.10 - added 'names' support
%       05.13.10 - added quick switch (internal) for using strict ROI definition
%       04.2014  - updated for UNR scripts


%% validate arguments
if nargin < 2
    format = 'vector';
end


%% define indicies and return names if requested
objROIs.idx.FFA = 1; %index into ROI BRIK
objROIs.idx.PPA = 2;
objROIs.idx.LOC = 3;
objROIs.idx.PFS = 4;
objROIs.idx.OFA = 5;
objROIs.idx.hMT = 6;
objROIs.idx.OBJ = [objROIs.idx.FFA objROIs.idx.PPA objROIs.idx.LOC objROIs.idx.PFS objROIs.idx.OFA];
objROIs.idx.FFAPPA = [objROIs.idx.FFA objROIs.idx.PPA];

if strcmp(format,'names')
    % just return the ROIs defined by this function
    objROIs = fieldnames(objROIs.idx);
    return
end


%% get data
datadir = ['../' subj_initials '/analysis/'];

% object ROIs
suffix = '_objectROIs';


[objROIs.err.lh, objROIs.data.lh, objROIs.info.lh, objROIs.errMessage.lh] = BrikLoad([datadir sdb.subj suffix '_lh+orig'],format);
if objROIs.err.lh
    error('error loading lh object ROIs file (see above).')
end
[objROIs.err.rh, objROIs.data.rh, objROIs.info.rh, objROIs.errMessage.rh] = BrikLoad([datadir sdb.subj suffix '_rh+orig'],format);
if objROIs.err.rh
    error('error loading rh object ROIs file (see above).')
end
