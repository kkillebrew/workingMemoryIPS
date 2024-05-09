function suffix = UNR_GetROIsSuffix()
% UNR_GetROIsSuffix()
%   suffix = UNR_GetROIsSuffix
%
% Shortcut for supplying same suffix for volume ROI files to all the GetXXXROIs functions.
%
% Arguments
%   none
%
% Output
%   suffix = suffix defining which volume ROI file to load.
%            e.g. ['_retROIs' suffix];
%                 suffix = '' or '_2mm' or '_unique100', etc

% REHBM 01.30.11 - created for use by all GetXXXROIs functions

suffix = '';          % straight surface->volume projection using MODE mapping function (as of Dec 2011)
%suffix = '_2mm';      % "strict" cortical definition based on using -f_p1_mm 2.0 and -f_pn_mm -2.0 as part of 3dSurf2Vol
%suffix = '_uniqueROI50'; % "unique" mapping function where MODE is accepted if it accounts for 50% of contributing nodes (ie, nodes with a defined ROI in the current file)
%suffix = '_uniqueCortical50'; % "unique" mapping function where MODE is accepted if it accounts for 50% of contributing (ie, nodes with a defined ROI in the current file) and cortical nodes
% etc.

%fprintf('+++ using ROI suffix = %s\n',suffix)