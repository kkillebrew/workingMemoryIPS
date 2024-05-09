#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/preprocessing/
#
# N.B. There are many options when preprocessing fMRI data.  This script
#      is not meant to be general enough to capture all possible combinations
#      of volume and/or surface based preprocessing steps (e.g., detrending,
#      smoothing, signal scaling or normalization).  Be sure to consider the
#      inclusion / exclusion of each preprocessing step for your data.
#
#      This script implements:
#      1. Volume Smoothing (optional, can be set to zero)
#      2. Volume Signal Scaling (i.e., normalization)
#      3. Copy final volume files to ../analysis/ (will be input to GLM)
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# variable validation
if ( $#VOLSMOOTHS == 0) then
    echo; echo "ERROR: you must define at least one level of smoothing in global_vars.sh, even if you set it to zero (0)."
    exit 1
endif


# local variables


# for quick commenting
if ( 0 ) then
endif
set curfunc = "_tscvr";


# loop over every run
set base_curfunc = $curfunc; # we'll need to reset this back to the original value for different smoothing levels
foreach run ($ALLRUNS)
    # loop over every volume-smoothing value (may be zero)
    foreach sm ($VOLSMOOTHS)
	set curfunc = $base_curfunc; # reset curfunc for new smoothing level

        #_________________________
	# smoothing at $VOLSMOOTHS
	echo; echo "== smoothing run $run with $sm mm kernal"
	if ( $sm == 0 ) then
	    # for zero smoothing, we still want a file with "sm0" appended, just to be explicit.
	    # so although this is really of a waste of disc space, we'll make a copy of the curfunc.
	    3dcopy ${SUBJ}_${run}${curfunc}+orig ${SUBJ}_${run}${curfunc}sm$sm+orig
	else
	    3dmerge \
		-1blur_fwhm $sm \
		-doall \
		-prefix ${SUBJ}_${run}${curfunc}sm$sm \
		${SUBJ}_${run}${curfunc}+orig
	endif
	set curfunc = ${curfunc}sm$sm

        #________________________________________________________
	# signal scaling - normalization to percent signal change
	# calculate voxel-wise mean (i.e., mean over time points)
	echo; echo "== signal scaling (normalization) for run $run ($sm mm sm)"
	3dTstat \
	    -prefix ${SUBJ}_${run}${curfunc}_mean \
	    ${SUBJ}_${run}${curfunc}+orig

	# Signal Scaling (i.e., Signal Normalization)
	# there are multiple options for percent signal normalization:
	#    -expr '100*a/b'                   ; basic percent of mean (mean = 100)
	#    -expr '100*(a-b)/b'               ; basic percent of mean (mean = 0)
	#    -expr '100*a/b*ispositive(b-200)' ; masked. voxels must have mean of at least 200 across a run
	#    -expr 'min(200, 100*a/b)'         ; cap normalized signal change at 200 (twice the mean)
	#                                        avoids very large, and presumably non-physiological, signal changes.
	#                                        also keeps range symmetrical since scanner values are always positive and thus, will never go below 0%.
	# For a GLM analysis, normalization shouldn't appreciably affect significance (although there may be minor effects from masking/capping)
	# and it is meant to make the outputted beta-weights in interpretable units (percent signal change)
	3dcalc \
	    -a ${SUBJ}_${run}${curfunc}+orig \
	    -b ${SUBJ}_${run}${curfunc}_mean+orig \
	    -prefix ${SUBJ}_${run}${curfunc}_norm \
	    -expr 'min(200, 100*a/b)';
	set curfunc = ${curfunc}_norm

	# move (copy?) normalized data to analysis directory
	mv -f ${SUBJ}_${run}${curfunc}+orig.* ../analysis/
    end
end


echo; echo "== $0 complete"
exit 0
