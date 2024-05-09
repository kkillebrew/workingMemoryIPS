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
#      1. Volume Detrending (optional, can be set to zero)
#         (this is in preparation for a TimeCourse analysis, so typically no smoothing as
#          smoothing will occur during ROI-extraction and averaging)
#      2. Copy final volume files to ../analysis/ (will be input to TC analysis in Matlab)
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# local variables
# TODO: polort -> global_vars.sh as POLORT
set polort = 1; # polort level for 3dDetrend (1 = linear; 2 = linear+quadratic; etc.)
set norm_dtm_data = 1; # boolean, should we normalize (i.e., percent signal change) the detrended data?  useful for timecourse analysis in which data for each voxel is normalized to each run's mean (instead of a baseline period before each block/event)


# for quick commenting
if ( 0 ) then
endif
set curfunc = "_tscvr";



#______________________________________________________________________________
# DETREND DATA - be sure to add back in mean for accuracy
echo; echo "== Detrending";
foreach run ($ALLRUNS)
    echo; echo "== detrending run ${run}"
    # Detrending
    3dDetrend \
	-prefix ${SUBJ}_${run}${curfunc}dt \
	-polort $polort \
	${SUBJ}_${run}${curfunc}+orig

    # add back in mean for accuracy (_mean may already exist from concurrently running scripts)
    if ( ! -e ${SUBJ}_${run}${curfunc}_mean.HEAD ) then
	3dTstat \
	    -prefix ${SUBJ}_${run}${curfunc}_mean \
	    ${SUBJ}_${run}${curfunc}+orig
    endif
		
    3dcalc \
	-a ${SUBJ}_${run}${curfunc}dt+orig \
	-b ${SUBJ}_${run}${curfunc}_mean+orig \
	-prefix ${SUBJ}_${run}${curfunc}dtm \
	-expr "a+b"
end
set curfunc = ${curfunc}dtm;



#________________________________________________________
# signal scaling - normalization to percent signal change
# calculate voxel-wise mean (i.e., mean over time points)
if ( $norm_dtm_data ) then
    foreach run ($ALLRUNS)
	echo; echo "== signal scaling (normalization) run $run"
	
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
    end
    set curfunc = ${curfunc}_norm;
endif



foreach run ($ALLRUNS)
    # move (copy?) detrended data to analysis directory
    mv -f ${SUBJ}_${run}${curfunc}+orig.* ../analysis/
end

echo; echo "== $0 complete"
exit 0
