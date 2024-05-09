#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/preprocessing/
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# local variables


# for quick commenting
if ( 0 ) then
endif
set curfunc = ""; # "" at begining of preprocessing (start of this script)


#______________________________________________________________________________
# get a copy of the original data.
# by keeping the original data in ./orig/, you can reset the preprocessing/analysis
# by deleting all files in ./analysis/, without fear of loosing the original data or
# having to re-convert from dicom.
cp ../orig/${SUBJ}_r??+orig.* ./


#______________________________________________________________________________
# de-spiking
# optional, only necessary if there are actually spikes in your dataset.
# N.B. AFNI default is to despike first.  "Most likely due to quick movement that would not be fixed by volume registration.  Despike first to minimize corruption of TRs that are adjacent to spike during time slice correction. -Rick Reynolds (10-08-08)
# N.B. A better option may be to censor bad timepoints from subsequent analyses (e.g., GLM)
#foreach run ($ALLRUNS)
#    echo; echo "== despiking run $run"
#    3dDespike \
#	-ssave ${SUBJ}_${run}_spikiness \
#	-nomask \
#	-ignore 4 \
#	-prefix ${SUBJ}_${run}${curfunc}_dspk \
#	${SUBJ}_${run}${curfunc}+orig
#end
#set curfunc = ${curfunc}_dspk; # only first preprocessing step should contain leading "_"


#----------------------
# time slice correction
# N.B. AFNI Message Board - You will have to do the slice timing correction first
#      because all the other transformations (deobliquing, motion correction, etc.) 
#      won't be able to keep track of the time each voxel was sampled. -Daniel Geln 12-13-07
# 2012.02.06 - changed from -Fourier to -quintic upon advice from AFNI gurus at Princeton AFNI bootcamp
# N.B. If you see initial scanner drift (i.e., scanner did not reach steady state and did
#      not have sufficient disdacqs), then add -ignore '4'
foreach run ($ALLRUNS)
    echo; echo "== time-slice correction for run $run"
    3dTshift \
	-slice '1' \
	-prefix ${SUBJ}_${run}_tsc \
	-quintic \
	${SUBJ}_${run}${curfunc}+orig
end
set curfunc = ${curfunc}_tsc; # only first preprocessing step should contain leading "_"


#----------------------------------------------
# volume registration (i.e., motion correction)
# N.B. do motion correction before fieldmap undistortion (and anatomical alignment)
#      so that the reference run/volume will match all other EPI volumes.
# 2012.02.06 - changed from -Fourier (default) to -cubic upon advice from AFNI gurus at Princeton AFNI bootcamp
foreach run ($ALLRUNS)
    echo; echo "== motion correction for run $run"
    3dvolreg \
       -zpad 4 \
       -prefix ${SUBJ}_${run}${curfunc}vr \
       -dfile ${SUBJ}_${run}${curfunc}vr.1D \
       -base ${SUBJ}_${REFRUN}${curfunc}+orig'['$REFBASE']' \
       -verbose \
       -cubic \
       ${SUBJ}_${run}${curfunc}+orig
end
set curfunc = ${curfunc}vr;

# concat motion params
set cat_cmd = "cat"
foreach run ($ALLRUNS)
   set cat_cmd = "$cat_cmd ${SUBJ}_${run}${curfunc}.1D"
end
set cat_cmd = "$cat_cmd > ${SUBJ}_mcparams.1D"
eval $cat_cmd; # concat all motion params together into one file

# plot motion parameters to screen...
1dplot -volreg -one ${SUBJ}_mcparams.1D'[1..6]' &
# ...and to a .png file (requires "the netpbm program pnmtopng be installed somewhere in your PATH")
1dplot -volreg -one -png ../${SUBJ}_mcparams_plot.png ${SUBJ}_mcparams.1D'[1..6]'

# move (copy?) motion params to analysis (for GLM detrending)
mv -f ${SUBJ}_mcparams.1D ../analysis/


echo; echo "== $0 complete"
exit 0
