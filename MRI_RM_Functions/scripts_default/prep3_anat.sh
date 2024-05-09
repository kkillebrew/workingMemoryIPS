#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/preprocessing/
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# local variables
set use_mean_epi_for_alignment = 0; # should we make a grand-mean EPI to use for alignment?  unlikely to make much of a difference either way, unless the individual EPI volumes are very noisy (e.g., very high-res data).  AFNI folks advise against this in most cases. default=0
set manual_nudge = 0; # should we open AFNI to allow for manual nudging of the surfvol before alignment?  otherwise, use a default "coarse" pass using 3dAllineate (see below). default=0
set deoblique_anat = 0; # should we de-oblique the anatomical as part of the coarse alignement?  sometimes this helps, but not always. default=0


# for quick commenting
if ( 0 ) then
endif
set curfunc = "_tscvr"
set curanat = "_surfvol"; # should usually be "_surfvol"; this script assumes that you have copied a version of the NoSkull_SurfVol to the ./orig/${SUBJ}_surfvol directory


#______________________________________________________________________________
# get a copy of the original data.
# by keeping the original data in ./orig/, you can reset the preprocessing/analysis
# by deleting all files in ./analysis/, without fear of loosing the original data or
# having to re-convert from dicom.
cp ../orig/${SUBJ}_surfvol+orig.* ./


#____________________________________
# get a single EPI volume to use for alignment and automask
if $use_mean_epi_for_alignment then
# get grand-average EPI
# since we've already aligned all of our functionals, we can average them all together and use
# a single average EPI frame (which visibly has better structure in some cases) for alignment with our
# anatomical.

# average over timepoints (do this first, since number of TPs may differ across runs)
echo; echo "== calculating grand mean EPI volume"
foreach run ($ALLRUNS)
3dTstat \
-prefix ${SUBJ}_${run}${curfunc}_mean+orig \
${SUBJ}_${run}${curfunc}+orig.HEAD
end

# average over runs (timepoint by timepoint, but there is only 1 TP/run now)
3dMean \
-prefix ${SUBJ}_allruns${curfunc}_mean \
${SUBJ}_r??${curfunc}_mean+orig.HEAD
set EPI_for_align = ${SUBJ}_allruns${curfunc}_mean
else
#____________________________________
# pull out reference EPI TR for alignment (same as used for motion correction)
# AFNI gurus advise against the mean EPI approach due to loss of structure from averaging
echo; echo "== extracting reference EPI volume"
3dTcat \
-prefix ${SUBJ}_refEPIvol+orig \
${SUBJ}_${REFRUN}${curfunc}+orig'['$REFBASE']'
set EPI_for_align = ${SUBJ}_refEPIvol
endif

# also create a "normalized" version of the alignment EPI for easier EPI-ANAT alignment visualization
# N.B. 3dUnifize is really designed for T1 images.  I adjust -Urad to account for the lower
#      resolution EPI voxels.  But there is no garuntee that this will work "well".  However,
#      I really just want a more uniform EPI image to judge the alignment of the EPI and ANAT,
#      so if it works for that, then the details don't really matter.
echo; echo "== unifizing EPI volume"
3dUnifize \
-prefix ${EPI_for_align}_uniform+orig \
-Urad 6.1 -GM \
${EPI_for_align}+orig


#----------------------------------------------------------
# surfvol to epi alignment
# implementing a double-alignment: first pass is "coarse" (and potentially manual), second-pass is "fine"
#
# N.B. all alignments should be rigid-body, 6-parameter, linear alignment
# N.B. If the following fails at the align_epi_anat.py stage, consider...
#      - for failure at the 3dSkullStrip stage for the EPI, try -epi_strip 3dAutomask (e.g., helps with very smaller coverage of brain)
#      - try -partial_coverage for smaller coverage
#      - see -AddEdge for some additional information to compare alignment (see align_epi_anat.py -help)

echo; echo "== starting coarse alignment (first-pass)"
if ($manual_nudge) then
# if the default automatic coarse alignment is not working well (and you don't want to spend time tweaking the default parameters)
# you can implement a manual coarse alignment using the NudgeDataset plugin of AFNI.  This code will create the file you want to
# manually nudge (_nudged), open AFNI, wait until you nudge the dataset and close AFNI, and then continue along with the fine
# pass alignment.

# make a copy of the surfvol that you will manually nudge
3dcopy ${SUBJ}_surfvol+orig ${SUBJ}_surfvol_nudged+orig.

setenv AFNI_DETACH NO; # don't let AFNI detach from terminal so we can nudge dataset before proceeding
echo; echo "*** manually nudge surfvol_nudged dataset in AFNI to be close to the $curfunc of one run.  close AFNI to continue ***"
afni
setenv AFNI_DETACH YES; # back to default
set curanat = ${curanat}_nudged
else
# the following is the default "coarse" alignment.  It uses 3dAllineate to get the surfvol close enough to the mean epi
# for the fine alignment to work.
# N.B. Previously, the coarse alignement was implemented using -giant_move of align_epi_anat.py.  But, because -giant_move
#      resets -Allineate_opts, using -giant_move utilizes a non-linear (12 parameter) warping of the anatomical, even with
#      -Allineate_opts "-warp shift_rotate"
#      Here, we are implementing something similar to -giant_move more directly, but not as elegantly and perhaps not over such
#      a large search space.
# N.B. equivilent -giant_move arguments in align_epi_anat.ph for 3dAllineate are:
#           -twobest 11 -twopass -VERB -maxrot 45 -maxshf 40 -fineblur 1 -source_automask+2
#      but using these caused poor alignment on some test datasets for reasons that I do not fully understand.
# N.B. you can emulate -giant_move and still get a rigid body transform with these options for align_epi_anat.py:
#           -Allineate_opts "-weight_frac 1.0 -maxrot 45 -maxshf 40 -VERB -warp shift_rotate" -big_move -cmass cmass
# N.B. If you ever do use -giant_move, make sure to include -master_anat SOURCE so that the anatomical is not cropped to the size
#      of the EPI, which will lead to subsequent failure of the alignment between the _al_al and SurfVol.

# sometimes deobliquing the anatomical helps.  it usually gets it closer aligned with the EPI, if the EPIs were
# collected at an oblique angle.  but this isn't always the case and usually we can deal with it without deobliquing
# the anatomy.  set deoblique_anat = 1 (above) if automatic alignment fails and surfvol and the epi are rotated far
# apart to begin with.
if ( $deoblique_anat ) then
3dWarp -verb -card2oblique ${SUBJ}_${REFRUN}+orig. -prefix ${SUBJ}${curanat}_obl+orig ${SUBJ}${curanat}+orig.
set curanat = ${curanat}_obl
endif

# the "coarse" alignment command
3dAllineate \
-prefix ${SUBJ}${curanat}_coarse+orig \
-base ${EPI_for_align}+orig \
-master ${SUBJ}${curanat}+orig \
-warp shift_rotate \
${SUBJ}${curanat}+orig
set curanat = ${curanat}_coarse
endif


# "fine" alignment
# now that the surfvol is "close" to the EPI, do another pass using the default align_epi_anat.py parameters
# N.B. Sometimes, for reasons unknown to me, this will fail miserably, and take a pretty-close alignment
#      between the EPI and the coarse-SurfVol and output something that is way off.  In that case, you can
#      try the following (but keep in mind that I haven't really had much success with the first two anyway):
#      1. If you know that the input EPI and ANAT to this fine-pass are very close (as they
#         should be - you can visually check in AFNI), then limit the range of startings points for the alignment
#         search by updating the -Allineate_opts as follows (but note that you may want/need to play with the exact
#         values for -maxrot and -maxshf):
#            -Allineate_opts "-weight_frac 1.0 -maxrot 4 -maxshf 8 -VERB -warp shift_rotate"
#      2. If you think that the alignment is failing because of a particular region (e.g., bad distortions in the
#         frontal cortex, perhaps due to a retainer or glasses), then you can try to provide an exclusion mask
#         to 3dAllineate with:
#            -Allineate_opts "-weight_frac 1.0 -maxrot 6 -maxshf 10 -VERB -warp shift_rotate -emask KK_surfvol_obl_coarse_exclude_frontal_mask+orig"
#         You can create such a mask using the DrawDataset plugin in the AFNI gui. Should be the same size as the
#         source (-anat) dataset in the command below.  If you draw it on the EPI, try using 3dresample to convert
#         to the appropriate grid/resolution.
#      3. Rely on a purely manual alignment using NudgeDataset in the AFNI gui.  You can set manual_nudge=1 (above)
#         and comment out the "fine" pass code.  In that case, the manual nudging will be the "final" alignmen
#         between EPI and ANAT
#
# N.B. the default -Allineate_opts to implement a rigid body transformation (and keep other align_epi_anat.py defaults) is:
#           -Allineate_opts "-weight_frac 1.0 -maxrot 6 -maxshf 10 -VERB -warp shift_rotate"
echo; echo "== starting fine alignment (second-pass)"
align_epi_anat.py \
-anat ${SUBJ}${curanat}+orig \
-epi ${EPI_for_align}+orig \
-epi_base 0 \
-volreg off \
-tshift off \
-deoblique off \
-anat_has_skull no \
-Allineate_opts "-weight_frac 1.0 -maxrot 6 -maxshf 10 -VERB -warp shift_rotate"
set curanat = ${curanat}_al



#______________________________________________________________________________
# align high-res surface anatomical to epi-aligned anatomical
# N.B. align_epi_anat.py skull-strips the anatomical, so use NoSkull hi-res for surface alignment
# output is ${SUBJ}_[NoSkull_]SurfVol_Alnd_Exp+orig
echo; echo "== aligning surface volume"
@SUMA_AlignToExperiment \
-align_centers \
-strip_skull neither \
-surf_anat ${SUMADIR}/${SUBJ}_NoSkull_SurfVol+orig \
-exp_anat ${SUBJ}${curanat}+orig

# copy experiment-aligned surface volumne to ./analysis
cp -f ${SUBJ}_NoSkull_SurfVol_Alnd_Exp+orig.* ../analysis/



#______________________________________________________________________________
# create automask - just after motion correction/undistortion AND anatomical alignment,
# BUT before smoothing/filtering/normalization etc.
# this can be used to speed up 3dDeconvolve by ignoring the out-of-brain voxels
# create binary mask from  EPI brik
echo; echo "== creating automask"
3dAutomask \
-prefix ${SUBJ}_automask \
-dilate 3 \
${EPI_for_align}+orig

# copy automask and EPI_for_align for checking alignment (now that we've used it for the automask)
cp -f ${SUBJ}_automask+orig.* ../analysis/
cp -f ${EPI_for_align}+orig.* ../analysis/
cp -f ${EPI_for_align}_uniform+orig.* ../analysis/


echo; echo "== $0 complete"
exit 0
