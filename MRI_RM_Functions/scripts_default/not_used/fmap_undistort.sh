#!/bin/tcsh

# fmap_undistort.sh
# expected to be run from a tmp directory one level above the data to be corrected
#   e.g., .../${SUBJ}/analysis/fmap/

##---------------------------------------
## So we can share this script with a variety of tasks/pipelines, we need to 
## input the task as an argument so we can source the appropriate _vars.sh script
#if ($#argv != 1) then
#    echo "Usage: $0 task"
#    exit
#endif
#set task_arg = $1;


#---------------------------------------
# 0. Local Variables
# N.B. This script assumes that FSLOUTPUTTYPE = NIFTI_GZ

# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../../scripts/global_vars.sh; # should give us variables: SUBJ, SCANNER, REFRUN, UNWARPDIR, FMAP_SM, etc.
set curfunc = "_tscvr"


# boolean, do you want to run this script with interactive checking/editing?  This is highly suggested...
set interactive = 1


# N.B. the $UNWARPDIR and $FMAP_SM parameters, defined in {experiment}_vars.sh will determine the final epi-aligned phase map's filename
#     e.g., phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz



# dwelll time (echo spacing) in s
switch ($SCANNER)
    case allegra:
	set dwell_s  = 0.00069; # dwelll time (echo spacing)in sec
	breaksw
    case skyra:
	set dwell_s  = 0.00075; # dwelll time (echo spacing)in sec
	breaksw
    case UCDavis_skyra:
	set dwell_s  = 0.00075; # dwelll time (echo spacing) in sec (*** NEED TO CHECK THIS ***)
	breaksw
    default:
	echo "ERROR: SCANNER must be 'skyra', 'allegra', or 'UCDavis_Skyra'"
	exit 1
endsw
# N.B. the "effective" echo spacing must take into account any acceleration factor
#      dwell_time_sec = (EPI_echo_spacing_ms / IPAT_FACTOR) / 1000
# N.B. to use the MATH command, add the following lines (excluding first # of each line) to your .cshrc file
## Set MATH alias - takes an arithmetic assignment statement                                                                  
## as argument, e.g., newvar = var1 + var2                                                                                    
## Separate all items and operators in the expression with blanks                                                             
#alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
if ( $IPAT_FACTOR > 0 ) then
    MATH dwell_s = $dwell_s / $IPAT_FACTOR
endif
echo; echo "*** dwell time (s) for $SCANNER and IPAT_FACTOR of $IPAT_FACTOR = $dwell_s\n"

set asym_s  = 0.00246; # fieldmap asym time in sec (default 0.00246 for Siemens?)



#---------------------------------------
# for skipping to the middle of the script if needed.
if ( 0 ) then
endif


#---------------------------------------
# 1. Prepare funcs, mag and phase images

# Create NIFTI versions of necessary files in fmap directory
# N.B. If the data is oblique, then conversion to NIFTI (and consequently back to BRIK in the final step)
#      will alter the centers/origin information (due to different formats for storing this info). The work
#      around, implemented in other scripts, is to use 3drefit to copy certain params from undistroted files.
##3dAFNItoNIFTI -float -prefix phase ../${SUBJ}_phz_shft+orig.BRIK; # need phase a floating-point precision, so must use 3dAFNItoNIFTI, not 3dcopy
cp ../${SUBJ}_phz.nii.gz phase.nii.gz
##3dAFNItoNIFTI -prefix magnitude ../${SUBJ}_mag_shft+orig.BRIK; # 3dcopy ../${SUBJ}_mag_shft+orig.BRIK magnitude.nii.gz
cp ../${SUBJ}_mag.nii.gz magnitude.nii.gz
3dTstat -prefix ./reference_meanTPs ../${SUBJ}_${REFRUN}${curfunc}+orig.BRIK; # refrence run averaged over timepoints
3dAFNItoNIFTI -prefix reference reference_meanTPs+orig; # 3dcopy ../${SUBJ}_${REFRUN}${curfunc}+orig.BRIK reference.nii.gz


# it's possible that you need to adjust the orientation of the mag and phase to match the EPI...
fslswapdim magnitude.nii.gz z x y magnitude.nii.gz
fslswapdim phase.nii.gz     z x y phase.nii.gz

if ( $interactive ) then
    # probably a good idea to check that the orientation of the phase, mag and reference image match
    echo ""
    echo "*******************************************"
    echo "Check orientation of phase, mag and functional...close fslview to continue"
    fslview magnitude.nii.gz phase.nii.gz &
    fslview reference.nii.gz
endif

# TODO: clean up old files that we'll never use, like phase.nii.gz

#-------------------------------------------------------------------------------------
# 1a. EPI: Skull strip your reference EPI using bet, check results and edit in FSLVIEW

# create binary mask, to be edited in FSLVIEW
# consider the following options (see bet help for more details)
#   -f 0.3 (default is 0.5, lower gets more brain)
#   -c x y z (center of mass in voxel space - get from fslview of reference.nii.gz and click to center of brain)
#   -r radius (approximate radius of brain in mm, count slices to edge of brain from center and multiply by voxel resolution)
#   -Z for z-plane padding (could be useful for partial coverage, but significantly increases processing time)
bet reference.nii.gz reference -n -m -f 0.3
if ( $interactive ) then
    # check mask and edit (interactive mode - will continue with script after exiting fslview)
    echo ""
    echo "*******************************************"
    echo "Check mask and edit...close fslview to continue"
    fslview reference.nii.gz reference_mask.nii.gz
endif
# apply mask to epi for skull stripping
fslmaths reference.nii.gz -mas reference_mask.nii.gz reference_brain.nii.gz


#---------------------------------------------------------------------------------
# 2. Create a brain mask from mag image, edit it, use it to create mag_brain image
# create mask
# you can try these optional arguments, -c/-f/-r, to try to get optimal results
bet magnitude.nii.gz magnitude -n -m -f 0.6
if ( $interactive ) then
    # check mask and edit (interactive mode - will continue with script after exiting fslview)
    echo ""
    echo "*******************************************"
    echo "Check mask and edit...close fslview to continue"
    fslview magnitude.nii.gz magnitude_mask.nii.gz
endif
# apply mask
fslmaths magnitude.nii.gz -mas magnitude_mask.nii.gz magnitude_brain


#---------------------------------------
# 3a. Convert Phase units to Radians/sec
#
# This next step is specific to the Siemens stock gre fieldmap sequence.
# It will be different with fieldmap scans at different scanners, depending on how 
# the phase information is encoded in the image by the fieldmap sequence.
# The Siemens stock gre fieldmap sequence provides you with a "single real fieldmap image"
# and its values may range from 0 to +4094 or from -4096 to +4092.
# The range depends on what tool you used to convert from to NIFTI.
# You must convert to rad/s (-1277 to +1277 range).
# It appears that dcm2nii scales the phase image to -4096 to +4092.

# guess which range we're working with based on min phase value
fslmaths phase.nii.gz -Xmin -Ymin -Zmin phase_min
3dmaskdump -noijk phase_min.nii.gz > phase_min.1D
set phase_min = `tail -1 phase_min.1D`;
if ( $phase_min < 0 ) then
    # assume range is -4069 to +4092
    fslmaths phase.nii.gz -mul 3.14159 -div 4096 -div $asym_s phase_rads.nii.gz
else
    # assume range is 0 to +4094
    fslmaths phase.nii.gz -sub 2047.5 -mul 3.14159 -div 2047.5 -div $asym_s phase_rads.nii.gz
endif

# 3b. Smooth (aka regularise) the phase image to reduce noise
# Smooth the phase image in 3D by $sm mm 
#   play with that value. I've seen 0.5 - 4mm
#   It doesn't effect the unwarping except it helps clean up around the brain edges a bit. 
# Note that while $UNWARPDIR is not used in this step, it is part of the file name so that
#   any script outside of the current script must be explicit in defining the smoothing and
#   unwarping directions to read the correct files.
fugue --loadfmap=phase_rads.nii.gz -s $FMAP_SM --savefmap=phase_rads_s${FMAP_SM}$UNWARPDIR


#--------------------------------------------------------------------------------------
# 4. Forward warp (aka distort) the brain-only mag image using the smoothed phase image.
# should you specify the mag brain mask here?  MGH script uses it, but comments that they are not sure what it does.
# MP tried both and doesn't see a difference, so we leave it out.
#
# Check that forward warped mag image looks to have similar distortion as your reference
# (i.e. check your unwarp direction)
fugue \
    --verbose \
    --in=magnitude_brain.nii.gz \
    --unwarpdir=$UNWARPDIR \
    --dwell=$dwell_s \
    --asym=$asym_s \
    --loadfmap=phase_rads_s${FMAP_SM}${UNWARPDIR}.nii.gz \
    --nokspace \
    --warp=magnitude_brain_warped_s${FMAP_SM}$UNWARPDIR


#---------------------------------------------------------
# 5. Register forward warped magnitude to reference_brain, check it.
# Register the magnitude_brain_warped to the reference_brain (a
#  skull-stripped, motion corrected, reference EPI volume)
# We do this so we can get the transform matrix to get the phase into register with data.
# Options include -nosearch to avoid large changes that can lead to flipping of directions
# (not sure if -paddingsize is used without -applyxfm)
flirt \
    -paddingsize 4 \
    -in magnitude_brain_warped_s${FMAP_SM}${UNWARPDIR}.nii.gz \
    -dof 6 \
    -ref reference_brain.nii.gz \
    -out magnitude_brain_warped_s${FMAP_SM}${UNWARPDIR}_to_reference_brain \
    -omat fieldmap2reference_brain.mat

# also align the magnitue_brain (true anatomy) to the the same reference brain
# so that we can compare how well the forwarded warped distortion match the "raw" epi
# if results are way off (and mag seems cut in half), might try to add -nosearch
flirt \
    -paddingsize 4 \
    -in magnitude_brain.nii.gz \
    -dof 6 \
    -ref reference_brain.nii.gz \
    -out magnitude_brain_to_reference_brain

if ( $interactive ) then
    # check alignment (interactive mode - will continue with script after exiting fslview)
    echo ""
    echo "*******************************************"
    echo "Check alignment...close fslview to continue"
    fslview magnitude_brain_warped_s${FMAP_SM}${UNWARPDIR}_to_reference_brain magnitude_brain_to_reference_brain reference_brain.nii.gz
endif


#----------------------------------------------------------------------------------------
# 6. Use above transformation matrix to align the smoothed phase image to reference_brain
# FSL
flirt \
    -paddingsize 4 \
    -in phase_rads_s${FMAP_SM}${UNWARPDIR}.nii.gz \
    -ref reference_brain.nii.gz \
    -applyxfm \
    -init fieldmap2reference_brain.mat \
    -out phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain

#----------------------------------------------------------------------------------------------
# 7- Unwarp the skull-stripped reference image (reference_brain) and check the result.
#    Also check the non skull-stripped epi too

# reference_brain (skull-stripped)
fugue \
    --verbose \
    --in=reference_brain.nii.gz \
    --icorr \
    --unwarpdir=$UNWARPDIR \
    --dwell=$dwell_s \
    --asym=$asym_s \
    --loadfmap=phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz \
    --unwarp=reference_brain_unwarped_s${FMAP_SM}${UNWARPDIR}

# reference (non skull-stripped)
fugue \
    --verbose \
    --in=reference.nii.gz \
    --icorr \
    --unwarpdir=$UNWARPDIR \
    --dwell=$dwell_s \
    --asym=$asym_s \
    --loadfmap=phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz \
    --unwarp=reference_unwarped_s${FMAP_SM}${UNWARPDIR}


# Now you can visually see which EPI, original or unwarped, provides a better match to the Anatomy (mag)
# using the same alignment method (anat->epi) that is standard in our fMRI preprocessing stream...
# ...and decide whether all this is worth it for you.
#
# Check results by registering:
# 1. magnitude_brain (unwarped skull-stripped mag = true anatomy) -> reference_brain_unwarped_s?? (post-undistortion, output of fugue command)
flirt \
    -paddingsize 4 \
    -in magnitude_brain.nii.gz \
    -dof 6 \
    -ref reference_brain_unwarped_s${FMAP_SM}${UNWARPDIR}.nii.gz \
    -out magnitude_brain_to_reference_brain_unwarped_s${FMAP_SM}${UNWARPDIR}

# 2. magnitude_brain (unwarped skull-stripped mag = true anatomy) -> reference_brain (pre-undistortion, "raw epi")
# N.B. we already did this registration above (when checking the forward warping of the mag)

# and visualize results...
if ( $interactive ) then
    # (interactive mode - will continue with script after exiting fslview)
    echo ""
    echo "*******************************************"
    echo "Check results...close fslview to continue"
    fslview magnitude_brain_to_reference_brain_unwarped_s${FMAP_SM}${UNWARPDIR} reference_brain_unwarped_s${FMAP_SM}${UNWARPDIR} magnitude_brain_to_reference_brain reference_brain
endif


#----------------------------
# 8. Unwarp all of your data.
#  i.e., your motion corrected EPI time series.

# N.B. this could be done in another script that checks to the presence of phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz

#foreach run ($funcruns)
#    echo; echo "${run}"
#
#    # convert copy to NIFTI
#    3dcopy ../${SUBJ}_${run}${curfunc}+orig.BRIK ${SUBJ}_${run}${curfunc}.nii.gz
#    
#    # unwarp
#    fugue \
#	--verbose \
#	--in=${SUBJ}_${run}${curfunc}.nii.gz \
#	--icorr \
#	--unwarpdir=$UNWARPDIR \
#	--dwell=$dwell_s \
#	--asym=$asym_s \
#	--loadfmap=phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz \
#	--unwarp=${SUBJ}_${run}${curfunc}ud.nii.gz
#    
#    CheckOutputFile.sh ${SUBJ}_${run}${curfunc}ud.nii; if ($status) exit 1
#
#    # copy back out of tmp directory
#    3dcopy ${SUBJ}_${run}${curfunc}ud.nii.gz ../${SUBJ}_${run}${curfunc}ud; # default is +orig.BRIK/HEAD
#end

#------------------------------------------------------------
# DONE !!!!!!!!!!!!!!!!!!!!
# Now you should only have to use 6 DOF when registering your T1 structural scan to you DTI or EPI data.
#------------------------------------------------------------

echo; echo "== $0 complete"
exit 0
