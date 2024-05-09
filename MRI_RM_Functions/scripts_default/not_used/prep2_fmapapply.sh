#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/preprocessing/
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# local variables
# dwelll time (echo spacing) in s
switch ($scanner)
    case allegra:
	set dwell_s  = 0.00069; # dwelll time (echo spacing)in sec
	breaksw
    case skyra:
	set dwell_s  = 0.00075; # dwelll time (echo spacing)in sec
	breaksw
    case UCDavis_skyra:
	set dwell_s  = 0.00075; # dwelll time (echo spacing)in sec (*** NEED TO CHECK THIS ***)
	breaksw
    default:
	echo "ERROR: scanner must be 'skyra', 'allegra', or 'UCDavis_Skyra'"
	exit 1
endsw
# N.B. the "effective" echo spacing must take into account any acceleration factor
#      dwell_time_sec = (EPI_echo_spacing_ms / iPAT_factor) / 1000
# N.B. to use the MATH command, add the following lines (excluding first # of each line) to your .cshrc file
## Set MATH alias - takes an arithmetic assignment statement                                                                  
## as argument, e.g., newvar = var1 + var2                                                                                    
## Separate all items and operators in the expression with blanks                                                             
#alias MATH 'set \!:1 = `echo "\!:3-$" | bc -l`'
if ( $iPAT_factor > 0 ) then
    MATH dwell_s = $dwell_s / $iPAT_factor
endif
echo; echo "== dwell time (s) for $scanner and iPAT_factor of $iPAT_factor = $dwell_s"; echo

set asym_s  = 0.00246; # fieldmap asym time in sec (default 0.00246 for Siemens?)


# N.B. the $UNWARPDIR and $FMAP_SM parameters, defined in {experiment}_vars.sh will determine the final epi-aligned phase map's filename
#     e.g., phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz


# for quick commenting
if ( 0 ) then
endif
set curfunc = "_tscvr"; # "" at begining of script


#______________________________________________________________________________
# get a copy of the original data.
# by keeping the original data in ./orig/, you can reset the preprocessing/analysis
# by deleting all files in ./analysis/, without fear of loosing the original data or
# having to re-convert from dicom.
cp ../orig/${SUBJ}_mag.nii.gz ../orig/${SUBJ}_phz.nii.gz ./


#____________________________________
# apply undistortion
set i = 0;
foreach run ($ALLRUNS)
    echo; echo "== processing run $run"

    # check that the aligned_phase file exists
    # the aligned_phase file should be created manually using the fmap_undistort.sh script
    # N.B. $UNWARPDIR and $FMAP_SM are defined in {task}_vars.sh
    set aligned_phase = ./fmap/phase_rads_s${FMAP_SM}${UNWARPDIR}_to_reference_brain.nii.gz
    if ( ! -e $aligned_phase ) then
	echo "ERROR: can't find aligned_phase file $aligned_phase.  You need to go to ./fmap/ (create if necessary) and run fmap_undistort.sh before proceeding..."
	exit 1
    endif

    # convert copy to NIFTI
    3dcopy ${SUBJ}_${run}${curfunc}+orig.BRIK ${SUBJ}_${run}${curfunc}.nii.gz
   
    # unwarp
    fugue \
	--verbose \
	--in=${SUBJ}_${run}${curfunc}.nii.gz \
	--icorr \
	--unwarpdir=$UNWARPDIR \
	--dwell=$dwell_s \
	--asym=$asym_s \
	--loadfmap=$aligned_phase \
	--unwarp=${SUBJ}_${run}${curfunc}ud.nii.gz
	    
    # copy back to BRIK, remove tmp NIFTI files (orig and unwarped)
    3dcopy ${SUBJ}_${run}${curfunc}ud.nii.gz ${SUBJ}_${run}${curfunc}ud; # default is +orig.BRIK/HEAD
    rm -f ${SUBJ}_${run}${curfunc}*.nii.gz

    # in case the data was oblique, re-apply the origin/oblique information from original data
    # (this information is lost during BRIK->NIFTI conversion necessary for unwarping)
    3dAttribute IJK_TO_DICOM_REAL ${SUBJ}_${run}${curfunc}+orig > epi_obl_matrix.1D
    3drefit -atrfloat IJK_TO_DICOM_REAL epi_obl_matrix.1D ${SUBJ}_${run}${curfunc}ud+orig
    3drefit -duporigin ${SUBJ}_${run}${curfunc}+orig ${SUBJ}_${run}${curfunc}ud+orig
end

# update curfunc
set curfunc = ${curfunc}ud

echo; echo "== $0 complete"
exit 0
