#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ./scripts/global_vars.sh


# for quick commenting
if (0) then
endif


# make directories for data extraction and subsequent preprocessing and analysis
mkdir orig preprocessing analysis


# setup a list of scan directories that will be tar-gzipped
set toarchive = ""


#-----------------------------
# functional scans
set i = 0;
foreach scan ($FUNC_SCANS)
    @ i = $i + 1;

    # set slice acquisition depending on number of slices (for Siemens scanners)
    set odd_slices = `expr $NSLICES[$i] % 2`;
    if ($odd_slices) then
	set sacq = alt+z;
    else
	set sacq = alt+z2;
    endif

    # get a two-digit run id
    if ( $i < 10) then
	set ii = 0$i
    else
	set ii = $i
    endif
    set run = r$ii

    # what is the directory for the dicoms?
    set this_dicomdir = ${scan}_$FUNC_SERIES_NAME[$i]
    set toarchive = "$toarchive $this_dicomdir"

    # convert from dicom to BRIK
    echo; echo "== run $run, scan $scan"    
    to3d \
	-oblique_origin \
	-time:zt $NSLICES[$i] $TPS[$i] $TR[$i] $sacq \
        -save_outliers ${SUBJ}_${run}_outliers.1D \
        -prefix ${SUBJ}_${run} \
        ${this_dicomdir}/*.dcm
    
    # move BRIK/HEAD and 1D (outlier) output
    mv ${SUBJ}_${run}* ./orig/
end


#-----------------------------
# fieldmaps
if ( ${%FMAP_SCANS} != 0 ) then
    set magscan = $FMAP_SCANS[1]
    echo; echo "== fieldmap magnitude, scan $magscan"    
    dcm2nii ${magscan}_gre_field_mapping/*.dcm
    mv ${magscan}_gre_field_mapping/*.nii.gz ./orig/${SUBJ}_mag.nii.gz

    set phzscan = $FMAP_SCANS[2]
    echo; echo "== fieldmap phase, scan $phzscan"
    dcm2nii ${phzscan}_gre_field_mapping/*.dcm
    mv ${phzscan}_gre_field_mapping/*.nii.gz ./orig/${SUBJ}_phz.nii.gz

    # add field map dicom dirs to archive list
    set toarchive = "$toarchive ${magscan}_gre_field_mapping ${phzscan}_gre_field_mapping"    
endif


#-----------------------------
# anatomical data
# get a copy of the NoSkull_SurfVol for this subject from the surfaces directory
3dcopy $SUMADIR/${SUBJ}_NoSkull_SurfVol+orig ./orig/${SUBJ}_surfvol+orig


#-----------------------------
# archive dicoms
echo; echo "== tar-gzipping dicom directories..."
tar -czf ${SUBJ}_dicom.tgz $toarchive
echo "done"

# remove raw dicom directories
rm -fR $toarchive


echo; echo "== $0 complete"
exit 0
