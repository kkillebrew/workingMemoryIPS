#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/analysis/runwise_glm
#
# This script will execute a volume-based GLM, separately for each run
#
# N.B. 3dREMLfit is the command to run a GLM accounting for temporal auto-correlation
#      in the data.  To use it, we first run the command we would use with 3dDeconvolve,
#      but we tell 3dDeconvolve to stop (just do the prep work and write out some files):
#          -x1D_stop
#      3dDeconvolve automatically creates a [prefix].REML_cmd txt file that you can
#      execute to run the appropriate 3dREMLfit command.
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../../scripts/global_vars.sh


# variable validation
if ( $#VOLSMOOTHS == 0) then
    echo; echo "ERROR: you must define at least one level of smoothing in global_vars.sh, even if you set it to zero (0)."
    exit 1
endif


# local variables
set overwrite = 0; # overwrite bucket output files if they already exist?  if not, just stop script


set sm = $VOLSMOOTHS[1]; # we only want to worry about 1 smoothing level for this quality check
set curfunc = "_tscvrsm${sm}_norm"; # fully specifies all preprocessing (e.g., volume or surface smoothing)


set i = 0;
foreach run ($ALLRUNS)
    @ i = $i + 1;

    # setup tmp stim_time files with just the desired run extracted
    foreach reg ( face house genobject scram )
	set tmp = `head -$i ../../scripts/stim_times/${reg}_times.1D | tail -1`
	if ( "$tmp" == "*" ) then
	    # the first run in a stim_times file must be ** (not *) to indicate no events for that run
	    set tmp = "**";
	endif
	echo "$tmp" > $run.${reg}_times.1D
    end

    # now run the GLM using only this run

    # get a list of all preprocessed files that will serve as input to the GLM
    set glminput = "../${SUBJ}_${run}${curfunc}+orig"

    # verify output bucket doesn't exist.  if it does, and local variable overwrite is 1, delete it.
    set output_prefix = ${SUBJ}_${run}${curfunc}_bucket; # do NOT include +orig
    set output = ${output_prefix}+orig
    if ( -e $output.HEAD ) then
	echo "== output file $output.HEAD exists"
	if ( $overwrite ) then
	    echo "   overwriting..."
	    rm -f $output.HEAD $output.BRIK*
	else
	    echo "   stopping (either delete $output.HEAD/BRIK or set overwrite to 1 in $0)"
	    exit 1
	endif
    endif

    
    # run the GLM
    # N.B. if you are only running one GLM at a time on a computer, you can
    #      speed things up taking advantage of parallel processing by adding
    #      -jobs 4
    #      This tells 3dDeconvolve to use 4 processors.
    # N.B. could add a mask using -mask ${SUBJ}_automask+orig, but only makes sense for volume data
    3dDeconvolve \
	-input $glminput \
	-CENSORTR '*:0-5' \
	-polort 2 \
	-mask ${SUBJ}_automask+orig \
	-num_stimts 10 \
	-basis_normall 1 \
	-local_times \
	-stim_file 1 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[1]' \
	-stim_label 1 roll \
	-stim_base 1 \
	-stim_file 2 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[2]' \
	-stim_label 2 pitch \
	-stim_base 2 \
	-stim_file 3 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[3]' \
	-stim_label 3 yaw \
	-stim_base 3 \
	-stim_file 4 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[4]' \
	-stim_label 4 dIS \
	-stim_base 4 \
	-stim_file 5 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[5]' \
	-stim_label 5 dRL \
	-stim_base 5 \
	-stim_file 6 ../../preprocessing/${SUBJ}_${run}_tscvr.1D'[6]' \
	-stim_label 6 dAP \
	-stim_base 6 \
	-stim_times 7 ./$run.face_times.1D 'BLOCK(15,1)' \
	-stim_label 7 face \
	-stim_times 8 ./$run.house_times.1D 'BLOCK(15,1)' \
	-stim_label 8 house \
	-stim_times 9 ./$run.genobject_times.1D 'BLOCK(15,1)' \
	-stim_label 9 genobject \
	-stim_times 10 ./$run.scram_times.1D 'BLOCK(15,1)' \
	-stim_label 10 scram \
	-gltsym 'SYM: 3*face -house -scram -genobject' \
	-glt_label 1 face_vs_all \
	-gltsym 'SYM: 3*house -face -scram -genobject' \
	-glt_label 2 house_vs_all \
	-gltsym 'SYM: 3*genobject -face -scram -house' \
	-glt_label 3 genobject_vs_all \
	-gltsym 'SYM: face  -house' \
	-glt_label 4 face_vs_house \
	-gltsym 'SYM: genobject  -scram' \
	-glt_label 5 genobject_vs_scram \
	-bucket $output \
	-x1D ${output_prefix}.xmat.1D \
	-tout \
	-jobs 4
end


# get a copy of the NoSkull_SurfVol_Alnd_Exp
3dcopy ../${SUBJ}_NoSkull_SurfVol_Alnd_Exp+orig ./


echo; echo "== $0 complete"
exit 0


