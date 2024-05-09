#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/analysis/
#
# This script will execute a volume-based GLM
#
# N.B. 3dREMLfit is the command to run a GLM accounting for temporal auto-correlation
#      in the data.  To use it, we first run the command we would use with 3dDeconvolve,
#      but we tell 3dDeconvolve to stop (just do the prep work and write out some files):
#          -x1D_stop
#      3dDeconvolve automatically creates a [prefix].REML_cmd txt file that you can
#      execute to run the appropriate 3dREMLfit command.
echo "== $0 starting"


# get global (i.e., "shared") variables for this subject (by convention, global variables are all UPPERCASE)
source ../scripts/global_vars.sh


# variable validation
if ( $#VOLSMOOTHS == 0) then
    echo; echo "ERROR: you must define at least one level of smoothing in global_vars.sh, even if you set it to zero (0)."
    exit 1
endif


# local variables
set overwrite = 0; # overwrite bucket output files if they already exist?  if not, just stop script


foreach sm ( $VOLSMOOTHS )
    echo; echo "== starting GLM for smoothing of $sm mm"

    set curfunc = "_tscvrsm${sm}_norm"; # fully specifies all preprocessing (e.g., volume or surface smoothing)

    # get a list of all preprocessed files that will serve as input to the GLM
    set glminput = ""
    foreach run ($ALLRUNS)
	set glminput = "$glminput ${SUBJ}_${run}${curfunc}+orig"
    end

    # verify output bucket doesn't exist.  if it does, and local variable overwrite is 1, delete it.
    set output_prefix = ${SUBJ}${curfunc}_bucket; # do NOT include +orig
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
	-num_stimts 22 \
	-basis_normall 1 \
	-local_times \
	-stim_file 1 ${SUBJ}_mcparams.1D'[1]' \
	-stim_label 1 roll \
	-stim_base 1 \
	-stim_file 2 ${SUBJ}_mcparams.1D'[2]' \
	-stim_label 2 pitch \
	-stim_base 2 \
	-stim_file 3 ${SUBJ}_mcparams.1D'[3]' \
	-stim_label 3 yaw \
	-stim_base 3 \
	-stim_file 4 ${SUBJ}_mcparams.1D'[4]' \
	-stim_label 4 dIS \
	-stim_base 4 \
	-stim_file 5 ${SUBJ}_mcparams.1D'[5]' \
	-stim_label 5 dRL \
	-stim_base 5 \
	-stim_file 6 ${SUBJ}_mcparams.1D'[6]' \
	-stim_label 6 dAP \
	-stim_base 6 \
    -stim_times 7 ../scripts/stim_times/${SUBJ}_let_left_cue.1D 'BLOCK(2,1)' \
    -stim_label 7 let_left_cue \
    -stim_times 8 ../scripts/stim_times/${SUBJ}_let_left.1D 'BLOCK(36,1)' \
    -stim_label 8 let_left \
    -stim_times 9 ../scripts/stim_times/${SUBJ}_let_right_cue.1D 'BLOCK(2,1)' \
    -stim_label 9 let_right_cue \
    -stim_times 10 ../scripts/stim_times/${SUBJ}_let_right.1D 'BLOCK(36,1)' \
    -stim_label 10 let_right \
    -stim_times 11 ../scripts/stim_times/${SUBJ}_let_both_cue.1D 'BLOCK(2,1)' \
    -stim_label 11 let_both_cue \
    -stim_times 12 ../scripts/stim_times/${SUBJ}_let_both.1D 'BLOCK(36,1)' \
    -stim_label 12 let_both \
    -stim_times 13 ../scripts/stim_times/${SUBJ}_let_pass_cue.1D 'BLOCK(2,1)' \
    -stim_label 13 let_pass_cue \
    -stim_times 14 ../scripts/stim_times/${SUBJ}_let_pass.1D 'BLOCK(36,1)' \
    -stim_label 14 let_pass \
    -stim_times 15 ../scripts/stim_times/${SUBJ}_ori_left_cue.1D 'BLOCK(2,1)' \
    -stim_label 15 ori_left_cue \
    -stim_times 16 ../scripts/stim_times/${SUBJ}_ori_left.1D 'BLOCK(36,1)' \
    -stim_label 16 ori_left \
    -stim_times 17 ../scripts/stim_times/${SUBJ}_ori_right_cue.1D 'BLOCK(2,1)' \
    -stim_label 17 ori_right_cue \
    -stim_times 18 ../scripts/stim_times/${SUBJ}_ori_right.1D 'BLOCK(36,1)' \
    -stim_label 18 ori_right \
    -stim_times 19 ../scripts/stim_times/${SUBJ}_ori_both_cue.1D 'BLOCK(2,1)' \
    -stim_label 19 ori_both_cue \
    -stim_times 20 ../scripts/stim_times/${SUBJ}_ori_both.1D 'BLOCK(36,1)' \
    -stim_label 20 ori_both \
    -stim_times 21 ../scripts/stim_times/${SUBJ}_ori_pass_cue.1D 'BLOCK(2,1)' \
    -stim_label 21 ori_pass_cue \
    -stim_times 22 ../scripts/stim_times/${SUBJ}_ori_pass.1D 'BLOCK(36,1)' \
    -stim_label 22 ori_pass \
    -num_glt 12 \
    -gltsym 'SYM: let_left ori_left -let_right -ori_right' \
    -glt_label 1 left_vs_right \
    -gltsym 'SYM: let_left let_right let_both -ori_left -ori_right -ori_both' \
    -glt_label 2 let_vs_ori \
    -gltsym 'SYM: let_left -let_right' \
    -glt_label 3 let_left_vs_right \
    -gltsym 'SYM: ori_left -ori_right' \
    -glt_label 4 ori_left_vs_right \
    -gltsym 'SYM: let_both ori_both -let_pass -ori_pass' \
    -glt_label 5 both_vs_pass \
    -gltsym 'SYM: let_both -let_pass' \
    -glt_label 6 let_both_vs_pass \
    -gltsym 'SYM: ori_both -ori_pass' \
    -glt_label 7 ori_both_vs_pass \
    -gltsym 'SYM: let_left let_right let_both ori_left ori_right ori_both' \
    -glt_label 8 allwm \
    -gltsym 'SYM: let_left let_right let_both ori_left ori_right ori_both -3*let_pass -3*ori_pass' \
    -glt_label 9 allwm_vs_pass \
    -gltsym 'SYM: let_left ori_left' \
    -glt_label 10 left \
    -gltsym 'SYM: let_right ori_right' \
    -glt_label 11 right \
    -gltsym 'SYM: let_both ori_both' \
    -glt_label 12 both \
	-bucket $output \
	-x1D ${output_prefix}.xmat.1D \
	-tout \
	-jobs 4
end

echo; echo "== $0 complete"
exit 0
