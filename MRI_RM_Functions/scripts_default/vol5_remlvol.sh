#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/$SUBJ/analysis/
#
# This script will execute a volume-based REMLfit (GLM accounting for temporal auto-correlation)
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
set overwrite = 1; # overwrite bucket output files if they already exist?  if not, just stop script
set use_trialwise_stimtimes = 0; # boolean, should we use the trial-wise stimtimes, which define each trial, instead of each block?
if ( $use_trialwise_stimtimes ) then
    set stimtime_suffix = "_TRIALS";
    set stimtime_hrf    = 'BLOCK(2,1)'; # HRF for duration of a single trial
else
    set stimtime_suffix = "";
    set stimtime_hrf    = 'BLOCK(36,1)'; # HRF for duration of an entire block
endif



foreach sm ( $VOLSMOOTHS )
    echo; echo "== starting GLM for smoothing of $sm mm"

    set curfunc = "_tscvrsm${sm}_norm"; # fully specifies all preprocessing (e.g., volume or surface smoothing)

    # get a list of all preprocessed files that will serve as input to the GLM
    set glminput = ""
    foreach run ($ALLRUNS)
	set glminput = "$glminput ${SUBJ}_${run}${curfunc}+orig"
    end

    # verify output bucket doesn't exist.  if it does, and local variable overwrite is 1, delete it.
    set output_prefix = ${SUBJ}${curfunc}${stimtime_suffix}_bucket; # do NOT include the _REML suffix that 3dREMLfit automatically adds (or +orig)
	if ( -e ${output_prefix}_REML+orig.HEAD ) then
	    echo "== output file ${output_prefix}_REML+orig.HEAD exists"
	    if ( $overwrite ) then
		echo "   overwriting..."
		rm -f ${output_prefix}_REML*+orig.HEAD ${output_prefix}_REML*+orig.BRIK*
	    else
		echo "   stopping (either delete ${output_prefix}_REML+orig.HEAD/BRIK or set overwrite to 1 in $0)"
		exit 1
	    endif
	endif

    
    # setup the GLM
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
    -stim_times 8 ../scripts/stim_times/${SUBJ}_let_left${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 8 let_left \
    -stim_times 9 ../scripts/stim_times/${SUBJ}_let_right_cue.1D 'BLOCK(2,1)' \
    -stim_label 9 let_right_cue \
    -stim_times 10 ../scripts/stim_times/${SUBJ}_let_right${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 10 let_right \
    -stim_times 11 ../scripts/stim_times/${SUBJ}_let_both_cue.1D 'BLOCK(2,1)' \
    -stim_label 11 let_both_cue \
    -stim_times 12 ../scripts/stim_times/${SUBJ}_let_both${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 12 let_both \
    -stim_times 13 ../scripts/stim_times/${SUBJ}_let_pass_cue.1D 'BLOCK(2,1)' \
    -stim_label 13 let_pass_cue \
    -stim_times 14 ../scripts/stim_times/${SUBJ}_let_pass${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 14 let_pass \
    -stim_times 15 ../scripts/stim_times/${SUBJ}_ori_left_cue.1D 'BLOCK(2,1)' \
    -stim_label 15 ori_left_cue \
    -stim_times 16 ../scripts/stim_times/${SUBJ}_ori_left${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 16 ori_left \
    -stim_times 17 ../scripts/stim_times/${SUBJ}_ori_right_cue.1D 'BLOCK(2,1)' \
    -stim_label 17 ori_right_cue \
    -stim_times 18 ../scripts/stim_times/${SUBJ}_ori_right${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 18 ori_right \
    -stim_times 19 ../scripts/stim_times/${SUBJ}_ori_both_cue.1D 'BLOCK(2,1)' \
    -stim_label 19 ori_both_cue \
    -stim_times 20 ../scripts/stim_times/${SUBJ}_ori_both${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 20 ori_both \
    -stim_times 21 ../scripts/stim_times/${SUBJ}_ori_pass_cue.1D 'BLOCK(2,1)' \
    -stim_label 21 ori_pass_cue \
    -stim_times 22 ../scripts/stim_times/${SUBJ}_ori_pass${stimtime_suffix}.1D $stimtime_hrf \
    -stim_label 22 ori_pass \
    -num_glt 18 \
    -glt_label 1 let_vs_ori \
    -gltsym 'SYM: let_left -let_right' \
    -gltsym 'SYM: let_left ori_left -let_right -ori_right' \
    -glt_label 2 left_vs_right \
    -gltsym 'SYM: let_left let_right let_both -ori_left -ori_right -ori_both' \
    -glt_label 3 let_left_vs_right \
    -gltsym 'SYM: ori_left -ori_right' \
    -glt_label 4 ori_left_vs_right \
    -gltsym 'SYM: let_both ori_both -let_pass -ori_pass' \
    -glt_label 5 both_vs_pass \
    -gltsym 'SYM: let_both -let_pass' \
    -glt_label 6 let_both_vs_pass \
    -gltsym 'SYM: ori_both -ori_pass' \
    -glt_label 7 ori_both_vs_pass \
    -gltsym 'SYM: let_left ori_left -let_pass -ori_pass' \
    -glt_label 8 left_vs_pass \
    -gltsym 'SYM: let_left -let_pass' \
    -glt_label 9 let_left_vs_pass \
    -gltsym 'SYM: let_right -let_pass' \
    -glt_label 10 let_right_vs_pass \
    -gltsym 'SYM: let_right ori_right -let_pass -ori_pass' \
    -glt_label 11 right_vs_pass \
    -gltsym 'SYM: ori_left -ori_pass' \
    -glt_label 12 ori_left_vs_pass \
    -gltsym 'SYM: ori_right -ori_pass' \
    -glt_label 13 ori_right_vs_pass \
    -gltsym 'SYM: let_left let_right let_both ori_left ori_right ori_both' \
    -glt_label 14 allwm \
    -gltsym 'SYM: let_left let_right let_both ori_left ori_right ori_both -3*let_pass -3*ori_pass' \
    -glt_label 15 allwm_vs_pass \
    -gltsym 'SYM: let_left ori_left' \
    -glt_label 16 left \
    -gltsym 'SYM: let_right ori_right' \
    -glt_label 17 right \
    -gltsym 'SYM: let_both ori_both' \
    -glt_label 18 both \
	-bucket $output_prefix \
	-x1D $output_prefix.xmat.1D \
	-x1D_stop \
	-tout
	
    # run 3dREMLfit command that was produced by 3dDeconvolve
    # N.B. 3dREMLfit will automatically run in parallel across all of the processors
    #      on the current computer.  If this is undesirable, you should uncomment the
    #      following line:
    #      setenv OMP_NUM_THREADS 1; # uncomment this line to NOT use multi-thread processing for 3dREMLfit
    tcsh $output_prefix.REML_cmd
end
    
echo; echo "== $0 complete"
exit 0

