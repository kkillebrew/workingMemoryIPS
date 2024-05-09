#!/bin/tcsh
# this script should be run from ~/MR_DATA/CLAB/[experiment]/
#
# This script is for doing a serial batch process on multiple subjects

echo "== $0 starting"

set subjs = "CB GC GG JV KK KM LS NS"; # MG

foreach s ( $subjs )
    echo "== batch iteration for $s"
    cd $s/analysis/
    ../../scripts_default/vol5_remlvol.sh
    cd ../../
end

echo; echo "== $0 complete"
exit 0

