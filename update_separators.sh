#!/bin/bash
# update_separators.sh - RETIRED, DO NOT RUN
#
# This script's global sed replacements (e.g. 's/\* *100/* 80/g') were intended
# for separator lines like '='*100, but they also corrupted physics constants in
# code/single_q_analysis.py on 2026-06-08:
#   - kT_THz = kT_cm * const.c * 100 / 1e12   became   * 80   (T wrong by 20%)
#   - longitudinal_signed = 100 * Q_dot_e     became   80 *   (wrong % scale)
# Those constants were restored on 2026-06-09. Running this script again would
# re-corrupt them, so it now refuses to run.
#
# If separators ever need restandardizing, do it with a pattern anchored to the
# quoted separator character, e.g.:  sed "s/'='\* *100/'='* 80/g"

echo "ERROR: update_separators.sh is retired - it corrupts numeric constants."
echo "See comments in this file and CODE_REVIEW.md for details."
exit 1
