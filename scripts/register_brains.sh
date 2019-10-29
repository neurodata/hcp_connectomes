PARTICIPANT_FILE := /mnt/ssd3/j1c/hcp1200/HCP_1200_participants.txt
SCRIPT_FILE := /mnt/ssd3/j1c/hcp_connectomes/scripts/run_reg.py

INPUT_PATH := /mnt/ssd3/j1c/hcp1200
OUTPUT_PATH := /mnt/ssd3/j1c/hcp1200_brain_volume/

cat ${PARTICIPANT_FILE} | parallel --jobs 52 --max-args=1 python ${SCRIPT_FILE} ${INPUT_PATH} ${OUTPUT_PATH} {1}