#!/bin/bash

#Select which runcards should be run, from runcards/runcard_**.yaml (must be 3 characters, include 00 if 001, 002, ...)
random_seed_end_list=(001)

#Select how many parallel sessions should run simultaneously (one screen session for each)
N_parallel=5


# Loop to start multiple screen sessions
for i in $(seq 1 $N_parallel); do
    echo "Starting screen session $i"
    
    # Run the command in a detached screen session
    screen -dm bash -c "
        for random_seed_end in ${random_seed_end_list[@]}; do  
            echo \"Running runcard_\$random_seed_end.yaml -R ${i}\$random_seed_end\"
            ../../bin/Sherpa -f runcards/runcard_\$random_seed_end.yaml -R ${i}\$random_seed_end
        done
    "
done
