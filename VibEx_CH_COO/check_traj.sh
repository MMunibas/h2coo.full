#!/bin/bash

# Function to check termination in files
check_termination_in_files() {
    output_filename=$1
    total_filename=$2
    filenames=(out*.out)
    ncount=0
    ecount=0

    for filename in "${filenames[@]}"; do
        echo "$filename" >> "$total_filename"
        last_line=$(tail -n 6 "$filename" | head -n 1 | tr -d '[:space:]')

        if [[ "$last_line" != *"NORMALTERMINATION"* ]]; then
            ecount=$((ecount + 1))
            job_index=$(echo "$filename" | grep -oE '[0-9]+')
            num_step=$(tail -n 100 "$filename" | grep "DYNA>" | tail -n 1 | grep -oE '[0-9]+'|head -n 1)
            echo $num_step
            echo "sbatch sub-${job_index}.sh # crashed at ${num_step}" >> $output_filename
        else
            ncount=$((ncount + 1))
        fi
    done

    echo "The total number of out files: ${#filenames[@]}"
    echo "All output files have been written in $total_filename"
    echo "The total number of normally terminated files: $ncount"
    echo "The total number of crashed files: $ecount"
    echo "All crashed files have been written in $output_filename"
}

echo "This script can be run only after all your CHARMM jobs are finished."

check_termination_in_files "resub.sh" "totalfile.txt"
chmod +x "resub.sh"

