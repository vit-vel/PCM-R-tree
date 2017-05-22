#!/usr/bin/env bash

HOME_PROJECT_DIR=".."

# args:
# $1 --- max number of elements inserted to r-tree

# clear files
> orig_results.csv
> pcm_results.csv

for ((i=10; i <= $1 ; i *= 10))
do
    valgrind --tool=exp-dhat --show-top-n=100 $HOME_PROJECT_DIR/run-rtree $i 0 &> tmp_orig_results.dhat
    # array[2] with bytes read and bytes written
    orig_stats=($($HOME_PROJECT_DIR/utils/parse_dhat_result.py tmp_orig_results.dhat))

    valgrind --tool=exp-dhat --show-top-n=100 $HOME_PROJECT_DIR/run-rtree $i 1 &> tmp_pcm_results.dhat
    pcm_stats=($($HOME_PROJECT_DIR/utils/parse_dhat_result.py tmp_pcm_results.dhat))


    # terminal output
    echo "===================================================="

    echo "original R-tree with $i elements"
    echo "bytes read: ${orig_stats[0]}"
    echo "bytes written: ${orig_stats[1]}"

    echo "----------------------------------------------------"

    echo "PCM R-tree with $i elements"
    echo "bytes read: ${pcm_stats[0]}"
    echo "bytes written: ${pcm_stats[1]}"

    echo "----------------------------------------------------"

    read_advantage=$(echo "${orig_stats[0]} / ${pcm_stats[0]}" |bc -l)
    write_advantage=$(echo "${orig_stats[1]} / ${pcm_stats[1]}" |bc -l)

    echo "PCM superiority over original"
    echo "read: $read_advantage"
    echo "write: $write_advantage"

    # save results in csv file
    echo "$i,${orig_stats[0]},${orig_stats[1]}" >> orig_results.csv
    echo "$i,${pcm_stats[0]},${pcm_stats[1]}" >> pcm_results.csv
done
