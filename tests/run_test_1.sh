#!/bin/bash

output_folder=output/test1
nextflow main.nf -profile test --output $output_folder
test -s $output_folder/TESTX_S1_L001/TESTX_S1_L001_final.result.txt || { echo "Missing test 1 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002/TESTX_S1_L002_final.result.txt || { echo "Missing test 2 output file!"; exit 1; }
