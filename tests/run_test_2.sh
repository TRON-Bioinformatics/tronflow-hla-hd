#!/bin/bash

output_folder=output/test2

nextflow main.nf -profile test,conda --output $output_folder --input_bams test_data/test_input_bams.txt

test -s $output_folder/TESTX_S1_L001/TESTX_S1_L001_final.result.txt || { echo "Missing test 1 output file!"; exit 1; }
test -s $output_folder/TESTX_S1_L002/TESTX_S1_L002_final.result.txt || { echo "Missing test 2 output file!"; exit 1; }
