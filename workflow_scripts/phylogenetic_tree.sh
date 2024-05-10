#!/bin/bash

./build/testNJ > workflow_scripts/temp_output.txt
python ./Alignment/MultipleSequenceAlignment/phylogenetic_tree.py workflow_scripts/temp_output.txt
