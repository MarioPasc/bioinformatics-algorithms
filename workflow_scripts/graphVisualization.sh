#!/bin/bash
TESTTYPE="readsNonEulerian"

./build/main testFleury "$TESTTYPE" > workflow_scripts/temp_output.txt
python ./Assembly/De_Brujin_Graphs/visualize_graph.py < workflow_scripts/temp_output.txt
