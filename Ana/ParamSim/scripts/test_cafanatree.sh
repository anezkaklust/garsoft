#!/bin/bash

#bin/mcp_skimmer --infile anatree_100.root --outfile skimmed.root --debug 0
#echo "done with running skimming module"
cafanatree_module --infile andytree.root --outfile caf.root --correct4origin 1

echo "done with running cafanatree module"
