#!/bin/bash

#bin/mcp_skimmer --infile anatree_100.root --outfile skimmed.root --debug 0
#echo "done with running skimming module"

#### For CDR production sample, where the TPC is not in (0, 0, 0), correction for mcp position wrt TPC center
#cafanatree_module --infile andytree.root --outfile caf.root --correct4origin 1 --originTPC 0.00012207 -150.473 1486

#### For CDR production sample, where the TPC is not in (0, 0, 0), no correction for mcp position wrt TPC center
#cafanatree_module --infile andytree.root --outfile caf.root --correct4origin 0 --originTPC 0.00012207 -150.473 1486

#### For any sample, where the TPC is in (0, 0, 0)
#cafanatree_module --infile andytree.root --outfile caf.root --correct4origin 0 --originTPC 0 0 0

echo "done with running cafanatree module"
