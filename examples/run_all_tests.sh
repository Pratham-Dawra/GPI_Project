#!/bin/bash

# list="1 2 3 4 5 6 7 8 3_fy 8_fy"
list="3 4 5 6 7 8 3_fy 8_fy"

for l in ${list} ; do
   ./run_single_test.sh ${l}
done
