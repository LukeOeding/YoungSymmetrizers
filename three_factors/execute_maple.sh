#!/bin/bash
for i in 2553 2607
 do
 bsub -q long -o out_sec4.11.$i "/tools/maple18/bin/maple <coreS4P3.d11.$i"
 done
