#!/bin/bash
for i in `seq 1 163`
 do
 `nice nohup maple <input.6.$i &> 10quad_out.6.$i &`
 done
