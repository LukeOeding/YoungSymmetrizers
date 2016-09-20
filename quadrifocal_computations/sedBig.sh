#!/bin/bash
for i in `seq 1 163`
do
`sed -i "s/place:=1:/place:=$i:/" input.6.$i`&
done
