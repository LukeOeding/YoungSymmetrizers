#!/bin/bash
for i in  `seq 1 2622`
do
`sed -i "s/place:=1:/place:=$i:/" coreS4P3.d11.$i` &
done
