#!/bin/bash

for i in $(ls 1a93_B.p22.mincrd)
do
	case=$(echo $i | cut -d \/ -f 2 | sed 's/.mincrd//g')
	./Run.testCase.min $case.mincrd $case.parm $case
done
