#!/bin/bash

for i in $(ls 1tsrb.inpcrd)
do
	case=$(echo $i | cut -d \/ -f 2 | sed 's/.inpcrd//g')
	./Run.testCase.min $case.inpcrd $case.prmtop $case
done
