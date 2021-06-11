#!/bin/bash
#pugetensis are 102, so jk with 101
#to make 10 pops w 101 indivs=total-1
bcfs=(`ls ../singlesample/*bcf`)
for p in {0..101};do
    newpop=()
    i=0
    for i in {0..101};do
	if [ $i != $p ]; then
	    newpop+=(${bcfs[$i]} )
	fi
    done
    echo ${newpop[@]}
    bcftools merge  ${newpop[@]} > pop${p}.vcf
done
