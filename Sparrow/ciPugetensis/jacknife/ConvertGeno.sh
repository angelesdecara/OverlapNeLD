#!/bin/bash
for i in {0..101};do
    ~/Downloads/jdk1.7.0_80/bin/java -jar ~/Downloads/PGDSpider_2.1.1.5/PGDSpider2-cli.jar -inputformat VCF -outputformat GENEPOP -inputfile pop${i}.vcf -outputfile pop${i}.geno -spid ../../ciNuttalli/template_VCF_GENEPOP.spid
done
