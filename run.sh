rm *.o
rm PatternMatch
make
rm *.o

./PatternMatch -name TEST -input ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf -ploidy 2 -ldthreshold 1.0