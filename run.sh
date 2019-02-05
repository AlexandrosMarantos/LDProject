rm *.o
rm PatternMatch
make
rm *.o

./PatternMatch -name TEST -input vcf_from_ms.vcf -ploidy 2 -ldthreshold 1.0