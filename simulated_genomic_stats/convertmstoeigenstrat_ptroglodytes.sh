path=~/$1
MSOUT=~/$1/$1.out
LENGTH=1000000
NSEQ=100
MUT=1.5e-8
REC=0.7e-8
GENTIME=25
OUTPREFIX=$1

echo "hello"
cat $MSOUT | ~/ms2vcf -length $LENGTH -ploidy 2 >$path/$OUTPREFIX.vcf
python ~/vcf2eigenstrat_RT.py -v $path/$OUTPREFIX.vcf -o $path/$OUTPREFIX

echo "seq_length $LENGTH" > $path/$OUTPREFIX.par
echo "n_seq $NSEQ" >> $path/$OUTPREFIX.par
echo "mut_rate $MUT" >> $path/$OUTPREFIX.par
echo "recomb_rate $REC" >> $path/$OUTPREFIX.par
echo "generation_time $GENTIME" >> $path/$OUTPREFIX.par
