demes=("verus" "ellioti" "troglo" "schwein")
n=${#demes[@]}

for (( i=0; i<=$n-1; i++ )); do
    deme1=${demes[i]}1;
    deme2=${demes[i]}2;
    echo $deme1;
    echo $deme2;
    python3 ~/sumstats.py \
    --bin_dir ~/ \
    -i ~/$1/$1 \
    -o ~/$1/$1\_$deme1\_$deme2 \
    -p ~/stats_noMAC_diplo.spar \
    --stats \
    --ceu $deme1 --yri $deme2 --vindija $deme1; \
    sed -i "s/CEU/$deme1/g" ~/$1/$1\_$deme1\_$deme2.stats ; \
    sed -i "s/YRI/$deme2/g" ~/$1/$1\_$deme1\_$deme2.stats ; \
    for (( j=i+1; j<=$n-1; j++ )); do
        deme1=${demes[i]}1;
        deme2=${demes[j]}1;
        echo $deme1;
        echo $deme2;
        python3 ~/sumstats.py \
        --bin_dir ~/ \
        -i ~/$1/$1 \
        -o ~/$1/$1\_$deme1\_$deme2 \
        -p ~/stats_noMAC_diplo.spar \
        --stats \
        --ceu $deme1 --yri $deme2 --vindija $deme1; \
        sed -i "s/CEU/$deme1/g" ~/$1/$1\_$deme1\_$deme2.stats ; \
        sed -i "s/YRI/$deme2/g" ~/$1/$1\_$deme1\_$deme2.stats ; \
    done;
done
