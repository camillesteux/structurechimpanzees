cd ~/general_model_mscommands/

demes=("verus" "ellioti" "troglo" "schwein")
n=${#demes[@]}

for model in general_nisland_model_*.ms ; do
    mkdir ~/stats/${model%.*};
    bash $model > ~/stats/${model%.*}/${model%.*}.out ;

    #Compute Fst (en pi)
    bash ~/convertmstoeigenstrat_ptroglodytes.sh ${model%.*} ;
    cp ~/pop_ptroglodytes.ind ~/${model%.*}/${model%.*}.ind ;
    bash ~/compute_sumstats_ptroglodytes.sh ${model%.*} ;
    for (( i=0; i<=$n-1; i++ )); do
        deme1=${demes[i]}1;
        deme2=${demes[i]}2;
        file=$(ls ~/${model%.*}/ | grep $deme1"_" | grep $deme2"\.");
        #cat ~/${model%.*}/$file | awk '$1~"pi_" {print $1, $3}' >> ~/${model%.*}/${model%.*}_pi.stats ;
        cat ~/${model%.*}/$file | awk -v deme1=$deme1 -v deme2=$deme2 '$1~"Fst" {print $1"_"deme1"_"deme2, $3}' >> ~/${model%.*}/${model%.*}_Fst.stats ;
        for (( j=i+1; j<=$n-1; j++ )); do
            deme1=${demes[i]}1;
            deme2=${demes[j]}1;
            file=$(ls ~/${model%.*}/ | grep $deme1"_" | grep $deme2"\.");
            #cat ~/${model%.*}/$file | awk '$1~"pi_" {print $1, $3}' >> ~/${model%.*}/${model%.*}_pi.stats ;
            cat ~/${model%.*}/$file | awk -v deme1=$deme1 -v deme2=$deme2 '$1~"Fst" {print $1"_"deme1"_"deme2, $3}' >> ~/${model%.*}/${model%.*}_Fst.stats ;
        done ;
    done
    #cat ~/${model%.*}/${model%.*}_pi.stats | sort -V | uniq > ~/${model%.*}/${model%.*}_pi_sorted.stats
    echo ${model%.*} >> ~/Stats_P_troglodytes.txt ;
    #cat ~/${model%.*}/${model%.*}_pi_sorted.stats >> ~/Stats_P_troglodytes.txt
    cat ~/${model%.*}/${model%.*}_Fst.stats >> ~/Stats_P_troglodytes.txt

    #Compute Ho
    for ind in {10..89} ; do
      res=$(cat ~/${model%.*}/${model%.*}.vcf | grep -v "#" | awk -v ind=$ind '{if ($ind=="0|1" || $ind=="1|0") print $0}' | wc -l);
      het=$(awk "BEGIN {print $res / 100000000}");
      echo $ind" "$res" "$het" "$model >> ~/${model%.*}/obs_het.csv;
    done
    sed -i 's/,/./g' ../obs_het.csv
done
