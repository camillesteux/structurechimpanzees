# rtournebize

mkdir -p sims

# Simulate data
python3 simulate_chimp.py \
-o sims/chimp \
--mut_model 'BinaryMutationModel(False)' \
--algo hudson \
--model chimp \
--seed 42 \
-i chimp.est \
-g 20 20 \
--samples \
Tw__2:15:0:Tw \
Tnc__2:15:0:Tnc \
Te__2:15:0:Te \
Tc__2:15:0:Tc \
P__2:15:0:P \
--mut_rate 1.2e-8 \
--rec_rate 0.7e-8 \
--generation_time 25

# Calculate the D-statistics between all requested group pairs
parallel --colsep " " -a pop_pairs.txt \
-j 6 \
python3 "sumstats_chimp.py \
-p chimp.spar \
--stats --SE \
--bin_dir bin \
-i sims/chimp \
-o sims/chimp.ceu_{1}.yri_{2} \
--ceu {1} --yri {2} --vindija P \
--assume_no_na"

# Plot
Rscript chimp.R

#___
