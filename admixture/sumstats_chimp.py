#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 18:24:33 2020
@author: rtournebize

2     011220   introduced PSMC + LD statistics; all statistics were checked against stats_v2.py results: OK
3     171220   introduced Sprime
4     211220   implemented the haploid version of CRF (as it should be)
5     140121   CRF: now outputs 1 for Neand geno if >=1 derived allele (example files say the contrary but the sup mat says so)
                 S' & CRF & LD : option for MAC-sensitive SNP-downsampling (to mimick 1000G properties), cf Sankararaman 14
                 S' : removed min_MatchRate and min_nSNP_per_segment filtering (they are not used in the paper in fact)
                 S' : outputs a 3rd line: the genome-wide match rate estimate (= single value for whole genome)
                 CRF : considers only the SNPs that are polymorphic in CEU (as done in Sankararaman 14)
                 CRF : introgression rate variation now computed across haploid individuals (instead of across chromosomes)
                 CRF : introduced option to filter out archaic segments < CRF_min_segment_length_cM
                 LD : now considers only the SNPs that are polymorphic in CEU (as done in Sankararaman 12)
6     190121   introduced the pi statistic (manually checked that pi=simulated_theta in a simple constant-No panmic model => OK)
7     220321   added uncertainty computation for the D statistic (added option --SE + option --block to specify block length)
                 added jackknife calculation for ancestry-LD (equivalent Sankararaman's `computedD`) + outputs *.cov.txt.gz only if --verbose
                 linked to remi_functions_v3 (instead of _v2) + outputs *.sprime.txt.gz only if --verbose
                 new ancestry-LD configuration was checked with simulations, cf. checks/LD => ok
                 added a sys.exit() in case of the presence of a 9 in the dataset -- we currently do not allow missing data
8     010421   added f3 calculation (outputs into *.stats)
                 changed Sprime records and outputs into a novel *.sprime2.stats output which contains different new statistics (eg. Detection Rate) and gives estimates
                 by varying the Minimum Match Rate for fragments to be considered introgressed (ie. to remove False Positives)
9     020521   changed the *.psmc.stats output to a single run: will systematically output the last iteration + few minor changes + added `--version`
10    190521   found that CRF could only work if the -o path was relative
                 CRF: reformatted the whole script (incl. debugged when no MAC downsampling (was reading single chromosome)) + no subset on *.crf.filtered.snps (already performed) + read gpos directly from *.SNP.snp
                 LD: added expfit with jacquelin inferred starting values + export all the jackknife runs
11    200521   stats: added a switch so that DCFS = NA if the number of target samples is greater or equal to `ndiploids_for_dcfs` (usually :=5)
                 LD: added the possibility to compute the single-sample LD statistic (Moorjani et al 2016) => introduced the --ld1 switch to call this mode
                 pseudodiploidization: functionality was moved out of the stats section (to be more general) => before, we could not do --ld without --stats
11.1  250521   bug in the condition leading to call pseudodiploidize(): `(do_stat)` instead of `(PAR["D_do_pseudodiploid"]=="YES" and do_stats)`
11.2  110621   LD: now outputs the number of SNP pairs and the NRMSD in the output (tried a convergence assessment approach, but not discriminative)
11.3  140621   can read only a subset of the SNPs using the `--chrom` option
                  LD: outputs on last line the sequence of age estimates for increasing number of chromosomes
11.4  180621   LD: added `--ld_Nea_ancestral` to ascertain SNPs where Nea has the ancestral allele (instead of derived, cf. Sankararaman 2012)
12.0  190621   major rewriting of the LD part: now --ld and --ld1 are in separate sections
               canceled the 11.4 change (ie removed the --ld_Nea_ancestral option)
               LD: now outputs 3 files corresponding to the three ascertainment schemes proposed in Sankararaman 2012
               LD: new parameter slot: "LD_sample_haploid_Neanderthal" in the stats parfile, which allows sampling a single chromosome from the geno1 file
                 |__. This is useful to exactly reproduce the Sankararaman's simulation procedure.
               stats: new parameter in par file "AFS_num_inds 10" (by default: 10) to restrict the number of inds on which to compute the AFSs
               general: added option `--bin_dir`, indicating the path to the dir containing the softwares (eg. sprime.jar) specified in the stat parfile
               LD & LD1: the block size is now defined as the seq length (instead of the number of SNPs polymorphic in CEU (v<12.0))
               checked by comparing results with sumstats 1.4: ok
               "psmc_samples_0based" parameter: changed the syntax: now write SampleIdx0Based:Label => this will write a file `output.Label.psmc2.stats`
12.1  240621   LD: bug in the Moorjani_0 ascertainment (missing variable name) + in --ld1, existed incremental_assessment but impossible for single sample LD
12.2  080921   Fst: following the Fst line in the *.stats output, now add one line: Fst calculated for SNPs only polymorphic in YRI (as in Bhatia 13)
                    added the option --Fst_ascertainment in relation to this (by default: False, meaning that the alternative ascertainment is not applied prior to Fst calculation)
13.0  120921   now able to handle chromosome of heterogeneous sizes -- in which case, the seq_length but be set to -1 in the *.par file
               the script will then use the last position of the chromosome as the chrom size -- note that improvement can be made setting chrom size as the last position in the genetic map (more accurate but in the end should change nearly nothing)
13.1  281121   now outputs also the median length for *.sprime2.stats
13.2  081221   was missing two nan to print per line in *.sprime2.stats, if no tracts detected
13.3  220322   now feeds read_fixedwidth with file handler instead of file name, because was throwing an error in latest pandas version (python3.10)
               generate_CRF_input: print the mean density of effective SNPs fed to CRF
14.0  270422   possibility now to specify the pops to analyze by argparse options
14.1  270422   included the --f4 option to calculate f4-ratios, this was checked in comparison with admixr results using simulations (also D-stats)
               checked 14.1 results vs. 13.3
14.2  020522   f4-ratio: moving window calculation with np.nansum instead of np.sum (in 14.1)
14.3  200622   now if `-v`, outputs the length of all S' segments
14.3.1  031122  microbug on get_sprime_stats() whereby script would through an error when there was no fragment at a particular min threshold (because forgot to return l, m in the exception)
15.0  111222   integrated the {extras.py} script v6.5 to avoid ext dependency
16.0  211222   added a new step to **subset sites** (either randomly to retain only `n` sites or using ascertainments);
               occurs _after_ the EIGENSTRAT data import step and _after_ the `na_rm` step (cf below)
17.0  291222   added option `--na_rm` which removes any SNP containing any `9` (missing genotype) and moved this step just before the SNP-subsetting step (`--subset`)
               added option `--assume_no_na` to speed up the calculation
               added a condition at the start of do_f4 to use the allele count matrices generated from do_stats, to speed up code
17.1  291222   removed some functions imported from previous extras.py script, which are useless here
18.0  291222   added option `--ancestral` to specify one single individual from the input *.ind file which specifies the ancestral genotypes
18.1  030122   debugged a small bug which occurred sometimes when generating PSMC output, when the very last bin was hetero (in prev version, would throw a stop error)
19.0  160123   introduced option `--ld_as_sanka12` which allows to calculate the ancestry-LD nearly exactly as in Sankararaman 2012; checked using simulations with computed vs. this script and agreement to 99%
               under `--ld_as_sanka12`, what changes:
                   - uses the biased covariance stat (divided by N) [instead of the unbiased (divided by N+1)]
                   - binned distances have a +1-incremented index [instead of non-incremented index]
                   - default value (for missing cov values) is 0 [instead of nan]
                   - do not set all cov values to nan if more than 50% of the bins have missing or infinite cov values [instead of setting them all to nan]
19.0.1  030523  added option to argparse to print out default values when using --help
19.0.2  040723  under --verbose mode, now keeps the intermediate *.out files produced by LD_1 (single_sample_covariance())
19.1.0  060723  implemented the statistics `--csfs` + added `SPECIFICITIES` subsection in this commented header

______________________

_____OBSERVATIONS_____
______________________

240521. MAC-sensitive downsampling should not be performed for single-sample LD stat. Pseudodiploidization has very minor impact on LD decay rate (for both single and multi-sample stat). MAC-sensitive downsampling for multi-sample LD stat has little impact on decay rate estimates. Filtering out monomorphs has little impact on the decay rate estimates.

______________________

________PAR FILE______
______________________
Options:

- snp_subsetting ::: optional ::: if not added (or if set as `NO`), will never subset SNPs
                                    otherwise, if equal to a number, will subset SNPs to this number
                                    otherwise, if a string, will subset SNPs using an ascertainment scheme:
                                        * Archaic_array : at least one Neanderthal allele differs from
                                                          the majority allele in a panel of 24 YRI samples
                                                          (Fu et al 2015, Moorjani et al 2016)
                                 *NOTE* The subsetting is down *after* the whole data importation & optional sequence subsetting (`--chom`)
                                        & *before* MAC subsetting and pseudodiplo conversion

______________________

____SPECIFICITIES_____
______________________


- MAC-filtering affects only:
        S'
        CRF
        LD
        CSFS

- for --csfs, the number of retained diploids is the one specified as `CSFS_num_inds` in the *.par file (-p)


"""

___version___ = '19.1.0'

import numpy as np
import sys, re, argparse, time, allel, os, gzip, copy
from os import path
import pandas as pd
from decimal import Decimal
from scipy.spatial.distance import pdist

np.seterr(divide = 'ignore', invalid = 'ignore')

############################################################
############################################################
############################################################
############################################################
############################################################

def global_params():
    ndiploids_for_dcfs = 5
    # { minor_allele_count: proba_of_SNP_acceptance } from Sankararaman et al 2014 (SI2.1)
    acceptance_rates = {0: 0,
                    1: 0.25,
                    2: 0.5,
                    3: 0.75,
                    4: 0.8,
                    5: 0.9,
                    6: 0.95,
                    7: 0.96,
                    8: 0.97,
                    9: 0.98,
                    10: 0.99 }
    return ndiploids_for_dcfs, acceptance_rates

###############################################################################

def parse_options():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', action='version', version='%(prog)s: '+str(___version___))
    parser.add_argument('-i', '--input_prefix', type=str, required=True, help="Prefix of the input genetic data files")
    parser.add_argument('-o', '--output_prefix', type=str, required=False, default=None, help="Prefix of the output files")
    parser.add_argument('-p', '--parameter_file', type=str, required=True, help="Path to the statistical parameter file")
    parser.add_argument('--stats', action="store_true", default=False, help="Compute classical summary statistics?")
    parser.add_argument('--psmc', action="store_true", default=False, help="Compute PSMC curves?")
    parser.add_argument('--csfs', action="store_true", default=False, help="Compute CSFS?")
    parser.add_argument('--ld', action="store_true", default=False, help="Compute the ancestry LD on a multi-sample mode.")
    parser.add_argument('--ld1', action="store_true", default=False, help="Compute the ancestry LD on a single-sample mode.")
    parser.add_argument('--sprime', action="store_true", default=False, help="Compute S' analysis?")
    parser.add_argument('--crf', action="store_true", default=False, help="Compute CRF analysis?")
    parser.add_argument('--SE', action="store_true", default=False, help="Compute block jackknife standard errors?")
    parser.add_argument('--f4', action="store_true", default=False, help="Compute f4-ratios?")
    parser.add_argument('--block', type=str, default="5000kb", required=False, help="Block length for jackknife, either *kb for length in kb or *snp for length in a fixed number of SNPs, where * is a numeric value. A block length of 5 Mbp is used for D calculation in Fu et al. 2015, SI 16.")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, help="Add this option to output additional files.")
    parser.add_argument('--chrom', type=int, nargs="*", default=[None], required=False, help="List of chromosomes to analyze.")
    parser.add_argument('--bin_dir', type=str, required=True, help="Path to the directory containing the binaries/softwares.")
    parser.add_argument('--Fst_ascertainment', action="store_true", default=False, help="Add this switch to calculate Hudson Fst only with SNPs polymorphic in YRI.")
    # Samples
    parser.add_argument('--ceu', type=str, required=False, help="CEU population label, as in the *.ind file; alternatively, in the param file: CEU_indices_0based", default = None)
    parser.add_argument('--yri', type=str, required=False, help="YRI population label, as in the *.ind file; alternatively, in the param file: YRI_indices_0based", default = None)
    parser.add_argument('--vindija', type=str, required=False, help="Vindija_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Vindija_indices_0based", default = None)
    parser.add_argument('--altai', type=str, required=False, help="Altai_Neanderthal population label, as in the *.ind file; alternatively, in the param file: Altai_indices_0based", default = None)
    parser.add_argument('--psmc_pops', type=str, required=False, nargs="*", help="Labels of the populations to analyze with PSMC, will select the first ind of each pop; alternatively, in the param file: psmc_samples_0based", default = None)
    parser.add_argument('--na_rm', action="store_true", default=False, help="Remove missing genotypes?") # v17.0
    parser.add_argument('--assume_no_na', action="store_true", default=False, help="Should the script assume there is no missing genotypes in the data?") # v17.1
    parser.add_argument('--ancestral', type=str, required=False, default = None, nargs=1, help="Population label, as in the *.ind file, of the single individual reference specifying the ancestral genotypes: note that we will polarize all alleles only if the ancestral genotypes are unambiguously 0 or 2") # v18.0
    parser.add_argument('--ld_as_sanka12', action='store_true', default=False, help="Add this switch to calculate the population-wise ancestry-LD with the same algorithm and parameters than Sankararaman et al. 12 original study, i.e. (1) using biased covariance stat (divide by N) instead of unbiased (divide by N+1); (2) binning distances with +1-incremented index instead of non-incremented index; (3) default value as 0 instead of nan; (4) do not set all cov values to nan if more than 50 percents of the bins have missing or infinite cov values") # v19.0
    return parser.parse_args()

def downsample(geno, snp, genetic_pos, acceptance_rates, generator, geno1 = None, columns = None):
    if columns is not None:
        nind = len(columns)
        MAC = np.sum(geno[:,columns], axis = 1).astype(int)
        MAC[MAC>nind] = 2*nind-MAC[MAC>nind]
    else:
        nind = geno.shape[1]
        MAC = np.sum(geno, axis = 1).astype(int)
        MAC[MAC>nind] = 2*nind-MAC[MAC>nind]
    MAC[MAC>10] = 10
    probs = np.array([acceptance_rates[x] for x in MAC])
    probs = generator.binomial(n = 1, p = probs)
    subset = np.array([False]*len(probs))
    subset[probs==1] = True
    print("Proportion of MAC-informed SNPs accepted:     "+str(np.round(np.true_divide(np.sum(subset), len(subset))*100,1))+"%")
    if geno1 is None:
        if genetic_pos is None:
            return geno[subset,:], snp[subset,:], None
        else:
            return geno[subset,:], snp[subset,:], genetic_pos[subset]
    else:
        if genetic_pos is None:
            return geno[subset,:], snp[subset,:], None,                geno1[subset,:]
        else:
            return geno[subset,:], snp[subset,:], genetic_pos[subset], geno1[subset,:]

def count(geno, sam_labels, label, ncol = None):
    # matrix of 2 columns
    # first column:  count of the derived alleles
    # second column: count of the ancestral alleles
    cols = np.where(sam_labels==label)[0]
    if ncol is not None:
        cols = cols[0:int(ncol)]
    x = np.sum(geno[:,cols], axis = 1)
    x = np.column_stack((x, 2*len(cols)-x))
    return x

def union(list_of_doublets):
    b = []
    for begin,end in sorted(list_of_doublets):
        if b and b[-1][1] >= begin - 1:
            b[-1][1] = max(b[-1][1], end)
        else:
            b.append([begin, end])
    return b

def P(x, n=6):
    f = "{0:."+str(int(n))+"f}"
    return f.format(x)

def populate(file, PAR):
    with open(file, "r") as FIN:
        for line in FIN:
            if line.startswith("#") or line=="\n":
                continue
            line = line.strip().replace("\t", " ")
            line = re.split("\s+", line)
            if ("_0based" in line[0]) and (line[0].startswith("psmc_")==False):
                if len(line)==2:
                    PAR[line[0]] = np.array([int(line[1])])
                else:
                    PAR[line[0]] = np.array([int(x) for x in line[1:]])
            elif ("_0based" in line[0]) and (line[0].startswith("psmc_")==True):
                if ":" not in line[1]:
                    sys.exit("New since v12: `psmc_samples_0based` must be of the form SampleIndex0Based:Label.")
                if len(line)==2:
                    PAR[line[0]] = np.array([line[1]])
                else:
                    PAR[line[0]] = np.array(line[1:])
            else:
                PAR[line[0]] = line[1]
    return PAR

def gpos_check(PAR, genetic_pos, snp, output_unit = "cM", verbose = True):
    if PAR["genetic_pos_unit"].lower().startswith("m"):
        if output_unit=="cM":
            unit = Decimal("100.0")
        else:
            unit = Decimal("1.0")
        if verbose and genetic_pos is not None:
            print("Unit of genetic pos:   Morgans")
    elif PAR["genetic_pos_unit"].lower().startswith("c"):
        if output_unit=="cM":
            unit = Decimal("1.0")
        else:
            unit = Decimal("0.01")
        if verbose and genetic_pos is not None:
            print("Unit of genetic pos:   centiMorgans")
    else:
        sys.exit("Unit of genetic positions not recognized.")

    # if no genetic_position provided in the snp file
    if genetic_pos is None:
        if verbose:
            print("No genetic positions provided: inferring using the uniform recombination rate r="+str(PAR["recomb_rate"])+".")
        rec = Decimal(str(PAR["recomb_rate"]))
        if output_unit=="cM":
            genetic_pos = np.array([Decimal("100")*(x*rec) for x in snp[:,1]])
        else:
            genetic_pos = np.array([x*rec for x in snp[:,1]])

    # if provided
    else:
        genetic_pos = unit*genetic_pos

    return genetic_pos

def rm_if_exists(file):
    if path.exists(file):
        os.remove(file)

def summary(x):
    a = ["n", "Min", "Q2.5p", "Q25p", "Median", "Mean", "Std", "Q75p", "Q97.5p", "Max"]
    if len(x)==0:
        b = [0]+["None"]*(len(a)-1)
    else:
        b = [len(x), np.min(x), np.quantile(x, 0.025), np.quantile(x, 0.25), np.median(x), np.mean(x),
             np.sqrt(np.var(x)), np.quantile(x, 0.75), np.quantile(x, 0.975), np.max(x)]
    return a, b

def NRMSD(y, yfit):
    rmsd = np.sqrt(np.mean((np.array(y) - np.array(yfit))**2))
    nrmsd = rmsd * (np.max(yfit) - np.min(yfit))**-1
    return nrmsd

def interweave(a, b):
    c = np.empty((a.size + b.size,), dtype=a.dtype)
    c[0::2] = a
    c[1::2] = b
    return c

def get_haplotypes(geno, geno1):
    geno2 = geno - geno1
    return interweave(geno1, geno2)

def get_pos(string):
    return int(string.split(":")[1])

def get_gpos(snp, gpos, seq, phys_pos):
    z = np.where(np.logical_and(snp[:,0]==seq, snp[:,1]==phys_pos))[0]
    if len(z)==0:
        sys.exit("Error multiple gpos for one ppos")
    elif len(z)==1:
        z = gpos[z]
    else:
        z = gpos[z[0]]
    return z

def pseudodiploidize(M, generator):
    geno = copy.deepcopy(M)
    rep = geno==1
    geno[rep] = generator.choice([0,2], size=np.sum(rep), replace=True).astype(np.int8)
    return geno

def main():
    ndiploids_for_dcfs, acceptance_rates = global_params()

    start_time = time.time()

    options = parse_options()
    prefix = options.input_prefix
    output = options.output_prefix
    parfile = options.parameter_file
    Chrom = options.chrom
    # statistic switches
    do_stats = options.stats
    do_psmc = options.psmc
    do_ld = options.ld
    do_ld1 = options.ld1
    if do_ld1 and do_ld:
        sys.exit("You can provide only one LD mode (either --ld or --ld1)")
    do_sprime = options.sprime
    do_crf = options.crf
    do_SE = options.SE
    do_f4 = options.f4
    bin_dir = options.bin_dir
    # block stats
    block_length = options.block
    # v12.2
    do_Fst_ascertainment = options.Fst_ascertainment
    # v17.0
    na_rm = options.na_rm
    assume_no_na = options.assume_no_na
    # v18.0
    ancestral = options.ancestral
    if assume_no_na and na_rm:
        sys.exit("You cannot specify altogether --na_rm and --assume_no_na")
    # v19.1.0
    do_csfs = options.csfs

    print("v"+str(___version___))
    print("allel: "+str(allel.__version__)+"\n")

    if output is None:
        output = prefix

    if not path.exists(prefix+".geno") and not path.exists(prefix+".geno.gz"):
        sys.exit(prefix+".geno(.gz) does not exist")
    if not path.exists(prefix+".snp") and not path.exists(prefix+".snp.gz"):
        sys.exit(prefix+".snp(.gz) does not exist")
    if not path.exists(prefix+".ind") and not path.exists(prefix+".ind.gz"):
        sys.exit(prefix+".ind(.gz) does not exist")
    if path.exists(prefix+".geno") and path.exists(prefix+".geno.gz"):
        sys.exit(prefix+".geno and .geno.gz both exist")
    if path.exists(prefix+".snp") and path.exists(prefix+".snp.gz"):
        sys.exit(prefix+".snp and .snp.gz both exist")
    if not path.exists(parfile):
        sys.exit(parfile+" does not exist")

    # read parameters
    PAR = {"CEU_indices_0based": None, "YRI_indices_0based": None, "Vindija_index_0based": None, "Altai_index_0based": None, "seed": None, "psmc_samples_0based": None}
    PAR = populate(parfile, PAR)
    PAR = populate(prefix+".par", PAR)
    PAR["ld_as_sanka12"] = options.ld_as_sanka12

    # checks (v14.0)
    IndFile = pd.read_csv(prefix+".ind", sep = "\s+", usecols = [2], engine = "c", dtype = str, na_filter = None, header = None).to_numpy()[:,0]
    ## all
    optSam = {"CEU": options.ceu, "YRI": options.yri, "Vindija": options.vindija, "Altai": options.altai}
    for key in optSam.keys():
        if key=="Vindija" or key=="Altai":
            sp = "index"
        else:
            sp = "indices"
        if PAR[key+"_"+sp+"_0based"] is None and optSam[key]==None and key!="Altai":
            sys.exit("Missing specification of "+key+" samples")
        elif PAR[key+"_"+sp+"_0based"] is not None and optSam[key]!=None:
            sys.exit("You cannot specify altogether "+key+"_"+sp+"_0based in the par file and in the command line.")
        elif PAR[key+"_"+sp+"_0based"] is None and optSam[key]!=None:
            PAR[key+"_"+sp+"_0based"] = np.where(IndFile==optSam[key])[0]
            if len(PAR[key+"_"+sp+"_0based"]) == 0:
                sys.exit("We found 0 samples for "+key)
    ## f4
    if do_f4: assert PAR["Altai_index_0based"] is not None and PAR["Vindija_index_0based"] is not None, "If --f4, you must specify two Neanderthal samples."
    if do_f4: assert do_stats, "If --f4, you must also add the switch --stats."
    ## psmc
    if options.psmc_pops is None and PAR["psmc_samples_0based"] is None and do_psmc:
        sys.exit("Missing specification of the PSMC samples")
    elif options.psmc_pops is not None and PAR["psmc_samples_0based"] is not None:
        sys.exit("You cannot specify altogether psmc_samples_0based in the par file and in the command line.")
    elif options.psmc_pops is not None and PAR["psmc_samples_0based"] is None:
        PAR["psmc_samples_0based"] = []
        for p in options.psmc_pops:
            print(p)
            PAR["psmc_samples_0based"].append(str(np.where(IndFile==p)[0][0])+":"+p)

    if PAR["seed"]=="None":
        seed = None
    else:
        seed = int(PAR["seed"])
    rng = np.random.default_rng(seed = seed)

    if not "LD_sample_haploid_Neanderthal" in PAR.keys():
        PAR["LD_sample_haploid_Neanderthal"] = "NO"

    nCEU, nYRI, nVindija = len(PAR["CEU_indices_0based"]), len(PAR["YRI_indices_0based"]), len(PAR["Vindija_index_0based"])
    if PAR["Altai_index_0based"] is None:
        nAltai = 0
    else:
        nAltai = len(PAR["Altai_index_0based"])
    samples1 = np.concatenate((PAR["CEU_indices_0based"], PAR["YRI_indices_0based"], PAR["Vindija_index_0based"]))
    if PAR["Altai_index_0based"] is not None:
        samples1 = np.concatenate((samples1, PAR["Altai_index_0based"]))
    samples2 = np.array([0]*nCEU + [1]*nYRI + [2]*nVindija + [3]*nAltai)
    # CEU: 0
    # YRI: 1
    # Vindija: 2
    # Altai: 3

    if ancestral is not None:
        ancestral_idx = np.where(IndFile==ancestral)[0]
        if len(ancestral_idx) != 1:
            sys.exit("Error. There must be only a single ancestral individual.")
        ancestral_idx = ancestral_idx[0]
        samples1 = np.concatenate((samples1, np.array([ancestral_idx])))
    else:
        ancestral_idx = None

    print("\n____SUMMARY STATISTICS___\n")
    print("Parameter file:        "+parfile)
    print("Input prefix:          "+prefix)
    print("Output prefix:         "+output)
    if Chrom[0] is not None:
        print("Sequences to analyze:  "+" ".join([str(x) for x in Chrom]))
    print("YRI idx:               "+"(n="+str(len(PAR["YRI_indices_0based"]))+") "+" ".join([str(x) for x in PAR["YRI_indices_0based"]]))
    print("CEU idx:               "+"(n="+str(len(PAR["CEU_indices_0based"]))+") "+" ".join([str(x) for x in PAR["CEU_indices_0based"]]))
    print("Vindija idx:           "+"(n="+str(len(PAR["Vindija_index_0based"]))+") "+" ".join([str(x) for x in PAR["Vindija_index_0based"]]))
    if nAltai > 0:
        print("Altai idx:             "+"(n="+str(len(PAR["Altai_index_0based"]))+") "+" ".join([str(x) for x in PAR["Altai_index_0based"]]))
    print("Ancestral idx:         "+str(ancestral_idx))
    print("Pseudodiploid for D:   "+str(PAR["D_do_pseudodiploid"]))
    print("Seed:                  "+str(seed))
    print("Calculating SEs:       "+str(do_SE))
    print("Block length:          "+block_length)
    print("Software directory:    "+bin_dir)
    print("Ascertained Fst:       "+str(do_Fst_ascertainment))
    print("Assume no NA:          "+str(assume_no_na))
    print("Remove SNPs with NA:   "+str(na_rm))
    print("")

    ind, snp, genetic_pos, geno = read_eigenstrat(prefix,
                             read_geno = True,
                             read_snp = True,
                             read_ind = True,
                             ancestral_allele_encoded_as_0 = True,
                             columns = samples1)

    # subset, if requested
    if Chrom[0] is not None:
        chrom_ok = np.in1d(snp[:,0], np.array(Chrom), assume_unique = False)
        snp = snp[chrom_ok,:]
        if genetic_pos is not None:
            genetic_pos = genetic_pos[chrom_ok]
        geno = geno[chrom_ok,:]

    seq_names = np.unique(snp[:,0], return_index = True)
    seq_end = snp[np.concatenate((seq_names[1][1:]-1, np.array([-1]))),:][:,1]
    seq_names = seq_names[0]
    n_seq = len(seq_names)

    # estimate sequence lengths, if heterogeneous
    if float(PAR["seq_length"]) == -1.:
        seq_lengths = { seq_names[sn]: float(seq_end[sn]) for sn in range(n_seq) }
    else:
        seq_lengths = { seq_names[sn]: float(PAR["seq_length"]) for sn in range(n_seq) }
    print("Sequence lengths:")
    for sn in seq_names:
        print(str(sn)+": "+str(np.round(seq_lengths[sn]*1e-6, 2))+" Mbp")

    # estimate genome size
    genome_size = float(sum(list(seq_lengths.values())))
    print("=========================")
    print(str(np.round(genome_size*1e-6, 2)),"Mbp ("+str(int(genome_size))+" bp)\n")

    if (do_ld==True) or (do_sprime==True) or (do_crf==True) or (do_ld1==True):
        genetic_pos = gpos_check(PAR, genetic_pos, snp, output_unit = "cM", verbose = True)

    if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
        print("Reading *.geno1.")
        if not path.exists(prefix+".geno1") and not path.exists(prefix+".geno1.gz"):
            sys.exit(prefix+".geno1(.gz) does not exist")
        if path.exists(prefix+".geno1") and path.exists(prefix+".geno1.gz"):
            sys.exit(prefix+".geno1 and .geno1.gz both exist")
        if path.exists(prefix+".geno1.gz"):
            FGENO = gzip.open(prefix+".geno1.gz", "rt")
        else:
            FGENO = open(prefix+".geno1", "r")
        geno1 = read_fixedwidth(FGENO, dtype = np.int8, columns = samples1)
        FGENO.close()
        if Chrom[0] is not None:
            geno1 = geno1[chrom_ok,:]
            del chrom_ok
        if geno.shape != geno1.shape:
            sys.exit("geno and geno1 must have the same dimensions.")

    print("Reading done.\n")

    # polarize alleles
    if ancestral is not None:
        print("Polarizing alleles according to ancestral reference:")
        anc_gt = geno[:,-1]
        anc_gt_0 = anc_gt==0
        anc_gt_2 = anc_gt==2
        #
        geno[anc_gt_2,0:-1] = 2-geno[anc_gt_2,0:-1]
        geno[geno<0] = 9
        ok = np.logical_or(anc_gt_0, anc_gt_2)
        nori = geno.shape[0]
        geno = geno[ok,0:-1]
        if genetic_pos is not None:
            genetic_pos = genetic_pos[ok]
        snp = snp[ok,:]
        if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
            geno1 = geno1[ok,0:-1]
            geno1[geno1<0] = 9
        del ok, anc_gt, anc_gt_0, anc_gt_2
        print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
        print("   Polarized: "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
        #
        samples1 = samples1[0:-1]

    # handle missing data in the *.geno file
    if assume_no_na:
        print("_/!\__ Assumes that there are NO missing genotypes in the input *.geno file.")
    else:
        # remove all SNPs containing at least one missing genotype (9)
        bads = np.any(geno==9, axis = 1)
        if np.any(bads==True):
            if na_rm:
                nori = geno.shape[0]
                ok = (bads==False)
                geno = geno[ok,:]
                if genetic_pos is not None:
                    genetic_pos = genetic_pos[ok]
                snp = snp[ok,:]
                if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                    geno1 = geno1[ok,:]
                del ok, bads
                print("Removing SNPs having missing genotypes:")
                print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
                print("   Retained:  "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
            else:
                sys.exit("Error. Presence of missing data. We currently do not handle missing genotypes.\n")
        else:
            pass

    # if requested, subset sites
    if "snp_subsetting" in PAR.keys():
        if PAR["snp_subsetting"].upper() == "NO":
            print("Not subsetting sites.\n")
            pass
        else:
            if PAR["snp_subsetting"].isdigit():
                if int(PAR["snp_subsetting"]) >= snp.shape[0]:
                    subSites = range(snp.shape[0])
                else:
                    subSites = np.sort(rng.choice(range(snp.shape[0]), int(PAR["snp_subsetting"])))
            elif PAR["snp_subsetting"] == "Archaic_array":
                if not 3 in samples2:
                    sys.exit("You must sample an Altai Neand __and__ a Vindija Neand for `Archaic` snp_subsetting")
                sys.exit()
            else:
                sys.exit("snp_subsetting not implemented.")
            #
            nori = geno.shape[0]
            geno = geno[subSites,:]
            if genetic_pos is not None:
                genetic_pos = genetic_pos[subSites]
            snp = snp[subSites,:]
            if (do_crf==True) or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                geno1 = geno1[subSites,:]
            del subSites
            print("Subsetting sites based on: `snp_subsetting="+PAR["snp_subsetting"]+"`")
            print("   Removed:   "+str(nori-geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(nori-geno.shape[0],nori),1))+"%)")
            print("   Retained:  "+str(geno.shape[0])+"   ("+str(np.round(100. * np.true_divide(geno.shape[0],nori),1))+"%)\n")
    else:
        print("Not subsetting sites.\n")

    # MAC-sensitive SNP-downsampling
    # To mimick 1000G, the MAC should be computed across sapiens populations only (excluding Neanderthal)
    if do_sprime or do_crf or do_ld or do_csfs:
        if (PAR["do_MAC_sensitive_SNP_downsampling"].upper() == "YES") or (do_csfs): # v19.1.0, for csfs we calculate separately with and without MAC-filtering
            columns = np.concatenate((np.where(samples2==0)[0], np.where(samples2==1)[0]))
            if do_crf or (do_ld==True and PAR["LD_sample_haploid_Neanderthal"].upper()=="YES"):
                xs_geno, xs_snp, xs_genetic_pos, xs_geno1 = downsample(geno, snp, genetic_pos, acceptance_rates, rng, geno1 = geno1, columns = columns)
            else:
                xs_geno, xs_snp, xs_genetic_pos = downsample(geno, snp, genetic_pos, acceptance_rates, rng, geno1 = None, columns = columns)

    # Pseudodiploidize
    if (do_stats) or (PAR["LD_do_pseudodiploid"]=="YES" and (do_ld or do_ld1)) or (do_csfs):
        geno_pd = pseudodiploidize(geno, rng)

    ################################################################################################
    if do_stats:

        if "AFS_num_inds" not in PAR.keys():
            PAR["AFS_num_inds"] = 10
        else:
            PAR["AFS_num_inds"] = int(PAR["AFS_num_inds"])
        print("Nb diploids for AFS:   "+str(PAR["AFS_num_inds"]))

        # count matrices have 2 columns: the first gives the count of the derived allele, the second the count of the ancestral allele
        ac_ceu = count(geno, samples2, 0)
        ac_yri = count(geno, samples2, 1)
        ac_nea = count(geno, samples2, 2)
        ac_ref = ac_nea*0
        ac_ref[:,1] = 2

        FOUT = open(output+".stats", "w")
        FOUT.write("#sumstats: "+str(___version___)+"\n")

        ################################################################################################
        # D-statistic

        # block length
        if do_SE:
            if block_length.endswith("kb"):
                bl = block_length.replace("kb", "")
                if not bl.isdigit():
                    sys.exit("D> Block length must be of the form FLOATkb or INTsnp")
                bl = float(bl)
                bl = int(bl * 1000. * np.true_divide(geno.shape[0], genome_size))
                print("Block size (#SNPs):    "+str(bl))
            elif block_length.endswith("snp"):
                bl = block_length.replace("snp", "")
                if not bl.isdigit():
                    sys.exit("Block length must be of the form FLOATkb or INTsnp")
                bl = int(bl)
            else:
                sys.exit("Block length must be of the form FLOATkb or INTsnp")

        # ABBA-BABA = D statistic
        if PAR["D_do_pseudodiploid"].upper() == "YES":
            label = "D_pseudodiploid"
            ac_ceu_pd = count(geno_pd, samples2, 0)
            ac_yri_pd = count(geno_pd, samples2, 1)
            ac_nea_pd = count(geno_pd, samples2, 2)
            if not do_SE:
                num, den = allel.patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref)
                D = np.true_divide(np.sum(num), np.sum(den))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(D)+"\n")
            else:
                num, den = allel.patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref)
                d1, d2, d3, d4, d5 = allel.average_patterson_d(ac_ceu_pd, ac_yri_pd, ac_nea_pd, ac_ref, bl)
                print("Number of blocks:      "+str(len(d4)))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(d1)+" "+P(d2)+" "+P(d3)+"\n")
        else:
            label = "D_diploid"
            if not do_SE:
                num, den = allel.patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref)
                D = np.true_divide(np.sum(num), np.sum(den))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(D)+"\n")
            else:
                num, den = allel.patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref)
                d1, d2, d3, d4, d5 = allel.average_patterson_d(ac_ceu, ac_yri, ac_nea, ac_ref, bl)
                print("Number of blocks:      "+str(len(d4)))
                FOUT.write(label+" "+str(int(np.sum(den)))+" "+P(d1)+" "+P(d2)+" "+P(d3)+"\n")

        ################################################################################################
        # Fst
        # _NOTA BENE_ The Fst measure is insensitive to the monomorphic SNPs
        # As a consequence, the Fst here calculated is identical as if calculated on the SNPs polymorphic in both YRI and CEU
        num, den = allel.hudson_fst(ac_ceu, ac_yri)
        fst = np.true_divide(np.sum(num), np.sum(den))
        FOUT.write("Fst "+str(int(len(den)))+" "+P(fst)+"\n")

        # v12.2:
        if do_Fst_ascertainment:
            # for SNP only polymorphic in YRI:
            ok = np.sum(ac_yri==0, axis = 1) == 0
            num, den = allel.hudson_fst(ac_ceu[ok,:], ac_yri[ok,:])
            fst = np.true_divide(np.sum(num), np.sum(den))
            FOUT.write("Fst_Polymorphic_YRI "+str(int(len(den)))+" "+P(fst)+"\n")

        ################################################################################################
        # Nucleotide diversity (pi) = proportion of average pairwise difference between chromosomes
        pi_ceu = np.true_divide(np.sum(allel.mean_pairwise_difference(ac_ceu)), genome_size)
        pi_yri = np.true_divide(np.sum(allel.mean_pairwise_difference(ac_yri)), genome_size)
        FOUT.write("pi_CEU_per_kb "+str(int(genome_size))+" "+str(P(pi_ceu*1e3))+"\n")
        FOUT.write("pi_YRI_per_kb "+str(int(genome_size))+" "+str(P(pi_yri*1e3))+"\n")

        FOUT.close()

    ################################################################################################
    print("\n____End____")
    print("Analyzed:   "+prefix)
    print("Output:     "+output)
    print(str(np.round((time.time()-start_time)/60, 1))+" min\n")

#_____________________________________________________________________________#
#_____________________________________________________________________________#
#_____________________________________________________________________________#
#                               EXTRAS v6.5                                   #
#_____________________________________________________________________________#
#_____________________________________________________________________________#
#_____________________________________________________________________________#

class HiddenPrints:
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout

# Replace genotypes encoded as 2 into 0, and encoded as 0 into 2
def invert_ancestrality(x):
    y = copy.deepcopy(x)
    dix = {0:2, 2:0, 1:1, 9:9}
    return np.vectorize(dix.get)(y)

# Quickly read fixed-width file (typically a .geno file)
def read_fixedwidth(file, dtype, columns):
    X = pd.read_fwf(file, colspecs=[(x,x+1) for x in columns], dtype = dtype, header = None, na_filter = False)
    X = X.to_numpy()
    return X

# Reshape EIGENSTRAT matrices
def reshape_EIGENSTRAT(x, isgeno = False):
    if len(x.shape) == 1:
        x = np.array([x])
    if isgeno == True and x.shape[0] == 1:
        x = np.transpose(x)
    return x

# Read an EIGENSTRAT file
def read_eigenstrat(prefix,
                    read_geno = True, read_snp = True, read_ind = True,
                    columns = None,
                    ancestral_allele_encoded_as_0 = True,
                    verbose = False):

    if type(columns) == type(None):
        read_specific_cols = False
    else:
        read_specific_cols = True
    if type(columns) == type(list()):
        columns = np.array(columns)

    if path.exists(prefix+".snp") and path.exists(prefix+".snp.gz"):
        sys.exit("There cannot be both a .snp and a .snp.gz file.")
    if path.exists(prefix+".geno") and path.exists(prefix+".geno.gz"):
        sys.exit("There cannot be both a .geno and a .geno.gz file.")
    if path.exists(prefix+".snp.gz"):
        snp_outfix = ".snp.gz"
    else:
        snp_outfix = ".snp"
    if path.exists(prefix+".geno.gz"):
        geno_outfix = ".geno.gz"
    else:
        geno_outfix = ".geno"

    if verbose:
        print(prefix)
        print("Read geno: "+str(read_geno))
        print("Read snp: "+str(read_snp))
        if read_specific_cols:
            print("Read ind: "+str(True)+" (to indicate selected sample names)")
        else:
            print("Read ind: "+str(read_ind))
        if read_specific_cols:
            print("Reading columns: "+", ".join([str(x) for x in columns]))

    if read_geno:
        if ancestral_allele_encoded_as_0 == True:
            print("[anc:=0] [na:=9]")
        else:
            print("[anc:=1 => anc:=0] [na:=9]")

    if read_specific_cols:
        read_ind = True

    no_ind = True
    if read_ind:
        ind = pd.read_csv(prefix+".ind", sep = "\s+", usecols = [0,2], engine = "c", dtype = str, na_filter = None, header = None).to_numpy()
        no_ind = False
    else:
        ind = None

    if read_specific_cols:
        inds = ind[columns,:]
        df = np.column_stack((columns, inds[:,0], inds[:,1]))
        df = pd.DataFrame(df)
        df.columns = ["0.Index.IndFile", "SampleName.IndFile", "Population"]
        print("_________________________________________________\n")
        print(df)
        print("_________________________________________________")
        del df

    if read_snp:
        snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,2,3], names = ["chrom", "gpos", "ppos"], engine = "c", dtype = {"chrom": str, "gpos": str, "ppos": str}, na_filter = None, header = None, float_precision = 'round_trip', nrows = 1)
        if snp["gpos"][0] == ".":
            snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,3], names = ["chrom", "ppos"], engine = "c", dtype = {"chrom": np.int32, "ppos": np.int32}, na_filter = None, header = None, float_precision = 'round_trip')
            genetic_pos = None
            snp = snp.to_numpy()
        else:
            snp = pd.read_csv(prefix+snp_outfix, sep = "\s+", usecols = [1,2,3], names = ["chrom", "gpos", "ppos"], engine = "c", dtype = {"chrom": np.int32, "gpos": str, "ppos": np.int32}, na_filter = None, header = None, float_precision = 'round_trip')
            genetic_pos = np.array([Decimal(x) for x in snp["gpos"]])
            snp = snp[["chrom", "ppos"]].to_numpy()
    else:
        genetic_pos, snp = None, None

    if read_geno:
        if no_ind:
            if geno_outfix.endswith(".gz"):
                with gzip.open(prefix+geno_outfix, "rt") as f:
                    first_line = f.readline().strip()
            else:
                with open(prefix+geno_outfix) as f:
                    first_line = f.readline().strip()
            nind = len(first_line)
        else:
            nind = ind.shape[0]
        if geno_outfix.endswith(".gz"):
            FGENO = gzip.open(prefix+geno_outfix, "rt")
        else:
            FGENO = open(prefix+geno_outfix, "r")
        if read_specific_cols:
            geno = read_fixedwidth(FGENO, dtype = np.int8, columns = columns)
        else:
            geno = read_fixedwidth(FGENO, dtype = np.int8, columns = list(range(nind)))
        FGENO.close()
        geno = reshape_EIGENSTRAT(geno, isgeno = True)
    else:
        geno = None

    if no_ind == False and read_specific_cols:
        ind = ind[columns,]

    if read_geno and ancestral_allele_encoded_as_0 == False:
        geno = invert_ancestrality(geno)

    # checks
    if read_geno==True and read_snp==True:
        if geno.shape[0] != snp.shape[0]:
            sys.exit("Not equal number of SNPs")
    if read_geno==True and read_ind==True:
        if geno.shape[1] != ind.shape[0]:
            sys.exit("Not equal number of individuals")

    return ind, snp, genetic_pos, geno

# unique() but keep the order of elements (by default unique sorts them)
def uniq(x):
    indexes = np.unique(x, return_index=True)[1]
    return np.array([x[index] for index in sorted(indexes)])

if __name__=="__main__":
    main()
