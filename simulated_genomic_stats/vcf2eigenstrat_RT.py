# Convert a vcf file to eigenstrat format
# removes multi-alleleic and indel sites. 
# usage: python vcf2eigenstrat.py -v vcf_file.vcf(.gz) -o out_root
# will generate out_root.[snp,ind,geno].
# removed multiallelic sites and indels
# Deals with haploid cases including mixed haploid/diplod like X as well. 
# -i option is a .ind file to get population names and sex. 

# 190721 : now exports *.geno1 file (haplotypes of the left chromosome)

from __future__ import division
import sys, getopt, gdc, pdb

################################################################################

def parse_options():
    """
    Options are described by the help() function
    """
    options ={ "vcf":None, "out":"out", "ref":None, "indAsPop":False, "indmap":None  }
	
    try:
        opts, args = getopt.getopt(sys.argv[1:], "v:o:r:i:", ["vcf", "out", "ref", "indmap", "indAsPop"])
        print opts, args
    except Exception as err:
        print str(err)
        sys.exit()

    for o, a in opts:
        print o,a
        if o in ["-v","--vcf"]:         options["vcf"] = a
        if o in ["-r","--ref"]:         options["ref"] = a
        if o in ["-i","--ind"]:         options["indmap"] = a
        if o in ["--indAsPop"]:         options["indAsPop"] = True
        elif o in ["-o","--out"]:       options["out"] = a

    print "found options:"
    print options
    return options

################################################################################

def main(options):
    """
    Convert vcf to eigenstrat format (ind, snp and geno files)
    """
    
    print("__/!\__ Beware! This script was tweaked with *reversed* EIGENSTRAT encoding, i.e. 0/0 is encoded as 0 and 1/1 is encoded as 2")
    print("__/!\__ The genetic positions are reported as \".\"")
    
    vcf=gdc.open2(options["vcf"])
    snp, ind, geno, geno1 = [open(options["out"]+x, "w") for x in [".snp", ".ind", ".geno", ".geno1"]]
    removed={"multiallelic":0, "indel":0}
    count=0
    
    if options["indmap"]:
        pop_map={}
        sex_map={}
        ind_map_file=open(options["indmap"], "r")
        for line in ind_map_file:
            bits=line[:-1].split()
            pop_map[bits[0]]=bits[2]
            sex_map[bits[0]]=bits[1]
        ind_map_file.close()
    
    for line in vcf:
        if line[:2]=="##":				  # Comment line
            next
        elif line[:6]=="#CHROM":			  # Header line
            inds=line.split()[9:]
            if options["ref"]:
                ind.write(options["ref"]+"\tU\tREF\n")
            
            if options["indmap"]:
                for indi in inds:
                    ind.write(indi+"\t"+sex_map.get(indi, "U")+"\t"+pop_map.get(indi, "POP")+"\n")
            elif options["indAsPop"]:
                for indi in inds:
                    ind.write(indi+"\tU\t"+indi+"\n")
            else:
                for indi in inds:
                    ind.write(indi+"\tU\tPOP\n")
                   
        else:							  # data
            bits=line.split()
            if "," in bits[4]:
                removed["indel"]+=1
                continue
            if len(bits[3])!=1 or len(bits[4])!=1:
                removed["multiallelic"]+=1
                continue
            else:
                if bits[2]==".":
                    bits[2]=bits[0]+":"+bits[1]
                snp.write(" ".join([bits[2], bits[0], ".", bits[1], "N", "N"])+"\n")
                geno_string=""
                geno1_string=""
                if options["ref"]:
                    geno_string="0"
                    geno1_string="0"
                for gt in bits[9:]:
                    rtgs=decode_gt_string(gt)
                    geno_string+=rtgs[0]
                    geno1_string+=rtgs[1]
                geno.write(geno_string+"\n")
                geno1.write(geno1_string+"\n")
                count+=1

    [f.close for f in [ind, snp, geno, geno1]]

    print "Done. Wrote "+str(count) + " sites"
    print "Excluded " + str(sum(removed.values())) + " sites"
    for key in removed:
        print "Excluded " + str(removed[key]) + " " + key
    return

################################################################################

def decode_gt_string(gt_string):
    """
    Tries to work out the genotype from a vcf genotype entry. 9 for missing [or not in {0,1,2}]
    """
    gt=gt_string.split(":")[0]
    if len(gt)==1:
        if gt=="0":                       # haploid
            return "0","0"
        elif gt=="1":
            return "2","1"
        else:
            return "9","9"
    elif len(gt)==3:
        if gt[0]=="0" and gt[2]=="0":
            return "0","0"
        if gt[0]=="0" and gt[2]=="1":
            return "1","0"
        if gt[0]=="1" and gt[2]=="0":
            return "1","1"
        if gt[0]=="1" and gt[2]=="1":
            return "2","1"
        else:
            return "9","9"

    raise Exception("Unknown genotype: "+gt)
        
################################################################################

if __name__=="__main__":
	options=parse_options()
	main(options)
	
