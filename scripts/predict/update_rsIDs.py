#!/usr/bin/env python
# update_rsIDs.py
# Laura Colbran, 5-25-17
# updates rsIDs in a dosage file to the ones in PrediXcan model by genomic location
# also filters by whether or not a location has an rsID
#
# Python 3
# Edited 06/28/22 Lin Poyraz to output dosage to specified directory. 
# Edited 07/26/22 Lin Poyraz to get column inputs 
#USAGE:
# python update_rsIDs.py FLAG RS_COL POS_COL REF_COL ALT_COL PATH/TO/PREDIXCAN/SNPS PATH/TO/OUTPUT PATH/TO/FILE/*
# FLAG can be bed or dosage, depending on format of File you want to convert
# RS_COL is an integer, and is the column in the Predixcan file you want to get the rsIDs from. 0-indexed.

import sys
import string
import gzip

#cols for snps in ref file. zero indexed

print ("printing") 
def dosageConvert(file_list,rs_col,ref_file,out_dir,pos_col,ref_col,alt_col):
    for fil in file_list:
        dict = {} #{loc:rsID}
        with gzip.open(fil, 'rt') as f:
            chr = f.readline().split('\t')[0].split('r')[-1] 
            if chr == "": chr = f.readline().split('\t')[0].split('r')[-1]
            f.seek(0)
            with gzip.open(ref_file,'rt') as g: #dict with PrediXcan rsIDs
                for l in g:
                    if l.startswith(chr + "\t") or l.startswith("chr" + chr + "\t"):
                        line = l.rstrip('\n')
                        dict[line.split("\t")[int(pos_col)]] = (line.split("\t")[int(rs_col)],line.split("\t")[int(ref_col)],line.split("\t")[int(alt_col)])
#           output
            print(len(dict.keys()))
            name = out_dir + "/chr" + chr + ".rs_updated.dos.txt"
            of = open(name,"w")
            for line in f:
                if line.startswith("#"):
                    of.write(line)
                    continue
                l = line.split("\t")
                if l[2] in dict:
                    # print(dict[l[2]])
                    alts= dict[l[2]][2].split(",") #possible to have multiple alts in ref file
                    # print(alts)
                    if l[3] == dict[l[2]][1] or l[3] in alts: #check whether predixcan ref allele is present
                        if l[4] == dict[l[2]][1] or l[4] in alts or l[4] == ".": #check whether predixcan alt allele is present
                            l[1]= dict[l[2]][0] #plug in predixcan ids where applicable and write line out
                            of.write("\t".join(l))
                    # sys.exit()
            of.close()

def bedConvert(file_list,rs_col,ref_file):
    dict = {} #{loc:rsID}
    with gzip.open(ref_file,'r') as g: #dict with PrediXcan rsIDs
        for line in g:
            if line.startswith('#'): continue
            dict['\t'.join([line.split('\t')[0],line.split("\t")[1]])] = (line.split("\t")[rs_col])
    for fil in file_list:
        with open(fil, 'r') as inf:
            with open("%s_rsupdated.bed" % (".".join(fil.split(".")[0:-1])),'w') as out:
                out.write("#chr\tpos-1\tpos\trsID\n")
                for line in inf:
                    if line.startswith('#'): continue
                    pos = line.split("\t")[0:3]
                    pos[0] = pos[0].split("r")[1]
                    try:
                        rs = dict["\t".join([pos[0],pos[2]])]
                    except:
                        rs = "."
                    out.write("chr%s\t%s\t%s\t%s\n" % (pos[0],pos[1],pos[2],rs))

def main():
    if sys.argv[1] == "dosage":
        dosageConvert(sys.argv[8:],sys.argv[2],sys.argv[6], sys.argv[7], sys.argv[3], sys.argv[4], sys.argv[5])
    elif sys.argv[1] == "bed":
        bedConvert(sys.argv[8:],sys.argv[2],sys.argv[6])
    else:
        print("ERROR: Specify format as bed or dosage!")

if __name__ == '__main__': main()
