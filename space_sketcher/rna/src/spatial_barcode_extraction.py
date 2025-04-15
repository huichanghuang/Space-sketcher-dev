import os
import gzip
import argparse
from collections import defaultdict
import dnaio
from itertools import product
import pandas as pd
from space_sketcher.tools.utils import add_log
import editdistance
"""
This function extract the umi and readcount of spatial barcode from oligo fastq,
and calculate the statistic of spatial barcode matching and saturation
"""

def get_whitelist(_cbwhitelist, chemistry):
    tempwhitelist = set()
    with open(_cbwhitelist, "r") as f:
        for line in f:
            line = line.strip().split()
            tempwhitelist.add(line[0])
    
    if chemistry == "leader_v1" and len(line[0]) == 10:
        whitelist = set()
        for cb1,cb2 in product(tempwhitelist, tempwhitelist):
            whitelist.add(cb1+cb2)
    else:
        whitelist = tempwhitelist

    return whitelist

def check_chemistry(chemistry):
    if chemistry == "10X":
        cbstart = 0
        cblen = 16
    elif chemistry == "leader_v1":
        cbstart = 0
        cblen = 20
    else:
        print("Please use the correct chemistry input, 10X or leader_v1.")
        exit
    return cbstart, cblen

def auto_detect_sb(seq, linker1, linker2, maxdis=2):
    """
    detect sb by linker position and match, it may takes a very long time.
    """
    linker1_pos = -1
    linker2_pos = -1
    #find linker1
    for i in range(len(seq) - len(linker1) + 1):
        if editdistance.eval(seq[i:i+len(linker1)], linker1) <= maxdis:
            linker1_pos = i
            l1match += 1
            break
    
    #find linker2
    for i in range(len(seq) - len(linker2) + 1):
        if editdistance.eval(seq[i:i+len(linker2)], linker2) <= maxdis:
            linker2_pos = i
            l2match += 1
            break

    if linker1_pos != -1 and linker2_pos != -1 and linker1_pos < linker2_pos:
        sb1 = seq[linker1_pos-6:linker1_pos]
        sb2 = seq[linker1_pos+len(linker1):linker1_pos+len(linker1)+8]
        sb3 = seq[linker2_pos+len(linker2):linker2_pos+len(linker2)+8]
        sbumi = seq[linker2_pos+len(linker2)+8:linker2_pos+len(linker2)+8+8]
    else:
        sb1, sb2, sb3, sbumi = "", "", "", ""

    sb = sb1+sb2+sb3
    return sb,sbumi

def pos_detect_sb(seq, linker1, linker2, sbstart):
    """
    extract sb base on position
    """
    sb1 = seq[sbstart:sbstart+6]
    sb2 = seq[sbstart+6+len(linker1):
              sbstart+6+len(linker1)+8]
    
    sb3 = seq[sbstart+6+len(linker1)+8+len(linker2):
              sbstart+6+len(linker1)+8+len(linker2)+8]
    
    sbumi = seq[sbstart+6+len(linker1)+8+len(linker2)+8:
                sbstart+6+len(linker1)+8+len(linker2)+8+8]
    sb = sb1+sb2+sb3

    return sb, sbumi

@add_log
def extract_sb_umis(oligoR1, oligoR2, linker1, linker2, sbstart, chemistry, cbwlfile, sbwlfile, outdir):
    
    cbstart, cblen = check_chemistry(chemistry)
    cbwhitelist = get_whitelist(cbwlfile, chemistry)
    sbwhitelist = get_whitelist(sbwlfile, "")

    totalreads, cbmatch = 0, 0
    spatial_umis = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    if sbstart == "auto":
        with dnaio.open(oligoR1, oligoR2) as reader:
            for r1, r2 in reader:
                totalreads += 1
                r1seq = r1.sequence
                r2seq = r2.sequence
                cb = r1seq[cbstart:cblen]
                if cb not in cbwhitelist: continue
                if chemistry == "leader_v1":
                    cb = cb[:10]+"_"+cb[10:20] ##to match the RNA cell barcode 
                cbmatch += 1
                sb, sbumi = auto_detect_sb(r2seq, linker1, linker2)
                if sb in sbwhitelist:
                    spatial_umis[cb][sbumi][sb] += 1
    else:
        sbstart = int(sbstart)
        with dnaio.open(oligoR1, oligoR2) as reader:
            for r1, r2 in reader:
                totalreads += 1
                r1seq = r1.sequence
                r2seq = r2.sequence
                cb = r1seq[cbstart:cblen]
                if cb not in cbwhitelist: continue
                if chemistry == "leader_v1":
                    cb = cb[:10]+"_"+cb[10:20] ##to match the RNA cell barcode 
                cbmatch += 1
                sb, sbumi = pos_detect_sb(r2seq, linker1, linker2, sbstart)
                if sb in sbwhitelist:
                    spatial_umis[cb][sbumi][sb] += 1

    spatial_umis_path = os.path.join(outdir, 'spatial_umis.csv.gz')
    with gzip.open(spatial_umis_path, 'wt') as f_umis:
        f_umis.write(",".join(["Cell_Barcode", "UMI", "Spatial_Barcode", "Read_Count"])+"\n")
        unique_sbs = set()
        for cb, sp_umis in spatial_umis.items():
            line = ''.join([",".join([cb, umi, sb, str(count)]) + "\n"
                            for umi, sb_count in sp_umis.items()
                            for sb, count in sb_count.items()
                            if (unique_sbs.add(sb) or True)])
            f_umis.write(line)

    return spatial_umis_path, totalreads, cbmatch


def Stat_spatial_barcodes(r1fastq, r2fastq, linker1, linker2, 
                          sbstart, chemistry, 
                          cbwhitelist, sbwhitelist, 
                          outdir):

    sb_umis_path, totalreads, cbmatch = extract_sb_umis(r1fastq, r2fastq, linker1, linker2, 
                                                        sbstart, chemistry, 
                                                        cbwhitelist, sbwhitelist, 
                                                        outdir)
                             
    summary = {"Total Spatial Reads": totalreads}
    summary["Spatial Reads with Valid Cellbarcode"] = cbmatch
    df_umis = pd.read_csv(f"{sb_umis_path}", compression="gzip")
    ### calculate saturation
    summary["Valid Spatial Reads"] = df_umis["Read_Count"].sum()
    summary["Fraction of Valid Spatial Reads"] = round(summary["Valid Spatial Reads"]/\
                                                       summary["Total Spatial Reads"], 4)
    summary["Valid Spatial UMIs"] = len(df_umis)
    summary["Spatial Barcode Saturation"] = round(1-(summary["Valid Spatial UMIs"]/summary["Valid Spatial Reads"]), 4)

    outstat = os.path.join(outdir, "sb_library_summary.temp.csv")
    with open(outstat, "wt") as outf:
        for k, v in summary.items():
            print(f"{k},{v}", file=outf)
    

def parse_args():
    parser = argparse.ArgumentParser(description='Counting the spatial umi saturation')
    parser.add_argument('-r1', '--r1fastq', 
        metavar='FILE', 
        type=str,
        help='The R1 fastq file, contain cell barcode'
        )
    parser.add_argument('-r2', '--r2fastq', 
        metavar='FILE', 
        type=str,
        help='The R2 fastq file, contain spatial barcode and linkers'
        )
    parser.add_argument('-sb', '--sbwhitelist', 
        metavar='FILE', 
        type=str,
        help='The spatial barcode whitelist files'
        )
    parser.add_argument('-cb', '--cbwhitelist', 
        metavar='FILE', 
        type=str,
        help='The cell barcode whitelist files'
        )
    parser.add_argument('-c', '--chemistry',
        metavar='STR',
        type=str,
        choices=['10X', 'leader_v1'],
        default='leader_v1',
        help='Chemistry version: 10X or leader_v1, [default: leader_v1]'
        )
    parser.add_argument('-o', '--outdir', 
        metavar='PATH', 
        type=str,
        help='The output directory'
        )
    parser.add_argument('-l1', '--linker1', 
        metavar='STRING', 
        type=str,
        help='The linker1 sequence',
        default='TCTTCAGCGTTCCCGAGATCGGACGATCATGGG'
        )
    parser.add_argument('-l2', '--linker2', 
        metavar='STRING', 
        type=str,
        help='The linker2 sequence',
        default='CAAGTATGCAGCGCGCTCAAGCACGTGGAT'
        )
    parser.add_argument('-st', '--sbstart', 
        metavar='STRING', 
        type=str,
        help='The spatial barcode start postion, default: auto',
        default='auto'
        )

if __name__=='__main__':
    args = parse_args()

    Stat_spatial_barcodes(args.r1fastq, args.r2fastq, args.linker1, args.linker2, 
                        args.sbstart,args.chemistry, 
                        args.cbwhitelist, args.sbwhitelist, 
                        args.outdir)