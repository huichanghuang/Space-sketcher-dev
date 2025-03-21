import os,gzip, argparse
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

def get_whitelist(_cbwhitelist, _ltype):
    tempwhitelist = set()
    with open(_cbwhitelist, "r") as f:
        for line in f:
            line = line.strip().split()
            tempwhitelist.add(line[0])
    
    if _ltype == "leader_v1" and len(line[0]) == 10:
        whitelist = set()
        for cb1,cb2 in product(tempwhitelist, tempwhitelist):
            whitelist.add(cb1+cb2)
    else:
        whitelist = tempwhitelist

    return whitelist


def check_chemistry(_ltype):
    if _ltype == "10X":
        cbstart = 0
        cblen = 16
    elif _ltype == "leader_v1":
        cbstart = 0
        cblen = 20
    else:
        print("Please use the correct LibraryType input, 10X or leader_v1.")
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
def extract_sb_umis(oligoR1, oligoR2, linker1, linker2, sbstart, _ltype, cbwlfile, sbwlfile, spatial_dir, samplename):
    
    cbstart, cblen = check_chemistry(_ltype)
    cbwhitelist = get_whitelist(cbwlfile, _ltype)
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
                cbmatch += 1
                sb, sbumi = pos_detect_sb(r2seq, linker1, linker2, sbstart)
                if sb in sbwhitelist:
                    spatial_umis[cb][sbumi][sb] += 1

    spatial_umis_path = os.path.join(spatial_dir, f'{samplename}.spatial_umis.csv.gz')
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

@add_log
def get_summary(summary, true_cbs, sb_umis, _ltype, spatial_dir, samplename):
    true_cells = get_whitelist(true_cbs, "")
    if _ltype == "leader_v1":
        ###structure of cell barcode in leader_v1 RNA matrix: AAAAAAAAAA_BBBBBBBBBB, must remove the "_"
        true_cells = list(map(lambda x: x[:10]+x[-10:], true_cells))

    df_umis = pd.read_csv(f"{sb_umis}", compression="gzip")

    ### calculate saturation
    summary["Valid_Spatial_Reads"] = df_umis["Read_Count"].sum()
    summary["Valid_Spatial_Reads_ratio"] = round(summary["Valid_Spatial_Reads"]/summary["Total_Reads"], 3)
    summary["Total_Spatial_UMIs"] = len(df_umis)
    summary["Spatial_Barcode_Saturation"] = round(1-(summary["Total_Spatial_UMIs"]/summary["Valid_Spatial_Reads"]), 3)

    truecell_umis = df_umis[df_umis["Cell_Barcode"].isin(true_cells)]
    summary["Valid_Spatial_Reads_in_cell"] = truecell_umis["Read_Count"].sum()
    summary["Total_Spatial_UMIs_in_cell"] = len(truecell_umis)
    summary["Total_Spatial_UMIs_in_cell_ratio"] = round(len(truecell_umis)/len(df_umis), 3)
    summary["Spatial_Barcode_Saturation_in_cell"] = round(1-(summary["Total_Spatial_UMIs_in_cell"]/summary["Valid_Spatial_Reads_in_cell"]),3)

    outfile = os.path.join(spatial_dir, f'{samplename}.trucells-spatial_umis.csv.gz')
    truecell_umis.to_csv(outfile, compression="gzip")

    return summary

def Stat_spatial_barcodes(r1fastq, r2fastq, linker1, linker2, 
                          sbstart, LibraryType, 
                          cellbarcode, spatialbarcode, 
                          outdir, sample, truecell):

    sb_umis_path, totalreads, cbmatch = extract_sb_umis(r1fastq, r2fastq, linker1, linker2, 
                                                        sbstart, LibraryType, 
                                                        cellbarcode, spatialbarcode, 
                                                        outdir, sample)
                             
    tempdict = {"Total_Reads": totalreads}
    tempdict = {"Cell_barcode_matched_reads": cbmatch}

    summary = get_summary(tempdict, truecell, sb_umis_path, LibraryType, outdir, sample)

    outstat = os.path.join(outdir, f"{sample}.sb_umis_saturation.txt")
    with open(outstat, "wt") as outf:
        for k, v in summary.items():
            print(f"{k}: {v}", file=outf)
    

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
    parser.add_argument('-sb', '--spatialbarcode', 
        metavar='FILE', 
        type=str,
        help='The spatial barcode whitelist files'
        )
    parser.add_argument('-cb', '--cellbarcode', 
        metavar='FILE', 
        type=str,
        help='The cell barcode whitelist files'
        )
    parser.add_argument('-t', '--truecell', 
        metavar='FILE', 
        type=str,
        help='The filtered cell list called by starsolo'
        )
    parser.add_argument('-l', '--LibraryType', 
        metavar='STRING', 
        type=str,
        choices=['10X','leader_v1'],
        help='LibraryType, can be: 10X, leader_v1, default: 10X',
        default='10X'
        )
    parser.add_argument('-o', '--outdir', 
        metavar='PATH', 
        type=str,
        help='The output directory'
        )
    parser.add_argument('-s', '--sample', 
        metavar='STRING', 
        type=str,
        help='Sample name'
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
                        args.sbstart,args.LibraryType, 
                        args.cellbarcode, args.spatialbarcode, 
                        args.outdir, args.sample, args.truecell)