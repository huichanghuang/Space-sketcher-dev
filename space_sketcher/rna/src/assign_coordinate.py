import os
import gzip
import argparse
import pandas as pd
from space_sketcher.tools.utils import csv2dict, add_log
"""
This function assign coordinate to each spatial barcode 
"""

def readlines(_file):
    lines = set()
    if _file.endswith(".gz"):
        with gzip.open(_file, "rt") as f:
            for line in f:
                line = line.strip().split()
                lines.add(line[0])
    else:
        with open(_file, "rt") as f:
            for line in f:
                line = line.strip().split()
                lines.add(line[0])        
    return lines

@add_log
def assign_coordinate(coordfile, sb_umis, true_cbs, summaryfile, sbwlfile, outdir):
   
    true_cells = readlines(true_cbs)    
    ###extract spatial barcode for those calling cell barcodes
    df_umis = pd.read_csv(f"{sb_umis}", compression="gzip")
    
    truecell_umis = df_umis[df_umis["Cell_Barcode"].isin(true_cells)]
    ###spatial barcode截取, 为了匹配puckfile中的barcode
    truecell_umis["SUB_SB"] = df_umis["Spatial_Barcode"].str[:10]+df_umis["Spatial_Barcode"].str[12:18]

    ###summary statistic
    summary = csv2dict(summaryfile)
    summary["Valid Spatial Reads in Cells"] = truecell_umis["Read_Count"].sum()
    summary["Fraction of Valid Spatial Reads in Cells"] = round(summary["Valid Spatial Reads in Cell"]/\
                                                               summary["Valid Spatial Reads"], 4)
    summary["Valid Spatial UMIs in Cell"] = len(truecell_umis)
    summary["Fraction of Valid Spatial UMIs in Cells"] = round(summary["Valid Spatial UMIs in Cell"]/\
                                                               summary["Valid Spatial UMIs"], 4)

    coordf = pd.read_csv(coordfile, header=0)
    coordf.columns = ["SUB_SB", "x", "y"]
    summary["Total Spatial Barcodes with Location in Chip"] = len(coordf)
    
    sbwhitelist = readlines(sbwlfile)
    sbwhitelist_modified = list(map(lambda x: x[:10]+x[12:18], sbwhitelist))
    filtered_coordf = coordf[coordf["SUB_SB"].isin(sbwhitelist_modified)]
    filtered_coordf_dedup = filtered_coordf.drop_duplicates(subset='SUB_SB') ##根据Spatial_Barcode去重
    summary["Unique Valid Spatial Barcodes with Location in Chip"] = len(filtered_coordf_dedup) ##unique and in whitelist
    summary["Fraction of Unique Valid Spatial Barcodes with Location in Chip"] = round(summary["Unique Valid Spatial Barcodes with Location in Chip"]/\
                                                                                       summary["Total Spatial Barcodes with Location in Chip"], 4)

    ###merge coordfile and sb_umis
    mergedf = pd.merge(truecell_umis, filtered_coordf_dedup, on="SUB_SB", how="inner")
    summary["Spatial Reads in Cells with Location in Chip"] = mergedf["Read_Count"].sum()
    summary["Fraction of Spatial Reads in Cells with Location in Chip"] = round(summary["Spatial Reads in Cell with Location in Chip"]/\
                                                                               summary["Valid Spatial Reads in Cell"],4)

    mergedf_dropped = mergedf.drop("Read_Count", axis=1)
    mergedf_dropped["UMI_count"] = mergedf_dropped.groupby(['Cell_Barcode', 'SUB_SB'])['UMI'].transform('nunique')
    df_sorted = mergedf_dropped.sort_values(by=['Cell_Barcode', 'UMI_count'], ascending=[True, False])
    df_sorted = df_sorted.drop("UMI", axis=1)
    df_sorted_dedup = df_sorted.drop_duplicates()
    df_sorted_dedup.columns = ["cb", "sb", "subsb", "xcoord", "ycoord", "umi_count"]

    outfile = os.path.join(outdir, "cb_sb_coord.txt")
    df_sorted_dedup.to_csv(outfile, index=False, header=True, sep="\t")

    outsummary = os.path.join(outdir, "sb_library_summary.csv")
    with open(outsummary, "wt") as outf:
        for k, v in summary.items():
            print(f"{k},{v}", file=outf)
    
    ###prepare barcode rank file for spatial knee plot
    cell_spatial_umi = df_umis.groupby(['Cell_Barcode', 'Spatial_Barcode'])['UMI'].nunique().reset_index()
    cell_spatial_umi = cell_spatial_umi.rename(columns={'UMI': 'spatial_UMI_count'})
    cell_spatial_umi = cell_spatial_umi.groupby('Cell_Barcode')['spatial_UMI_count'].sum().reset_index()
    cell_spatial_umi = cell_spatial_umi.rename(columns={'spatial_UMI_count': 'UMI'})
    cell_spatial_umi['is_cell_barcode'] = cell_spatial_umi['Cell_Barcode'].isin(true_cells).astype(int)
    cell_spatial_umi = cell_spatial_umi.sort_values('UMI', ascending=False)
    cell_spatial_umi['rank'] = range(1, len(cell_spatial_umi)+1)
    final_output = cell_spatial_umi.rename(columns={'Cell_Barcode': 'barcode'})
    final_output = final_output[['barcode', 'UMI', 'is_cell_barcode', 'rank']]
    final_output.to_csv(os.path.join(outdir, 'cell_sb_umi.rank.txt'), sep='\t', index=False)



def parse_args():
    parser = argparse.ArgumentParser(description='Assign coordinate to each spatial barcode.')
    parser.add_argument('-c', '--coordfile', 
        metavar='FILE', 
        type=str,
        help='The coordinate file, in which each spatial barcode with a exact x,y coordinate.'
        )
    parser.add_argument('-i', '--infile', 
        metavar='FILE', 
        type=str,
        help='Input file generated by spatial_barcode_extraction, contain Cell_Barcode,UMI,Spatial_Barcode,Read_Count.'
        )
    parser.add_argument('-sm', '--summary', 
        metavar='FILE', 
        type=str,
        help='Temporary summary file generated by spatial_barcode_extraction.'
        )
    parser.add_argument('-sb', '--sbwhitelist', 
        metavar='FILE', 
        type=str,
        help='The spatial barcode whitelist files.'
        )
    parser.add_argument('-tc', '--truecells', 
        metavar='FILE', 
        type=str,
        help='The filtered cell barcode files generated by STARsolo mapping.'
        )
    parser.add_argument('-o', '--outdir', 
        metavar='PATH', 
        type=str,
        help='The output directory'
        )
    

if __name__=='__main__':
    args = parse_args()

    assign_coordinate(args.coordfile, args.infile, args.truecells, 
                      args.summary, args.sbwhitelist, args.outdir)