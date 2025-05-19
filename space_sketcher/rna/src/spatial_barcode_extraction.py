import os
import argparse
import dnaio
import pandas as pd
import editdistance
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import numpy as np
from space_sketcher.tools.utils import add_log

def get_whitelist(_cbwhitelist, chemistry):
    with open(_cbwhitelist, "r") as f:
        tempwhitelist = {line.strip().split()[0] for line in f}
    
    if chemistry == "leader_v1" and len(next(iter(tempwhitelist))) == 10:
        arr = np.array(list(tempwhitelist))
        ##生成两个白名单条形码（barcode）的所有可能组合（笛卡尔积）
        whitelist = set(map(''.join, np.array(np.meshgrid(arr, arr)).T.reshape(-1, 2)))
    else:
        whitelist = tempwhitelist
    return whitelist

# def check_chemistry(chemistry):
#     chemistry = {
#         "10X": (0, 16),
#         "leader_v1": (0, 20)
#     }
#     return chemistry.get(chemistry, (None, None))

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

def auto_detect_sb(seq, linker1, linker2):
    """优化编辑距离计算，提前终止"""
    maxdis=1 ##最大只能设置1, 否则可能找到更多错误的linker,导致sb位置不准确
    def find_linker(seq, linker, maxdis):
        len_linker = len(linker)
        for i in range(len(seq) - len_linker + 1):
            if editdistance.eval(seq[i:i+len_linker], linker) <= maxdis:
                return i
        return -1
    
    linker1_pos = find_linker(seq, linker1, maxdis)
    linker2_pos = find_linker(seq, linker2, maxdis)

    if linker1_pos != -1 and linker2_pos != -1 and linker1_pos < linker2_pos:
        sb1 = seq[linker1_pos-6:linker1_pos]
        sb2 = seq[linker1_pos+len(linker1):linker1_pos+len(linker1)+8]
        sb3 = seq[linker2_pos+len(linker2):linker2_pos+len(linker2)+8]
        sbumi = seq[linker2_pos+len(linker2)+8:linker2_pos+len(linker2)+16]
        return sb1+sb2+sb3, sbumi
    return "", ""

def pos_detect_sb(seq, linker1, linker2, sbstart):
    """向量化位置检测"""
    sbstart = int(sbstart)
    parts = [
        seq[sbstart:sbstart+6],
        seq[sbstart+6+len(linker1):sbstart+6+len(linker1)+8],
        seq[sbstart+6+len(linker1)+8+len(linker2):sbstart+6+len(linker1)+8+len(linker2)+8]
    ]
    sbumi = seq[sbstart+6+len(linker1)+8+len(linker2)+8:sbstart+6+len(linker1)+8+len(linker2)+16]
    return ''.join(parts), sbumi


def process_chunk(records, chemistry, cbwhitelist, sbwhitelist, linker1, linker2, sbstart, cbstart, cblen):
    """处理数据块的并行函数，返回DataFrame"""
    chunk_data = []
    total = cbmatch = 0
    
    for r1, r2 in records:
        total += 1
        cb = r1.sequence[cbstart:cbstart+cblen]
        if cb not in cbwhitelist:
            continue
        cbmatch += 1
        if sbstart == "auto":
            sb, sbumi = auto_detect_sb(r2.sequence, linker1, linker2)
        else:
            sb, sbumi = pos_detect_sb(r2.sequence, linker1, linker2, sbstart)
        if sb in sbwhitelist:
            if chemistry == "leader_v1":
                cb = cb[:10]+"_"+cb[10:] ###匹配RNA文库比对结果
            chunk_data.append([cb, sbumi, sb, 1])  # 初始计数为1
    
    # 创建DataFrame并聚合计数
    if chunk_data:
        df = pd.DataFrame(chunk_data, columns=["Cell_Barcode", "UMI", "Spatial_Barcode", "Read_Count"])
        df = df.groupby(["Cell_Barcode", "UMI", "Spatial_Barcode"]).sum().reset_index()
    else:
        df = pd.DataFrame(columns=["Cell_Barcode", "UMI", "Spatial_Barcode", "Read_Count"])
    
    return df, total, cbmatch

@add_log
def extract_sb_umis(oligoR1, oligoR2, linker1, linker2, sbstart, chemistry, cbwlfile, sbwlfile, outdir, n_jobs=4):
    """并行化处理FASTQ文件"""

    cbstart, cblen = check_chemistry(chemistry)
    cbwhitelist = get_whitelist(cbwlfile, chemistry)
    sbwhitelist = get_whitelist(sbwlfile, "")
      
    with dnaio.open(oligoR1, oligoR2) as reader:
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = []
            chunk = []
            chunk_size = 10000000  # 每个任务处理10M reads
            
            for r1, r2 in tqdm(reader, desc="Preparing chunks"):             
                chunk.append((r1, r2))
                if len(chunk) >= chunk_size:
                    futures.append(executor.submit(
                        process_chunk,
                        chunk.copy(), chemistry, cbwhitelist, sbwhitelist,
                        linker1, linker2, sbstart, cbstart, cblen
                    ))
                    chunk = []
            
            if chunk:  # 处理剩余记录
                futures.append(executor.submit(
                    process_chunk,
                    chunk, chemistry, cbwhitelist, sbwhitelist,
                    linker1, linker2, sbstart, cbstart, cblen
                ))
            
            # 合并DataFrame结果
            result_dfs = []
            totalreads = cbmatch = 0
            
            for future in tqdm(futures, desc="Processing chunks results"):
                chunk_df, chunk_total, chunk_match = future.result()
                totalreads += chunk_total
                cbmatch += chunk_match
                result_dfs.append(chunk_df)
            
    # 合并所有DataFrame
    if result_dfs:
        final_df = pd.concat(result_dfs, ignore_index=True)
        # 再次聚合以合并相同条目的计数
        final_df = final_df.groupby(["Cell_Barcode", "UMI", "Spatial_Barcode"])["Read_Count"].sum().reset_index()
    else:
        final_df = pd.DataFrame(columns=["Cell_Barcode", "UMI", "Spatial_Barcode", "Read_Count"])

    # 保存结果
    spatial_umis_path = os.path.join(outdir, 'spatial_umis.csv.gz')
    final_df.to_csv(spatial_umis_path, index=False, compression='gzip')        

    return spatial_umis_path, totalreads, cbmatch


def Stat_spatial_barcodes(r1fastq, r2fastq, linker1, linker2, 
                          sbstart, chemistry, 
                          cbwhitelist, sbwhitelist, 
                          outdir, n_jobs):

    sb_umis_path, totalreads, cbmatch = extract_sb_umis(r1fastq, r2fastq, linker1, linker2, 
                                                        sbstart, chemistry, 
                                                        cbwhitelist, sbwhitelist, 
                                                        outdir, n_jobs)
                             
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
    parser.add_argument('-j', '--n_jobs', 
        metavar='INT', 
        type=int,
        help='Number of parallel workers, default: 4',
        default=4
        )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    Stat_spatial_barcodes(args.r1fastq, args.r2fastq, 
                          args.linker1, args.linker2, 
                          args.sbstart, args.chemistry, 
                          args.cbwhitelist, args.sbwhitelist, 
                          args.outdir, args.n_jobs)