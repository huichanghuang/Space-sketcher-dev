import os,collections
import argparse
import time
from typing import List
from space_sketcher.__init__ import __root_dir__


class Runpipe:
    def __init__(self, args):
        self.name = args.name
        self.outdir = os.path.abspath(args.outdir)
        self.rna1 = args.rna1
        self.rna2 = args.rna2
        self.oligor1 = args.oligor1
        self.oligor2 = args.oligor2
        self.genomeDir = os.path.abspath(args.genomeDir)
        self.coordfile = args.coordfile
        self.chemistry = args.chemistry        
        self.mapparams = args.mapparams
        self.threads = args.threads        
        self.calling_method = args.calling_method
        self.linker1 = args.linker1
        self.linker2 = args.linker2
        self.sbstart = args.sbstart
        self.maxumi = args.maxumi
        self.minumi = args.minumi
        self.eps = args.eps
        self.min_samples = args.min_samples
        self.minfeatures = args.minfeatures
        self.mincells = args.mincells
        self.ndims = args.ndims
        self.nvariables = args.nvariables
        self.resolution = args.resolution
        self.dev = args.dev
        self.nobam = args.nobam

        self.cbwhitelist = (
            args.cbwhitelist
            if args.cbwhitelist is not None
            else os.path.join(__root_dir__, 
                              f"data/cbwhitelist/{self.chemistry}/cbwhitelist.txt")
        )
        self.sbwhitelist = (
            args.sbwhitelist
            if args.sbwhitelist is not None
            else os.path.join(__root_dir__, 
                              "data/sbwhitelist/sbwhitelist.txt")
        )
        self.reference = (
            args.reference
            if args.reference is not None
            else os.path.basename(os.path.basename(self.genomeDir))
        )

    def runpipe(self):

        ### import lib
        from space_sketcher.tools.utils import (
            str_mkdir, 
            judgeFilexits,
            execute_and_log,
            bin_path,
            rm_temp,
        )
        
        print("test run")
        ### run       
        judgeFilexits(
            self.rna1,
            self.rna2,
            self.oligor1,
            self.oligor2,
            self.genomeDir,
            self.coordfile,
            self.cbwhitelist,
            self.sbwhitelist
            )

        print(bin_path())
        count_cmd = [
            f"{bin_path()}/space-sketcher rna count",
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--rna1 {self.rna1}",
            f"--rna2 {self.rna2}",
            f"--chemistry {self.chemistry}",
            f"--cbwhitelist {self.cbwhitelist}",
            f"--mapparams {self.mapparams}",
            f"--threads {self.threads}",
            f"--genomeDir {self.genomeDir}",
            f"--calling_method {self.calling_method}",
        ]
        count_cmd  = ' '.join(count_cmd) 

        oligo_cmd = [
            f"{bin_path()}/space-sketcher rna oligo",
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--oligor1 {self.oligor1}",
            f"--oligor2 {self.oligor2}",
            f"--sbwhitelist {self.sbwhitelist}",
            f"--cbwhitelist {self.cbwhitelist}",
            f"--chemistry {self.chemistry}",
            f"--linker1 {self.linker1}",
            f"--linker2 {self.linker2}",
            f"--threads {self.threads}",
            f"--coordfile {self.coordfile}",
            f"--sbstart {self.sbstart}",
            f"--maxumi {self.maxumi}",
            f"--minumi {self.minumi}",
            f"--eps {self.eps}",
            f"--min_samples {self.min_samples}",
        ]
        oligo_cmd  = ' '.join(oligo_cmd) 

        analysis_cmd = [
            f'{bin_path()}/space-sketcher rna analysis',
            f"--name {self.name}",
            f"--outdir {self.outdir}",
            f"--minfeatures {self.minfeatures}",
            f"--mincells {self.mincells}",
            f"--ndims {self.ndims}",
            f"--nvariables {self.nvariables}",
            f"--resolution {self.resolution}",
        ]
        analysis_cmd  = ' '.join(analysis_cmd)
        
        
        report_cmd = [
            f'{bin_path()}/space-sketcher rna report',
            f'--name {self.name}',
            f"--outdir {self.outdir}",
            f"--reference {self.reference}",
            f"--chemistry {self.chemistry}",
            f"--dev {self.dev}",
        ]
        report_cmd = ' '.join(report_cmd)
        
        cmdlist = collections.OrderedDict()
        cmdlist['count'] = count_cmd
        cmdlist['oligo'] = oligo_cmd
        cmdlist['analysis'] = analysis_cmd
        cmdlist['report'] = report_cmd

        logdir = os.path.join(self.outdir,self.name)
        str_mkdir(logdir)
        start_time = time.time()
        for pipe, pipecmd in cmdlist.items():
            execute_and_log(pipecmd, pipe, logdir)
        
        end_time = time.time()
        analysis_time = end_time - start_time
        analysis_time_minutes, analysis_time_seconds = divmod(analysis_time, 60)
        analysis_time_hours, analysis_time_minutes = divmod(analysis_time_minutes, 60)

        print(f'\nAnalysis Finished')
        print(f'Elapsed Time: {int(analysis_time_hours)} hours {int(analysis_time_minutes)} minutes {int(analysis_time_seconds)} seconds')

        ###remove bamfile if no bam
        bamfile = os.path.join(self.outdir,self.name,"01.count/Aligned.sortedByCoord.out.bam")
        if self.nobam and os.path.exists(bamfile):
            rm_temp(bamfile)

def run(args):
    Runpipe(args).runpipe()

def helpInfo_run(parser):
    parser.add_argument(
        '-n','--name', 
        metavar='STR',
        help='Sample name.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '-o', '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '-r1', '--rna1',
        metavar='FASTQ',
        help='RNA R1 fastq file, use commas to separate multiple files.', 
        required=True
        )
    parser.add_argument(
        '-r2', '--rna2',
        metavar='FASTQ',
        help='RNA R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with rna1.', 
        required=True
        )
    parser.add_argument(
        '-o1', '--oligor1', 
        metavar='FASTQ',
        help='oligo R1 fastq file, use commas to separate multiple files.', 
        required=True
        )
    parser.add_argument(
        '-o2', '--oligor2', 
        metavar='FASTQ',
        help='oligo R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with cDNAfastq1.', 
        required=True
        )
    parser.add_argument(
        '-g', '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '-cf', '--coordfile', 
        metavar='FILE', 
        type=str,
        required=True,
        help='The coordinate file, in which each spatial barcode with a exact x,y coordinate.'
        )
    parser.add_argument(
        '-cb', '--cbwhitelist',
        metavar='FILE',
        help='Path to the cell barcode whitelist file.'
        )
    parser.add_argument(
        '-sb','--sbwhitelist', 
        metavar='FILE',
        help='Path to the spatial barcode whitelist file.'
        ) 
    parser.add_argument(
        '--chemistry',
        metavar='STR',
        choices=["10X","leader_v1","other"],
        help='Chemistry version, can be "10X", "leader_v1", "other". If set to other, needs to be used with --mapparams. [default: leader_v1].',
        default='leader_v1'
        )
    parser.add_argument(
        '--mapparams',
        metavar='STR',
        help='Additional STAR mapping parameters. must be provide while setting chemistry to other.',
        default='None'
        )
    parser.add_argument(
        '-t', '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '--calling_method',
        metavar='STR',
        choices=["CellRanger2.2","EmptyDrops_CR"],
        help='Cell calling method, Choose from CellRanger2.2 and EmptyDrops_CR, [default: EmptyDrops_CR].', 
        default='EmptyDrops_CR'
        )
    parser.add_argument(
        '--linker1',
        metavar='STR',
        help='linker1 sequence between spatial barcode part1 and part2, [default: TCTTCAGCGTTCCCGAGATCGGACGATCATGGG].',
        default='TCTTCAGCGTTCCCGAGATCGGACGATCATGGG'
        )
    parser.add_argument(
        '--linker2',
        metavar='STR',
        help='linker2 sequence between spatial barcode part2 and part3, [default: CAAGTATGCAGCGCGCTCAAGCACGTGGAT].',
        default='CAAGTATGCAGCGCGCTCAAGCACGTGGAT'
        )
    parser.add_argument(
        '--sbstart', 
        metavar='STRING', 
        type=str,
        help='The spatial barcode start postion, default: auto',
        default='auto'
        )
    parser.add_argument(
        '--maxumi',
        metavar='INT',
        type=int,
        default=5000,
        help='Maximum UMI count threshold for spatial barcodes filtering, [default: 5000]'
        )
    parser.add_argument(
        '--minumi',
        metavar='INT',
        type=int,
        default=2,
        help='Minimum UMI count threshold for spatial barcodes filtering, [default: 2]'
        )
    parser.add_argument(
        '--eps',
        metavar='FLOAT',
        type=float,
        default=150.0,
        help='DBSCAN epsilon parameter (maximum distance between points), [default: 150]'
        )
    parser.add_argument(
        '--min_samples',
        metavar='INT',
        type=int,
        default=6,
        help='DBSCAN min_samples parameter (minimum points to form cluster), [default: 6]'
        )    
    parser.add_argument(
        '--minfeatures',
        type=int, 
        metavar='INT',
        help='Minimum features per cell. default: 5',
        default=5
        )
    parser.add_argument(
        '--mincells',
        type=int, 
        metavar='INT',
        help='Minimum cells per gene. default: 3',
        default=3
        )
    parser.add_argument(
        '--ndims',
        type=int, 
        metavar='INT',
        help='PCA dimensions. default: 30',
        default=30
        )
    parser.add_argument(
        '--nvariables',
        type=int, 
        metavar='INT',
        help='Number of variable genes. default: 2000',
        default=2000
        )
    parser.add_argument(
        '--resolution',
        type=float, 
        metavar='FLOAT',
        help='Clustering resolution. default: 0.5',
        default=0.5
        )
    parser.add_argument(
        '--reference',
        type=str, 
        metavar='STR',
        help='Reference name. default: directory name of genomedir'
        )
    parser.add_argument(
        '--dev',
        type=lambda x: x.lower() == 'true',
        choices=[True, False],
        default=True,
        help='Development mode, selected from [False, True],(default: True)'
    )
    parser.add_argument(
        '--nobam',
        type=lambda x: x.lower() == 'true',
        choices=[True, False],
        default=True,
        help='remove bamfile after finish running, selected from [False, True], default=True.'
    )
    parser.add_argument(
        '--velo',
        type=lambda x: x.lower() == 'true',
        choices=[True, False],
        default=True,
        help='Run STARsolo Velocyto, selected from [False, True], default=True.'
    )
    return parser