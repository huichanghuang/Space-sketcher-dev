import os
import sys
from loguru import logger
from pathlib import Path

class Oligo:
    def __init__(self, args):
        self.name = args.name
        self.outdir = os.path.abspath(os.path.join(args.outdir, args.name))
        self.oligor1 = args.oligor1
        self.oligor2 = args.oligor2
        self.sbwhitelist = args.sbwhitelist
        self.cbwhitelist = args.cbwhitelist
        self.chemistry = args.chemistry        
        self.linker1 = args.linker1
        self.linker2 = args.linker2
        self.threads = args.threads
        self.coordfile = args.coordfile
        self.sbstart = args.sbstart
        self.maxumi = args.maxumi
        self.minumi = args.minumi
        self.eps = args.eps
        self.min_samples = args.min_samples
               

    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits, get_formatted_time
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.spatial_barcode_extraction import Stat_spatial_barcodes
        from space_sketcher.rna.src.assign_coordinate import assign_coordinate
        from space_sketcher.rna.src.dbscan_filter import dbscan_filter
        print("test oligo")
        
        oligodir = os.path.join(self.outdir, "02.oligo")
        str_mkdir(oligodir)

        ###Extract spatial barcodes
        # judge file exits
        judgeFilexits(
            self.oligor1,
            self.oligor2,
            self.cbwhitelist,
            self.sbwhitelist
            )
        # run Stat_spatial_barcodes
        if not os.path.exists(f"{oligodir}/spatial_umis.csv.gz"):
            logger.info(f'\n{get_formatted_time()}\n'
                f'Extracting spatial barcode information.')
            Stat_spatial_barcodes(self.oligor1, self.oligor2, 
                                self.linker1, self.linker2, 
                                self.sbstart, self.chemistry, 
                                self.cbwhitelist, self.sbwhitelist,
                                oligodir, self.threads)
        else:
            logger.info(f"{oligodir}/spatial_umis.csv.gz already exits, skip Extracting spatial barcode information.")
            
        ###Assign spatial barcodes coordinate
        # judge file exits
        judgeFilexits(
            self.coordfile,
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
            f"{oligodir}/spatial_umis.csv.gz",
            f"{oligodir}/sb_library_summary.temp.csv",
            )

        #run assign_coordinate
        logger.info(f'\n{get_formatted_time()}\n'
            f'Assigning coordinate for each spatial barcode.')
        assign_coordinate(self.coordfile, f"{oligodir}/spatial_umis.csv.gz",
                          f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
                          f"{oligodir}/sb_library_summary.temp.csv",
                          self.sbwhitelist, oligodir)
        
        ###Filter cell barcode by dbscan clustering
        # judge file exits
        judgeFilexits(
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered",
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/CellReads.stats",
            f"{oligodir}/cb_sb_coord.txt",
            )

        #run dbscan_filter
        logger.info(f'\n{get_formatted_time()}\n'
            f'Performing dbscan filtering.')
        dbscan_filter(f"{oligodir}/cb_sb_coord.txt",
                      oligodir, 
                      self.maxumi, self.minumi, 
                      f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered", 
                      f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/CellReads.stats", 
                      self.eps, self.min_samples, self.threads)        

        (Path(self.outdir) / ".oligo.done").touch()

def oligo(args):
    Oligo(args).run()

def helpInfo_oligo(parser):
    parser.add_argument(
        '-n', '--name', 
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
        '-sb','--sbwhitelist', 
        metavar='FILE',
        help='Path to the spatial barcode whitelist file.',
        required=True
        )
    parser.add_argument(
        '-cb', '--cbwhitelist',
        metavar='FILE',
        help='Path to the cell barcode whitelist file.',
        required=True
        )
    parser.add_argument(
        '-c', '--chemistry',
        metavar='STR',
        choices=["10X","leader_v1","other"],
        help='Chemistry version, can be "10X", "leader_v1", "other". If set to other, needs to be used with --mapparams. [default: leader_v1].',
        default='leader_v1'
        )
    parser.add_argument(
        '-l1', '--linker1',
        metavar='STR',
        help='linker1 sequence between spatial barcode part1 and part2, [default: TCTTCAGCGTTCCCGAGATCGGACGATCATGGG].',
        default='TCTTCAGCGTTCCCGAGATCGGACGATCATGGG'
        )
    parser.add_argument(
        '-l2', '--linker2',
        metavar='STR',
        help='linker2 sequence between spatial barcode part2 and part3, [default: CAAGTATGCAGCGCGCTCAAGCACGTGGAT].',
        default='CAAGTATGCAGCGCGCTCAAGCACGTGGAT'
        )
    parser.add_argument(
        '-t', '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '-cf', '--coordfile', 
        metavar='FILE', 
        type=str,
        help='The coordinate file, in which each spatial barcode with a exact x,y coordinate.'
        )
    parser.add_argument(
        '-st', '--sbstart', 
        metavar='STRING', 
        type=str,
        help='The spatial barcode start postion, default: auto',
        default='auto'
        )
    parser.add_argument(
        '-max', '--maxumi',
        metavar='INT',
        type=int,
        default=5000,
        help='Maximum UMI count threshold for spatial barcodes filtering, [default: 5000]'
        )
    parser.add_argument(
        '-min','--minumi',
        metavar='INT',
        type=int,
        default=2,
        help='Minimum UMI count threshold for spatial barcodes filtering, [default: 2]'
        )
    parser.add_argument(
        '-e', '--eps',
        metavar='FLOAT',
        type=float,
        default=150.0,
        help='DBSCAN epsilon parameter (maximum distance between points), [default: 150]'
        )
    parser.add_argument(
        '-ms', '--min_samples',
        metavar='INT',
        type=int,
        default=6,
        help='DBSCAN min_samples parameter (minimum points to form cluster), [default: 6]'
        )
    return parser
