import os
import sys

class Oligo:
    def __init__(self, args):
        self.oligor1 = args.oligor1
        self.oligor2 = args.oligor2
        self.threads = args.threads
        self.name = args.name
        self.sbwhitelist = args.sbwhitelist
        self.cbwhitelist = args.cbwhitelist
        self.linker1 = args.linker1
        self.linker2 = args.linker2
        self.sbstart = args.sbstart
        self.library = args.library
        self.coordfile = args.coordfile
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        # self.genomeDir = args.genomeDir
        # self.gtf = args.gtf
        # self.chrMT = args.chrMT
        # self.no_introns = args.no_introns
        # self.outunmappedreads = args.outunmappedreads
        # self.end5= args.end5


    # def seqStructure(self) -> str:
    #     """
    #     Determines the sequence structure based on the input fastq files and the provided arguments.

    #     Returns:
    #         str: The sequence structure string.
    #     """
    #     ### import lib
    #     from dnbc4tools.tools.chemistydark import designated_chem_dark
    #     ### 
        
    #     oligoconfiglist = designated_chem_dark(
    #             type = 'oligo',
    #             fastqR1list = self.oligor1, 
    #             fastqR2list = self.oligor2,
    #             chemistry = self.chemistry,
    #             darkreaction = self.darkreaction,
    #             customize = self.customize, 
    #         )
    #     if self.end5:
    #         cDNAconfig = designated_chem_dark(
    #             type = 'rna5p',
    #             fastqR1list = self.cDNAr1, 
    #             fastqR2list = self.cDNAr2,
    #             chemistry = self.chemistry,
    #             darkreaction = self.darkreaction,
    #             customize = self.customize, 
    #         )
    #         oligoconfig = oligoconfiglist
    #     else:
    #         oligoconfig = oligoconfiglist
    #         cDNAconfig = designated_chem_dark(
    #             type = 'rna',
    #             fastqR1list = self.cDNAr1, 
    #             fastqR2list = self.cDNAr2,
    #             chemistry = self.chemistry,
    #             darkreaction = self.darkreaction,
    #             customize = self.customize, 
    #         )
        
    #     return cDNAconfig, oligoconfig

    def run(self):
        print("test oligo")
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.spatial_barcode_extraction import Stat_spatial_barcodes
        from space_sketcher.rna.src.assign_coordinate import assign_coordinate
        
    #     ### import lib
    #     from dnbc4tools.rna.src.staranno import process_libraries
    #     from dnbc4tools.rna.src.oligo_filter import oligo_combine_pl
    #     from dnbc4tools.tools.utils import rm_temp,str_mkdir,judgeFilexits,logging_call, get_formatted_time, change_path
    #     from dnbc4tools.__init__ import __root_dir__
    #     from multiprocessing import Pool


        
        str_mkdir('%s/02.oligo'%self.outdir)
        # str_mkdir('%s/log'%self.outdir)

        ###Extract spatial barcodes
        # judge file exits
        judgeFilexits(
            self.oligor1,
            self.oligor2,
            self.cbwhitelist,
            self.sbwhitelist
            )
        # run Stat_spatial_barcodes
        Stat_spatial_barcodes(self.oligor1, self.oligor2, self.linker1, self.linker2, 
                            self.sbstart, self.library, 
                            self.cbwhitelist, self.sbwhitelist,
                            self.outdir)
        
        ###Assign spatial barcodes coordinate
        # judge file exits
        judgeFilexits(
            self.coordfile,
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
            f"{self.outdir}/02.oligo/spatial_umis.csv.gz",
            f"{self.outdir}/02.oligo/sb_umis_summary-1.temp.csv",
            )        

        #run assign_coordinate
        assign_coordinate(self.coordfile, f"{self.outdir}/02.oligo/spatial_umis.csv.gz",
                          f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
                          f"{self.outdir}/02.oligo/sb_umis_summary-1.temp.csv",
                          self.sbwhitelist, self.outdir)
        
        ###Filter cell barcode by dbscan clustering
        # judge file exits
        judgeFilexits(
            self.coordfile,
            f"{self.outdir}/01.count/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
            f"{self.outdir}/02.oligo/spatial_umis.csv.gz",
            f"{self.outdir}/02.oligo/sb_umis_summary-1.temp.csv",
            )

        
    #     print(f'\n{get_formatted_time()}\n'
    #         f'Conduct quality control for cDNA library barcoding, perform alignment.')
    #     print(f'\n{get_formatted_time()}\n'
    #         f'Perform quality control for oligo library barcodes.')
    #     sys.stdout.flush()

    #     rm_temp(f'{self.outdir}/01.data/Log.final.out',
    #         f'{self.outdir}/01.data/Log.progress.out',
    #         f'{self.outdir}/01.data/Log.out')
        


    #     mission = [[oligo_qc_cmd], [cDNA_star_cmd]]

    #     with Pool(2) as pool:
    #         async_results = []
    #         for i in mission:
    #             async_result = pool.apply_async(logging_call, args=(i,'data',self.outdir,))
    #             async_results.append(async_result)

    #         for result in async_results:
    #             try:
    #                 result.get(timeout=None)
    #             except Exception as e:
    #                 pool.terminate()
    #                 pool.join()
    #                 raise e

    #     if os.path.exists(f'{self.outdir}/01.data/Log.final.out'):
    #         print(f'\n{get_formatted_time()}\n'
    #             f'Annotate gene regions for the aligned BAM.')
    #         logging_call(cDNA_anno_cmd, 'data', self.outdir)
    #     else:
    #         print('\033[0;31;40mUnable to complete cDNA mapping!\033[0m')
    #         raise Exception('Unable to complete cDNA mapping!')

    #     final_sort_cmd = [
    #         f"{__root_dir__}/software/samtools",
    #         f"sort -@ {self.threads}",
    #         f"{self.outdir}/01.data/final.bam",
    #         f"-o {self.outdir}/01.data/final_sorted.bam"
    #     ]

    #     final_sort_cmd_str = " ".join(final_sort_cmd)

    #     logging_call(final_sort_cmd_str, 'data', self.outdir)
        
    #     rm_temp('%s/01.data/Aligned.out.bam'%self.outdir)
    #     rm_temp('%s/01.data/final.bam'%self.outdir)

    #     oligo_combine_pl(f"{__root_dir__}/config/cellbarcode/oligo_type.json",
    #                   f"{self.outdir}/01.data/oligo.reads.fq.gz",
    #                   f"{self.outdir}/01.data",
    #                   f"{__root_dir__}/software",
    #                   f"{self.threads}",
    #                   10,
    #                   logdir=f'{self.outdir}'
    #                   )
    #     rm_temp('%s/01.data/temp'%self.outdir)
    #     rm_temp('%s/01.data/total_reads.xls'%self.outdir)

def oligo(args):
    Oligo(args).run()

def helpInfo_oligo(parser):
    parser.add_argument(
        '--name', 
        metavar='STR',
        help='Sample name.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--cDNAfastq1', 
        metavar='FASTQ',
        help='cDNA R1 fastq file, use commas to separate multiple files.', 
        required=True
        )
    parser.add_argument(
        '--cDNAfastq2', 
        metavar='FASTQ',
        help='cDNA R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with cDNAfastq1.', 
        required=True
        )
    parser.add_argument(
        '--oligofastq1', 
        metavar='FASTQ',
        help='oligo R1 fastq file, use commas to separate multiple files.',
        required=True
        )
    parser.add_argument(
        '--oligofastq2', 
        metavar='FASTQ',
        help='oligo R2 fastq file, use commas to separate multiple files, the files order needs to be consistent with oligofastq1.',
        required=True
        )
    parser.add_argument(
        '--chemistry',
        metavar='STR',
        choices=["scRNAv1HT","scRNAv2HT","scRNAv3HT","auto"],
        help='Chemistry version. Automatic detection is recommended. If setting, needs to be used with --darkreaction, can be "scRNAv1HT", "scRNAv2HT", [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--darkreaction',
        metavar='STR',
        help='Sequencing dark reaction. Automatic detection is recommended. If setting, needs to be used with --chemistry, use comma to separate cDNA and oligo, can be "R1,R1R2", "R1,R1", "unset,unset", [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--customize',
        metavar='STR',
        help='Customize files for whitelist and readstructure in JSON format for cDNA and oligo, use comma to separate cDNA and oligo.'
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--gtf',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--chrMT',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Not include intronic reads in count.'
        )
    parser.add_argument(
        '--outunmappedreads',
        action='store_true',
        help='Output of unmapped reads.'
        )
    parser.add_argument(
        '--end5', 
        action='store_true',
        help='Analyze 5-end single-cell transcriptome data.'
        )
    return parser
