import os
import subprocess
from loguru import logger
from pathlib import Path

class Count:
    def __init__(self,args):
        self.rna1 = args.rna1
        self.rna2 = args.rna2
        self.threads = args.threads
        self.name = args.name
        self.chemistry = args.chemistry
        self.whitelist = args.whitelist
        self.mapparams = args.mapparams
        self.genomeDir = args.genomeDir
        self.gtf = args.gtf
        self.chrMT = args.chrMT
        self.no_introns = args.no_introns
        self.outunmappedreads = args.outunmappedreads
        self.end5= args.end5
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        self.minumi = args.minumi
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))

    def Prepare_mapping_params(self) -> str:
        """
        Check the chemistry and prepare mapping parameters.
        Returns:
            str: The mapping parameters.
        """
        ###load function
        from space_sketcher.tools.utils import gunzip

        mapping_params = ""
        if self.chemistry == "10X":
            mapping_params += "--soloType CB_UMI_Simple "
            mapping_params += "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 "
            mapping_params += "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloBarcodeReadLength 0 "
        elif self.chemistry == "leader_v1":
            mapping_params += "--soloType CB_UMI_Complex "
            mapping_params += "--soloCBposition 0_0_0_9 0_10_0_19 --soloUMIposition 0_20_0_29 "
            mapping_params += "--soloCBmatchWLtype EditDist_2 "
        elif self.chemistry == "other":
            if self.mapparams == "":
                print("Please check if chemistry and mapparams provided in proper way!")
        else:
            print("Not available chemistry")

        mapping_params += self.mapparams
        ###to avoid STAR memery error
        if self.whitelist.endswith(".gz"):
            self.whitelist = gunzip(self.whitelist)

        ###add whitelist
        if "--soloType CB_UMI_Simple " in mapping_params:
            mapping_params += f"--soloCBwhitelist {self.whitelist} "
        elif "--soloType CB_UMI_Complex " in mapping_params:
            mapping_params += f"--soloCBwhitelist {self.whitelist} {self.whitelist} "
        else:
            ##log
            print("Please check if chemistry and mapparams provided in proper way!")

        return mapping_params
        
    def run(self):
        print("test count")
        ### import lib
        from space_sketcher.tools.utils import str_mkdir,logging_call,judgeFilexits, get_formatted_time,rm_temp, create_index,gunzip
        # from dnbc4tools.rna.src.singlecell_summary import cut_umi,generateCellSummary
        # from dnbc4tools.tools.cal_saturation import sub_sample_cDNA_rna
        # from dnbc4tools.tools.cell_calling import cell_calling
        # from dnbc4tools.tools.plot_draw import merge_graph
        # from dnbc4tools.tools.combineBeads import similarity_droplet_file,barcodeTranslatefile
        from space_sketcher.__init__ import __root_dir__

        ### run
        judgeFilexits(
            self.rna1,
            self.rna2,
            self.whitelist,
            self.genomeDir,
            self.gtf
            )
           
        str_mkdir('%s/01.count'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        str_mkdir('%s/log/.temp'%self.outdir)
        os.environ[ 'MPLCONFIGDIR' ] = '%s/log/.temp'%self.outdir
        os.environ[ 'NUMBA_CACHE_DIR' ] = '%s/log/.temp'%self.outdir
        
        star_version = subprocess.check_output(f"{__root_dir__}/software/STAR --version", shell=True)
        logger.info(f"STAR 版本号：{star_version.decode('utf8')}")

        mapping_pars = self.Prepare_mapping_params()
        STAR_cmd = (
            f"{__root_dir__}/software/STAR "
            f"--runMode alignReads "
            "--soloFeatures GeneFull_Ex50pAS "
            "--quantMode GeneCounts "
            "--soloCellFilter  EmptyDrops_CR"
            "--outFilterScoreMin 30 "
            "--soloStrand Unstranded "
            "--readFilesCommand zcat "
            "--soloCellReadStats Standard "
            "--soloMultiMappers EM "
            "--soloUMIdedup 1MM_CR "
            "--soloUMIfiltering MultiGeneUMI_CR "
            "--clipAdapterType CellRanger4 "
            "--outSAMtype BAM SortedByCoordinate "
            "--outSAMattributes CR CY UR UY NH HI nM AS GX GN gx gn CB UB sS sQ sM "
            f"{mapping_pars} "
            f"--readFilesIn {self.rna2} {self.rna1} "
            f"--genomeDir {self.genomeDir} "
            f"--outFileNamePrefix {self.gtf}/ "
            f"--runThreadN {self.threads} "
            "--outBAMsortingThreadN 6 "
        )

        subprocess.check_call(STAR_cmd, shell=True)
        ###change output directory permissions
        chmod_cmd = f"chmod -R a+r {self.outdir}/Solo.out && chmod a+x $(find {self.outdir}/Solo.out -type d)"
        subprocess.check_call(chmod_cmd, shell=True)
        
        # 对 STARsolo 文件进行压缩
        cmd = f"gzip {self.outdir}/Solo.out/GeneFull_Ex50pAS/raw/* && gzip {self.outdir}/Solo.out/GeneFull_Ex50pAS/filtered/*"
        subprocess.check_call(cmd, shell=True)

        (Path(self.outdir) / ".STAR.done").touch()
        # ## filter oligo
        # print(f'\n{get_formatted_time()}\n'
        #     f'Calculating bead similarity and merging beads within the same droplet.')
        
        # cut_umi(
        #     f"{self.outdir}/01.data/beads_stat.txt",100,
        #     '%s/02.count'%self.outdir
        #     )

        # similiarBeads_cmd = [
        #     f"{__root_dir__}/software/similarity",
        #     f"-n {self.threads}",
        #     f"{self.name}",
        #     f"{self.outdir}/01.data/CB_UB_count.txt",
        #     f"{self.outdir}/02.count/beads.barcodes.umi100.txt",
        #     f"{__root_dir__}/config/cellbarcode/oligo_type.txt",
        #     f"{self.outdir}/02.count/similarity.all.csv",
        #     f"{self.outdir}/02.count/similarity.droplet.csv",
        #     f"{self.outdir}/02.count/similarity.dropletfiltered.csv"
        # ]
        # similiarBeads_cmd_str = " ".join(similiarBeads_cmd)  
        # logging_call(similiarBeads_cmd_str,'count',self.outdir)

        # ### merge beads list
        # similarity_droplet_file('%s/02.count/similarity.droplet.csv'%self.outdir,
        #                         '%s/02.count/beads.barcodes.umi100.txt'%self.outdir,
        #                         '%s/02.count/combined.list.tsv'%self.outdir,
        #                         0.4,
        #                         1,
        #                         logdir=f'{self.outdir}')

        # barcodeTranslatefile(
        #     f"{self.outdir}/02.count/combined.list.tsv", 
        #     f"{self.outdir}/01.data/beads_stat.txt", 
        #     f"{self.outdir}/02.count/barcodeTranslate.txt",
        #     f"{self.outdir}/02.count/cell.id",
        #     f"{self.outdir}/02.count/barcodeTranslate.hex.txt",
        #     logdir=f'{self.outdir}'
        #     )
    
        # ### add DB tag for bam
        # # print(f'\n{get_formatted_time()}\t'
        # #     f'Generating anno decon sorted bam.')
        # tagAdd_cmd = [
        #     f"{__root_dir__}/software/tagAdd",
        #     f"-n {self.threads}",
        #     f"-bam {self.outdir}/01.data/final_sorted.bam",
        #     f"-file {self.outdir}/02.count/barcodeTranslate.hex.txt",
        #     f"-out {self.outdir}/02.count/anno_decon_sorted.bam",
        #     "-tag_check CB:Z:",
        #     "-tag_add DB:Z:"
        # ]
        # tagAdd_cmd_str = " ".join(tagAdd_cmd)
        # logging_call(tagAdd_cmd_str,'count',self.outdir)

        # ### get bam index
        # create_index(self.threads,'%s/02.count/anno_decon_sorted.bam'%self.outdir,self.outdir)

        # print(f'\n{get_formatted_time()}\n'
        #     f'Generating the raw expression matrix.')
        # str_mkdir('%s/02.count/raw_matrix'%self.outdir)
        # PISA_countRaw_cmd = [
        #     f"{__root_dir__}/software/PISA",
        #     "count",
        #     "-one-hit",
        #     f"-@ {self.threads}",
        #     "-cb DB",
        #     "-anno-tag GN",
        #     "-umi UB",
        #     f"-list {self.outdir}/02.count/cell.id",
        #     f"-outdir {self.outdir}/02.count/raw_matrix",
        #     f"{self.outdir}/02.count/anno_decon_sorted.bam"
        # ]
        # PISA_countRaw_cmd_str = " ".join(PISA_countRaw_cmd)
        # logging_call(PISA_countRaw_cmd_str,'count',self.outdir)


        # ## cell calling using DropletUtils
        # if self.forcecells:
        #     cell_bc, count_num = cell_calling(
        #         f"{self.outdir}/02.count/raw_matrix/", 
        #         force_cell_num = int(self.forcecells), 
        #         type = "rna",
        #         logdir=f'{self.outdir}')

        # else:
        #     ### using high min_umi to only get higher umi solution
        #     cell_bc, count_num = cell_calling(
        #         f"{self.outdir}/02.count/raw_matrix/", 
        #         expected_cell_num = int(self.expectcells),
        #         method = self.calling_method,
        #         min_umi = int(self.minumi) ,
        #         type = "rna",
        #         logdir=f'{self.outdir}')
            

        # generateCellSummary(
        #     f"{self.outdir}/01.data/beads_stat.txt", 
        #     f"{self.outdir}/02.count/barcodeTranslate.txt",
        #     f"{self.outdir}/02.count/raw_matrix",
        #     cell_bc,
        #     f"{self.outdir}/02.count"
        # )    
        
        # print(f'\n{get_formatted_time()}\n'
        #     f'Generating the filtered expression matrix.')
        # str_mkdir('%s/02.count/filter_matrix'%self.outdir)
        # PISA_countFilter_cmd = [
        #     f"{__root_dir__}/software/PISA",
        #     "count",
        #     "-one-hit",
        #     f"-@ {self.threads}",
        #     "-cb DB",
        #     "-anno-tag GN",
        #     "-umi UB",
        #     f"-list {self.outdir}/02.count/beads_barcodes.txt",
        #     f"-outdir {self.outdir}/02.count/filter_matrix",
        #     f"{self.outdir}/02.count/anno_decon_sorted.bam"
        # ]
        # PISA_countFilter_cmd_str = " ".join(PISA_countFilter_cmd)
        # logging_call(PISA_countFilter_cmd_str,'count',self.outdir)  

        # merge_graph(
        #     f"{self.outdir}/02.count/beads_barcodes.txt", 
        #     f"{self.outdir}/02.count"
        #     )     
        
        # # print(f'\n{get_formatted_time()}\t'
        # #     f'Calculate saturation.')
        # sub_sample_cDNA_rna(
        #     '%s/02.count/anno_decon_sorted.bam'%self.outdir,
        #     '%s/02.count/beads_barcodes.txt'%self.outdir,
        #     '%s/02.count'%self.outdir,
        #     threads = self.threads,
        #     quality=20,
        #     logdir=f'{self.outdir}'
        #     )

        # rm_temp(
        #     f'{self.outdir}/02.count/cell.id',
        #     f"{self.outdir}/02.count/cell_count_detail.xls",
        #     f"{self.outdir}/02.count/similarity.dropletfiltered.csv",
        #     f"{self.outdir}/02.count/beads.barcodes.umi100.txt",
        #     f"{self.outdir}/02.count/combined.list.tsv",
        # )

def count(args):
    Count(args).run()


def helpInfo_data(parser):
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
        '-c', '--chemistry',
        metavar='STR',
        choices=["10X","leader_v1","other"],
        help='Chemistry version, can be "10X", "leader_v1", "other". If set to other, needs to be used with --mapparams. [default: leader_v1].',
        default='leader_v1'
        )
    parser.add_argument(
        '-w', '--whitelist',
        metavar='STR',
        help='Path to the cell barcode whitelist file.',
        required=True
        )
    parser.add_argument(
        '-m', '--mapparams',
        metavar='STR',
        help='Additional STAR mapping parameters. must be provide while setting chemistry to other.',
        default=''
        )
    parser.add_argument(
        '-t', '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '-G', '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '-g', '--gtf',
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
    parser.add_argument(
        '--calling_method',
        metavar='STR',
        help='Cell calling method, Choose from barcoderanks and emptydrops, [default: emptydrops].', 
        default='emptydrops'
        )
    parser.add_argument(
        '--expectcells',
        metavar='INT',
        help='Expected number of recovered beads, used as input to cell calling algorithm, [default: 3000].', 
        default=3000
        )
    parser.add_argument(
        '--forcecells',
        help='Force pipeline to use this number of beads, bypassing cell calling algorithm.',
        metavar='INT',
        )
    parser.add_argument(
        '--minumi',
        metavar='INT',
        help='The min umi for use emptydrops, [default: 1000].', 
        default=1000
        )
    return parser