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
        self.cbwhitelist = args.cbwhitelist
        self.mapparams = args.mapparams
        self.genomeDir = args.genomeDir
        self.calling_method = args.calling_method
        self.velo = args.velo
        # self.expectcells = args.expectcells
        # self.forcecells = args.forcecells
        # self.minumi = args.minumi
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))

    def Prepare_mapping_params(self) -> str:
        """
        Check the chemistry and prepare mapping parameters.
        Returns:
            str: The mapping parameters.
        """
        ###load function
        from space_sketcher.tools.utils import gunzip
        ###to avoid STAR memery error
        if self.cbwhitelist.endswith(".gz"):
            self.cbwhitelist = gunzip(self.cbwhitelist)

        mapping_params = ""
        if self.velo:
            mapping_params += "--soloFeatures GeneFull_Ex50pAS Velocyto "
        else:
            mapping_params += "--soloFeatures GeneFull_Ex50pAS "
                    
        if self.chemistry == "10X":
            mapping_params += "--soloType CB_UMI_Simple "
            mapping_params += "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 "
            mapping_params += "--soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloBarcodeReadLength 0 "
            mapping_params += f"--soloCBwhitelist {self.cbwhitelist} "
        elif self.chemistry == "leader_v1":
            mapping_params += "--soloType CB_UMI_Complex "
            mapping_params += "--soloCBposition 0_0_0_9 0_10_0_19 --soloUMIposition 0_20_0_29 "
            mapping_params += "--soloCBmatchWLtype EditDist_2 "
            mapping_params += f"--soloCBwhitelist {self.cbwhitelist} {self.cbwhitelist} "
        elif self.chemistry == "other":
            if self.mapparams == 'None':
                logger.info("Please check if chemistry and mapparams provided in proper way!")
        else:
            logger.info("Not available chemistry")

        if not self.mapparams == 'None':
            mapping_params += self.mapparams

        return mapping_params
        
    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.saturation import count_saturation
        ### run
        judgeFilexits(
            self.rna1,
            self.rna2,
            self.cbwhitelist,
            self.genomeDir,
            )

        rnadir = os.path.join(self.outdir, "01.count")
        str_mkdir(rnadir)
        
        star_version = subprocess.check_output(f"{__root_dir__}/software/STAR --version", shell=True)
        logger.info(f"STAR 版本号：{star_version.decode('utf8')}")

        mapping_pars = self.Prepare_mapping_params()
        STAR_cmd = (
            f"{__root_dir__}/software/STAR "
            f"--runMode alignReads "
            "--quantMode GeneCounts "
            f"--soloCellFilter {self.calling_method} "
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
            "--soloOutFileNames Solo.out/ genes.tsv barcodes.tsv matrix.mtx "
            f"{mapping_pars} "
            f"--readFilesIn {self.rna2} {self.rna1} "
            f"--genomeDir {self.genomeDir} "
            f"--outFileNamePrefix {rnadir}/ "
            f"--runThreadN {self.threads} "
            "--outBAMsortingThreadN 6 "
        )

        if not os.path.exists(Path(self.outdir) / ".STAR.done"):
            subprocess.check_call(STAR_cmd, shell=True)
            ###change output directory permissions
            chmod_cmd = f"chmod -R a+r {rnadir}/Solo.out && chmod a+x $(find {rnadir}/Solo.out -type d)"
            subprocess.check_call(chmod_cmd, shell=True)
            (Path(self.outdir) / ".STAR.done").touch()
        else:
            print(".STAR.done exits, skip running STAR.")

        ###calculate saturation
        count_saturation(rnadir, self.threads)

        ###Prepare barcode rank file for RNA knee plot
        judgeFilexits(
            f"{rnadir}/Solo.out/GeneFull_Ex50pAS/raw",
            f"{rnadir}/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv",
            )
        
        import scanpy as sc
        import pandas as pd

        truecells = pd.read_csv(f"{rnadir}/Solo.out/GeneFull_Ex50pAS/filtered/barcodes.tsv", header=None)[0].tolist()
        adata = sc.read_10x_mtx(
            f"{rnadir}/Solo.out/GeneFull_Ex50pAS/raw",
            var_names='gene_symbols',
            cache=True
        )
        umi_counts = pd.DataFrame({
            'barcode': adata.obs_names,
            'UMI': adata.X.sum(axis=1).A1  # 转换为numpy数组
        })
        umi_counts['is_cell_barcode'] = umi_counts['barcode'].isin(truecells).astype(int)
        umi_counts = umi_counts.sort_values('UMI', ascending=False)
        umi_counts['rank'] = range(1, len(umi_counts)+1)
        umi_counts.to_csv(os.path.join(rnadir, 'cell_rna_umi.rank.txt'), sep='\t', index=False)

        (Path(self.outdir) / ".count.done").touch()


def count(args):
    Count(args).run()


def helpInfo_count(parser):
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
        '-w', '--cbwhitelist',
        metavar='FILE',
        help='Path to the cell barcode whitelist file.',
        required=True
        )
    parser.add_argument(
        '-m', '--mapparams',
        metavar='STR',
        help='Additional STAR mapping parameters. must be provide while setting chemistry to other.',
        default='None',
        )
    parser.add_argument(
        '-t', '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Number of threads to use, [default: 4].'
        )
    parser.add_argument(
        '-g', '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path to the directory where genome files are stored.',
        required=True
        )
    parser.add_argument(
        '-ca', '--calling_method',
        metavar='STR',
        choices=["CellRanger2.2","EmptyDrops_CR"],
        help='Cell calling method, Choose from CellRanger2.2 and EmptyDrops_CR, [default: EmptyDrops_CR].', 
        default='EmptyDrops_CR'
        )
    parser.add_argument(
        '-v', '--velo',
        type=lambda x: x.lower() == 'true',
        choices=[True, False],
        default=True,
        help='Run STARsolo Velocyto, selected from [False, True], default=True.'
    )
    # parser.add_argument(
    #     '--expectcells',
    #     metavar='INT',
    #     help='Expected number of recovered beads, used as input to cell calling algorithm, [default: 3000].', 
    #     default=3000
    #     )
    # parser.add_argument(
    #     '--forcecells',
    #     help='Force pipeline to use this number of beads, bypassing cell calling algorithm.',
    #     metavar='INT',
    #     )
    # parser.add_argument(
    #     '--minumi',
    #     metavar='INT',
    #     help='The min umi for use emptydrops, [default: 1000].', 
    #     default=1000
    #     )
    return parser