import os,shutil

class Report:
    def __init__(self,args):
        self.name = args.name
        self.species = args.species
        self.threads = args.threads
        self.no_introns = args.no_introns
        self.end5 = args.end5
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))

    def run(self):
        print("test report")
        # ### import lib
        # from dnbc4tools.tools.utils import str_mkdir,logging_call,\
        #     judgeFilexits, get_formatted_time, rm_temp
        # from dnbc4tools.tools.io import read_anndata
        # from dnbc4tools.__init__ import __root_dir__
        # from dnbc4tools.rna.src.generate_report import write_param_to_template

        # ### run
        # rm_temp(
        #     f'{self.outdir}/04.report/div',
        #     f'{self.outdir}/04.report/table'
        # )

        # str_mkdir(f'{self.outdir}/04.report/div')
        # str_mkdir(f'{self.outdir}/04.report/table')
        # str_mkdir(f'{self.outdir}/log')

        # judgeFilexits(
        #     f'{self.outdir}/01.data/cDNA.barcode.counts.txt',
        #     f'{self.outdir}/01.data/oligo.sequencing.report.csv',
        #     f'{self.outdir}/01.data/alignment_report.csv',
        #     f'{self.outdir}/01.data/anno_report.csv',
        #     f'{self.outdir}/02.count/singlecell.csv',
        #     f'{self.outdir}/02.count/saturation_cDNA.xls',
        #     f'{self.outdir}/02.count/raw_matrix',
        #     f'{self.outdir}/02.count/filter_matrix',
        #     f'{self.outdir}/03.analysis/cluster.csv',
        #     f'{self.outdir}/03.analysis/raw_qc.xls',
        #     f'{self.outdir}/03.analysis/marker.csv'   
        # )

        # ### plotly
        # print(f'\n{get_formatted_time()}\n'
        #     f'Statistical analysis and report generation for results.')

        # ### html 
        # intron = 'False' if self.no_introns else 'True' 
        # end5 = " 5\'" if self.end5 else ""
        # report_file, metrics_summary_df = write_param_to_template(
        #     ['%s/config/template/template_scRNA.html'%__root_dir__,
        #      '%s/config/template/template_scRNA_noanno.html'%__root_dir__],
        #     self.name,
        #     self.outdir,
        #     intron,
        #     self.species,
        #     end5
        #     )
        # htmlFile = open('%s/04.report/%s_scRNA_report.html'%(self.outdir,self.name),'w')
        # htmlFile.write(report_file)
        # htmlFile.close()
        # metrics_summary_df.to_csv('%s/04.report/metrics_summary.xls'%self.outdir,sep='\t',index=None)

        # ### copy file
        # str_mkdir('%s/output'%self.outdir)
        # if os.path.exists('%s/04.report/metrics_summary.xls'%self.outdir):
        #     shutil.copy("%s/04.report/metrics_summary.xls"%self.outdir,'%s/output'%self.outdir)
        # if os.path.exists('%s/04.report/%s_scRNA_report.html'%(self.outdir,self.name)):
        #     shutil.copy("%s/04.report/%s_scRNA_report.html"%(self.outdir,self.name),'%s/output'%self.outdir)
        # if os.path.exists('%s/02.count/singlecell.csv'%self.outdir):
        #     shutil.copy("%s/02.count/singlecell.csv"%self.outdir,'%s/output'%self.outdir)
        # if os.path.exists('%s/03.analysis/filter_feature.h5ad'%self.outdir):
        #     shutil.copy('%s/03.analysis/filter_feature.h5ad'%self.outdir,'%s/output'%self.outdir)
        # if os.path.exists("%s/02.count/filter_matrix"%self.outdir):
        #     if os.path.exists('%s/output/filter_matrix'%self.outdir):
        #         shutil.rmtree('%s/output/filter_matrix'%self.outdir)
        #     shutil.copytree("%s/02.count/filter_matrix"%self.outdir,'%s/output/filter_matrix'%self.outdir,dirs_exist_ok=True)
        # if os.path.exists("%s/02.count/raw_matrix"%self.outdir):
        #     if os.path.exists('%s/output/raw_matrix'%self.outdir):
        #         shutil.rmtree('%s/output/raw_matrix'%self.outdir)
        #     shutil.copytree("%s/02.count/raw_matrix"%self.outdir,'%s/output/raw_matrix'%self.outdir,dirs_exist_ok=True)

        # adata = read_anndata(f"{self.outdir}/output/filter_matrix/")
        # adata.write(f"{self.outdir}/output/filter_feature.h5ad")
        
        # ### move bam
        # if os.path.exists('%s/02.count/anno_decon_sorted.bam'%self.outdir):
        #     if(os.path.exists("%s/output/anno_decon_sorted.bam"%self.outdir)):
        #         os.remove("%s/output/anno_decon_sorted.bam"%self.outdir)
        #     shutil.move("%s/02.count/anno_decon_sorted.bam"%self.outdir,'%s/output'%self.outdir)
        # if os.path.exists('%s/02.count/anno_decon_sorted.bam.bai'%self.outdir):
        #     if(os.path.exists("%s/output/anno_decon_sorted.bam.bai"%self.outdir)):
        #         os.remove("%s/output/anno_decon_sorted.bam.bai"%self.outdir)
        #     shutil.move("%s/02.count/anno_decon_sorted.bam.bai"%self.outdir,'%s/output'%self.outdir)
        # if os.path.exists('%s/02.count/anno_decon_sorted.bam.csi'%self.outdir):
        #     if(os.path.exists("%s/output/anno_decon_sorted.bam.csi"%self.outdir)):
        #         os.remove("%s/output/anno_decon_sorted.bam.csi"%self.outdir)
        #     shutil.move("%s/02.count/anno_decon_sorted.bam.csi"%self.outdir,'%s/output'%self.outdir)

        # ### other matrix
        # splice_matrix_cmd = [
        #     f"{__root_dir__}/software/PISA",
        #     f"count -one-hit -@ {self.threads}",
        #     f"-cb DB -ttype E,S -anno-tag GN -umi UB",
        #     f"-list {self.outdir}/02.count/beads_barcodes.txt",
        #     f"-outdir {self.outdir}/output/attachment/splice_matrix",
        #     f"{self.outdir}/output/anno_decon_sorted.bam"
        # ]

        # splice_matrix_cmd_str = " ".join(splice_matrix_cmd)

        # RNAvelocity_matrix_cmd = [
        #     f"{__root_dir__}/software/PISA",
        #     f"count -one-hit -@ {self.threads}",
        #     f"-cb DB -velo -anno-tag GN -umi UB",
        #     f"-list {self.outdir}/02.count/beads_barcodes.txt",
        #     f"-outdir {self.outdir}/output/attachment/RNAvelocity_matrix",
        #     f"{self.outdir}/output/anno_decon_sorted.bam"
        # ]

        # RNAvelocity_matrix_cmd_str = " ".join(RNAvelocity_matrix_cmd)


        # if not self.no_introns:
        #     str_mkdir('%s/output/attachment/splice_matrix'%self.outdir)
        #     str_mkdir('%s/output/attachment/RNAvelocity_matrix'%self.outdir)
        #     logging_call(splice_matrix_cmd_str,'report',self.outdir)
        #     logging_call(RNAvelocity_matrix_cmd_str,'report',self.outdir)
        # else:
        #     if os.path.exists('%s/output/attachment/splice_matrix'%self.outdir):
        #         shutil.rmtree('%s/output/attachment/splice_matrix'%self.outdir)
        #     if os.path.exists('%s/output/attachment/RNAvelocity_matrix'%self.outdir):
        #         shutil.rmtree('%s/output/attachment/RNAvelocity_matrix'%self.outdir)

def report(args):
    Report(args).run()

def helpInfo_report(parser):
    parser.add_argument(
        '--name',
        required=True,
        help='Sample name.'
        )
    parser.add_argument(
        '--species',
        type=str, 
        metavar='STR',
        help='species.',
        required=True
        )
    parser.add_argument(
        '--outdir',
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='INT',
        default=4,
        help='Analysis threads, [default: 4].'
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Not include intronic reads in count.'
        )
    parser.add_argument(
        '--end5', 
        action='store_true',
        help='5-end single-cell transcriptome.'
        )
    return parser
