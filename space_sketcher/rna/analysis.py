import os

class Analysis:
    def __init__(self,args):
        self.name = args.name
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        self.species = args.species
        self.mtgenes = args.mtgenes
    
    def run(self):
        print("test analysis")

        ### import lib
        # from dnbc4tools.tools.utils import str_mkdir,judgeFilexits,\
        #     rm_temp, get_formatted_time,setup_logging
        # from dnbc4tools.tools.io import read_anndata
        # from dnbc4tools.tools.plot_draw import violinrnac_plot
        # from dnbc4tools.__init__ import __root_dir__
        # import matplotlib
        # matplotlib.use('Agg')
        # import warnings
        # import pandas as pd
        # warnings.filterwarnings("ignore")
        # import matplotlib.pyplot as plt
        # import scanpy as sc
        # from contextlib import redirect_stdout
        # from dnbc4tools.rna.src.scanpy_cluster import CellDataAnalyzer,run_AnnoMCA_HCL,\
        #     mtgene_file_read,get_markers,get_cluster
        
        # ### run
        # judgeFilexits(
        #     f'{self.outdir}/02.count/filter_matrix',
        #     )
        # str_mkdir(f'{self.outdir}/03.analysis')
        # str_mkdir(f'{self.outdir}/log')
        # str_mkdir(f'{self.outdir}/log/.temp')

        # os.environ[ 'MPLCONFIGDIR' ] = f'{self.outdir}/log/.temp'
        # os.environ[ 'NUMBA_CACHE_DIR' ] = f'{self.outdir}/log/.temp'
        
        # rm_temp(
        #     f'{self.outdir}/03.analysis/cluster_annotation.png'
        #     )
        
        # print(f'\n{get_formatted_time()}\n'
        #     f'Conducting dimensionality reduction and clustering.')
        

        # adata = read_anndata(
        #     '%s/02.count/filter_matrix/'%self.outdir
        #     )
        # adata.write(f"{self.outdir}/03.analysis/filter_feature.h5ad")

        # scanpy_object = CellDataAnalyzer(adata)
        # if self.mtgenes != 'None':
        #     mtgenelist = mtgene_file_read(self.mtgenes)
        #     scanpy_object._pp_qc(mtgenelist)
        # else:
        #     scanpy_object._pp_qc()
        # fig1 = violinrnac_plot(scanpy_object.adata)
        # fig1.savefig(f'{self.outdir}/03.analysis/raw_QCplot.png', dpi = 150)
        
        # scanpy_object.adata.obs.to_csv(
        #     f"{self.outdir}/03.analysis/raw_qc.xls",index=True,sep='\t'
        #     )
        # scanpy_object._pp_basicfilter()
        
        # if len(scanpy_object.adata.obs.index) > 100:
        #     logger = setup_logging('analysis', self.outdir)
        #     with open(logger.handlers[0].baseFilename, "a") as log_file:
        #         with redirect_stdout(log_file):
        #             scanpy_object._pp_scanpy_doubletdetect()
        
        # fig2 = violinrnac_plot(scanpy_object.adata)
        # fig2.savefig(f'{self.outdir}/03.analysis/filter_QCplot.png', dpi = 150)
        # scanpy_object._pp_hvgs()
        # scanpy_object._pp_reduce()
        # scanpy_object._pp_cluster()
        # scanpy_object._pp_deg()

        # if 'rank_genes_groups' in scanpy_object.adata.uns:
        #     marker_table = get_markers(scanpy_object.adata,cutoff=0.1)
        #     marker_table.to_csv(
        #         f"{self.outdir}/03.analysis/marker.csv",
        #         index=False
        #         )
        # else:
        #     with open(f"{self.outdir}/03.analysis/marker.csv", 'w') as markerfile:
        #         markerfile.write('cluster,gene,score,avg_log2FC,p_val,p_val_adj,pct.1,pct.2')

        # if self.species in ["human","Human","hg19","hg38","Homo_sapiens"]:
        #     try:
        #         # scanpy_object._pp_autoanno("Human")
        #         # sc.pl.umap(scanpy_object.adata, color = ['majority_voting', 'louvain'], legend_loc = 'on data')
        #         ref_df = pd.read_csv(
        #             f'{__root_dir__}/config/rna/HCL.ref.csv.gz', 
        #             sep=',', 
        #             index_col=0
        #             ).T
        #         scanpy_object.adata = run_AnnoMCA_HCL(
        #             ref_df,
        #             scanpy_object.adata,
        #             4
        #             )
        #         sc.pl.umap(
        #             scanpy_object.adata, 
        #             color = ['celltype', 'leiden'], 
        #             legend_loc = 'on data'
        #             )
        #         plt.savefig(f"{self.outdir}/03.analysis/cluster_annotation.png")
        #     except Exception as e:
        #         sc.pl.umap(
        #             scanpy_object.adata, 
        #             color = ['leiden'], 
        #             legend_loc = 'on data'
        #             )
        #         plt.savefig(f"{self.outdir}/03.analysis/cluster.png")
        
        # elif self.species in ["mouse","Mouse","mm10","Mus_musculus"]:
        #     try:
        #         # scanpy_object._pp_autoanno("Mouse")
        #         # sc.pl.umap(scanpy_object.adata, color = ['majority_voting', 'louvain'], legend_loc = 'on data')
        #         ref_df = pd.read_csv(
        #             f'{__root_dir__}/config/rna/MCA.ref.csv.gz', 
        #             sep=',', 
        #             index_col=0
        #             ).T
        #         scanpy_object.adata = run_AnnoMCA_HCL(
        #             ref_df,
        #             scanpy_object.adata,
        #             4
        #             )
        #         sc.pl.umap(
        #             scanpy_object.adata, 
        #             color = ['celltype', 'leiden'], 
        #             legend_loc = 'on data'
        #             )
        #         plt.savefig(f"{self.outdir}/03.analysis/cluster_annotation.png")
        #     except Exception as e:
        #         sc.pl.umap(
        #             scanpy_object.adata, 
        #             color = ['leiden'], 
        #             legend_loc = 'on data'
        #             )
        #         plt.savefig(f"{self.outdir}/03.analysis/cluster.png")
        # else:
        #     sc.pl.umap(
        #         scanpy_object.adata, 
        #         color = ['leiden'], 
        #         legend_loc = 'on data'
        #         )
        #     plt.savefig(f"{self.outdir}/03.analysis/cluster.png")
        
        # cluster_table = get_cluster(scanpy_object.adata)
        # cluster_table.to_csv(
        #     f"{self.outdir}/03.analysis/cluster.csv",
        #     index=True
        #     )
        # scanpy_object.adata.write(f"{self.outdir}/03.analysis/QC_Clutser.h5ad")

def analysis(args):
    Analysis(args).run()

def helpInfo_analysis(parser):
    parser.add_argument(
        '--name',
        metavar='NAME',
        required=True,
        help='Sample name.'
        )
    parser.add_argument(
        '--species',
        type=str, 
        metavar='STR',
        help='Species.',
        required=True
        )
    parser.add_argument(
        '--mtgenes',
        type=str, 
        metavar='PATH',
        help='Path to the mtgenes files.',
        required=True
        )
    parser.add_argument(
        '--outdir',
        metavar='DIR',
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    return parser