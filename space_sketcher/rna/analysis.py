import os
from pathlib import Path

class Analysis:
    def __init__(self,args):
        self.name = args.name
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        self.minfeatures = args.minfeatures
        self.mincells = args.mincells
        self.ndims = args.ndims
        self.nvariables = args.nvariables
        self.resolution = args.resolution

    
    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir, judgeFilexits, get_formatted_time
        from space_sketcher.__init__ import __root_dir__
        from space_sketcher.rna.src.matrix_QC import perform_matrix_qc
        print("test analysis")
        
        anadir = os.path.join(self.outdir, "03.analysis")
        str_mkdir(anadir)

        ###Extract spatial barcodes
        # judge file exits
        matrixdir = os.path.join(self.outdir, "02.oligo/SCST")
        judgeFilexits(matrixdir)
        # run perform_matrix_qc
        print(f'\n{get_formatted_time()}\n'
            f'Performing matrix QC.')
        perform_matrix_qc(matrixdir, anadir, self.minfeatures, self.mincells,
                        self.nvariables, self.ndims, self.resolution)
        
        (Path(self.outdir) / ".analysis.done").touch()

def analysis(args):
    Analysis(args).run()

def helpInfo_analysis(parser):
    parser.add_argument(
        '-n','--name',
        metavar='NAME',
        required=True,
        help='Sample name.'
        )
    parser.add_argument(
        '-o','--outdir',
        metavar='DIR',
        help='output dir, [default: current directory].',
        default=os.getcwd()
        )
    parser.add_argument(
        '-f', '--minfeatures',
        type=int, 
        metavar='INT',
        help='Minimum features per cell. default: 5',
        default=5
        )
    parser.add_argument(
        '-c', '--mincells',
        type=int, 
        metavar='INT',
        help='Minimum cells per gene. default: 3',
        default=3
        )
    parser.add_argument(
        '-d', '--ndims',
        type=int, 
        metavar='INT',
        help='PCA dimensions. default: 30',
        default=30
        )
    parser.add_argument(
        '-v', '--nvariables',
        type=int, 
        metavar='INT',
        help='Number of variable genes. default: 2000',
        default=2000
        )
    parser.add_argument(
        '-r', '--resolution',
        type=float, 
        metavar='FLOAT',
        help='Clustering resolution. default: 0.5',
        default=0.5
        )
    return parser