import os,collections
import argparse
import time
from typing import List

class Runpipe:
    def __init__(self, args):
        self.name = args.name
        self.cDNAr1 = args.cDNAfastq1
        self.cDNAr2 = args.cDNAfastq2
        self.oligor1 = args.oligofastq1
        self.oligor2 = args.oligofastq2
        self.genomeDir = os.path.abspath(args.genomeDir)
        
        self.outdir = os.path.abspath(args.outdir)
        self.threads = args.threads
        self.chemistry = args.chemistry
        self.darkreaction = args.darkreaction
        self.customize = args.customize
        self.calling_method = args.calling_method
        self.expectcells = args.expectcells
        self.forcecells = args.forcecells
        
        self.process: List[str] = args.process.split(",")
        self.no_introns = args.no_introns
        self.end5 = args.end5
        self.minumi = args.minumi
        self.outunmappedreads = args.outunmappedreads
        
    def runpipe(self):
        print("test run")
        ### import lib
        # from dnbc4tools.tools.utils import str_mkdir,judgeFilexits,read_json,start_print_cmd,bin_path
        # from dnbc4tools.__init__ import __root_dir__

        # ### run
        # judgeFilexits('%s/ref.json'%self.genomeDir)
        # indexConfig = read_json('%s/ref.json'%self.genomeDir)
        # genomeDir = indexConfig['genomeDir']
        # gtf = indexConfig['gtf']
        # chrmt = indexConfig['chrmt']
        # species = indexConfig['species']
        # mtgenes = indexConfig['mtgenes']
        # judgeFilexits(self.cDNAr1,self.cDNAr2,self.oligor1,self.oligor2,genomeDir,gtf)
        # if mtgenes != "None":
        #     judgeFilexits(mtgenes)

        # data_cmd = [
        #     f"{bin_path()}/dnbc4tools rna data",
        #     f"--cDNAfastq1 {self.cDNAr1}",
        #     f"--cDNAfastq2 {self.cDNAr2}",
        #     f"--oligofastq1 {self.oligor1}",
        #     f"--oligofastq2 {self.oligor2}",
        #     f"--threads {self.threads}",
        #     f"--name {self.name}",
        #     f"--chemistry {self.chemistry}",
        #     f"--darkreaction {self.darkreaction}",
        #     f"--outdir {self.outdir}",
        #     f"--genomeDir {genomeDir}",
        #     f"--gtf {gtf}",
        #     f"--chrMT {chrmt}"
        # ]
        # if self.customize:
        #     data_cmd += ['--customize %s'%self.customize]
        # if self.no_introns:
        #     data_cmd += ['--no_introns']
        # if self.end5:
        #     data_cmd += ['--end5']
        # if self.outunmappedreads:
        #     data_cmd += ['--outunmappedreads']
        # data_cmd = ' '.join(data_cmd)

        # count_cmd = [
        #     f"{bin_path()}/dnbc4tools rna count",
        #     f"--name {self.name}",
        #     f"--calling_method {self.calling_method}",
        #     f"--expectcells {self.expectcells}", 
        #     f"--threads {self.threads}",
        #     f"--outdir {self.outdir}"
        #     ]
        # if self.forcecells:
        #     count_cmd += ['--forcecells %s'%self.forcecells]
        # if self.minumi:
        #     count_cmd += ['--minumi %s'%self.minumi]
        # count_cmd  = ' '.join(count_cmd)
                

        # analysis_cmd = [
        #     f'{bin_path()}/dnbc4tools rna analysis',
        #     f"--name {self.name}",
        #     f"--outdir {self.outdir}",
        #     f"--species {species}",
        #     f"--mtgenes {mtgenes}"
        # ]
        # analysis_cmd  = ' '.join(analysis_cmd)
        
        
        # report_cmd = [
        #     f'{bin_path()}/dnbc4tools rna report',
        #     f'--name {self.name}',
        #     f"--species {species}",
        #     f"--threads {self.threads}",
        #     f"--outdir {self.outdir}"
        # ]
        # if self.no_introns:
        #     report_cmd += ['--no_introns']
        # if self.end5:
        #     report_cmd += ['--end5']
        # report_cmd = ' '.join(report_cmd)

       
        # # Assemble the list of commands to execute based on requested stages
        # # Append corresponding commands to cmdlist based on keywords present in `self.process`
        # allowed_processes = {'data', 'count', 'analysis', 'report'}
        # if not all(process in allowed_processes for process in self.process):
        #     raise ValueError("Invalid process name(s) detected. Allowed names are: data, count, analysis, report.")
        
        # cmdlist = collections.OrderedDict()
        # if 'data' in self.process:
        #     cmdlist['data'] = data_cmd
        # if 'count' in self.process:
        #     cmdlist['count'] = count_cmd
        # if 'analysis' in self.process:
        #     cmdlist['analysis'] = analysis_cmd
        # if 'report' in self.process:
        #     cmdlist['report'] = report_cmd

        # start_time = time.time()
        # str_mkdir('%s/log'%os.path.join(self.outdir,self.name))
        # for pipe,pipecmd in cmdlist.items():
        #     #logging_call(pipecmd,pipe,os.path.join(self.outdir,self.name))
        #     start_print_cmd(pipecmd,os.path.join(self.outdir,self.name))
        # end_time = time.time()

        # analysis_time = end_time - start_time
        # analysis_time_minutes, analysis_time_seconds = divmod(analysis_time, 60)
        # analysis_time_hours, analysis_time_minutes = divmod(analysis_time_minutes, 60)

        # print(f'\nAnalysis Finished')
        # print(f'Elapsed Time: {int(analysis_time_hours)} hours {int(analysis_time_minutes)} minutes {int(analysis_time_seconds)} seconds')



def run(args):
    Runpipe(args).runpipe()

def helpInfo_run(parser):
    parser.add_argument(
        '--name', 
        metavar='<SAMPLE_ID>',
        help='User-defined sample ID.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--cDNAfastq1', 
        metavar='<FQ1FILES>',
        help='Paths to the raw R1 FASTQ files of the cDNA library.', 
        required=True
        )
    parser.add_argument(
        '--cDNAfastq2', 
        metavar='<FQ2FILES>',
        help='Paths to the raw R2 FASTQ files of the cDNA library.', 
        required=True
        )
    parser.add_argument(
        '--oligofastq1', 
        metavar='<FQ1FILES>',
        help='Paths to the raw R1 FASTQ files of the oligo library.',
        required=True
        )
    parser.add_argument(
        '--oligofastq2', 
        metavar='<FQ2FILES>',
        help='Paths to the raw R2 FASTQ files of the oligo library.',
        required=True
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='<DATABASE>',
        help='Path to the directory containing genome files.',
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='<OUTDIR>',
        help='Output directory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--threads',
        type=int, 
        metavar='<CORENUM>',
        default=4,
        help='Number of threads used for analysis, [default: 4].'
        )
    parser.add_argument(
        '--calling_method',
        metavar='<CELLCALLING>',
        choices=["barcoderanks","emptydrops"],
        help='Cell calling method, choose from barcoderanks and emptydrops, [default: emptydrops].', 
        default='emptydrops'
        )
    parser.add_argument(
        '--expectcells',
        metavar='<CELLNUM>',
        help='Expected number of recovered cells, [default: 3000].', 
        default=3000
        )
    parser.add_argument(
        '--forcecells',
        metavar='<CELLNUM>',
        help='Force pipeline to use a specific number of cells.', 
        )
    parser.add_argument(
        '--minumi',
        metavar='<MINUMI>',
        help=argparse.SUPPRESS,
        )
    parser.add_argument(
        '--chemistry',
        metavar='<CHEMISTRY>',
        choices=["scRNAv1HT","scRNAv2HT","scRNAv3HT","auto"],
        help='Chemistry version. Automatic detection recommended , [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--darkreaction',
        metavar='<DARKCYCLE>',
        help='Sequencing dark cycles. Automatic detection recommended, [default: auto].', 
        default='auto'
        )
    parser.add_argument(
        '--customize',
        metavar='<STRUCTURE>',
        help='Customize whitelist and readstructure files in JSON format for cDNA and oligo.'
        )
    parser.add_argument(
        '--process', 
        metavar='<ANALYSIS_STEPS>',
        help='Custom analysis steps to skip unnecessary processes, [default: data,count,analysis,report].',
        type=str,
        default='data,count,analysis,report'
        )
    parser.add_argument(
        '--no_introns', 
        action='store_true',
        help='Intron reads are not included in the expression matrix.'
        )
    parser.add_argument(
        '--end5', 
        action='store_true',
        help='Perform 5\'-end single-cell transcriptome analysis.'
        )
    parser.add_argument(
        '--outunmappedreads',
        action='store_true',
        help=argparse.SUPPRESS,
        )
    return parser