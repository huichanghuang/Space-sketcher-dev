import os
import sys
import math
from subprocess import check_call
from space_sketcher.tools.utils import judgeFilexits, change_path
from space_sketcher.__init__ import __root_dir__

def count_chromosomes(genome_file):
    chromosomes = set()
    with open(genome_file, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                chromosome = line.strip().lstrip(">")
                chromosomes.add(chromosome)
    return len(chromosomes)

def star_index(fasta,gtf,genomeDir,star_program,limitram,threads):
    if not os.path.exists(genomeDir):
        os.system('mkdir -p %s'%genomeDir)
    genome_size = os.path.getsize(fasta)
    SAindexNbases = min(int(math.log2(genome_size)/2 - 1),14)
    n_entries = count_chromosomes(fasta)
    chr_bins = min(18, int(math.log2(max(int(genome_size / n_entries),100))))
    
    star_cmd = [
        star_program,
        '--runMode genomeGenerate',
        f'--runThreadN {threads}',
        f'--genomeDir {genomeDir}',
        f'--genomeFastaFiles {fasta}',
        f'--sjdbGTFfile {gtf}',
        '--sjdbOverhang 99',
        f'--limitGenomeGenerateRAM {limitram}',
        f'--genomeSAindexNbases {SAindexNbases}',
        f'--genomeChrBinNbits {chr_bins}'
    ]
    star_cmd_str = ' '.join(star_cmd)

    print('STAR verison: 2.7.11b')
    print('runMode: genomeGenerate')
    print('runThreadN: %s'%threads)
    print('limitGenomeGenerateRAM: %s'%limitram)
    print('genomeSAindexNbases: %s'%SAindexNbases)
    print('genomeChrBinNbits: %s'%chr_bins)
    print('genomeDir: %s'%genomeDir)
    print('fasta: %s'%fasta)
    print('gtf: %s'%gtf)

    sys.stdout.flush()
    check_call(star_cmd_str,shell=True)


class Ref:
    def __init__(self, args):
        self.ingtf = os.path.abspath(args.ingtf)
        self.fasta = os.path.abspath(args.fasta)
        self.species = args.species
        self.genomeDir = os.path.abspath(args.genomeDir)
        self.limitram = args.limitram
        self.threads = args.threads

    def run(self):
        change_path()
        judgeFilexits(self.ingtf,self.fasta)
        self.genomeDir = os.path.abspath(self.genomeDir)
        starbin = os.path.join(__root_dir__, "software/STAR")
        if not self.noindex:
            star_index(self.fasta,
                       self.ingtf,
                       self.genomeDir,
                       starbin, 
                       self.limitram, 
                       self.threads)
        print("\033[0;32;40mAnalysis Complete\033[0m")

def mkref(args):
    Ref(args).run()

def helpInfo_mkref(parser):
    parser.add_argument(
        '-f', '--fasta',
        metavar='<FASTA>',
        help='Path to the genome file in FASTA format.'
        )
    parser.add_argument(
        '-i','--ingtf', 
        metavar='<GTF>' ,
        help='Path to the genome annotation file in GTF format.'
        )
    parser.add_argument(
        '-s', '--species',
        metavar='<SPECIES>',
        default='undefined',
        help='Species name, [default: undefined].'
        )
    parser.add_argument(
        '-g', '--genomeDir',
        metavar='<DATABASE>',
        default=os.getcwd(),
        help='Path to the directory where database files will be stored., [default: current dir].'
        )
    parser.add_argument(
        '-l', '--limitram',
        metavar='<MEMORY>',
        help='Maximum available RAM (bytes) for genome index generation.[default: 125000000000]',
        default=125000000000
        )
    parser.add_argument(
        '-t', '--threads',
        metavar='<CORENUM>', 
        default=4,
        help='Number of threads used for analysis, [default: 4].'
        )
    return parser

