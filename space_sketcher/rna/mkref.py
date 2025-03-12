import os
import sys
import math
import re
import json
import collections
from subprocess import check_call
# from dnbc4tools.tools.utils import judgeFilexits, change_path
# from dnbc4tools.__init__ import __root_dir__
# from dnbc4tools.tools.mkgtf import *
from space_sketcher.tools.utils import judgeFilexits, change_path
from space_sketcher.__init__ import __root_dir__
from space_sketcher.tools.mkgtf import *

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

def mtgenes_list(chrM,genomeDir,gtf):
    if chrM == 'auto':
        chrmtlist = ['chrM','MT','chrMT','mt','Mt']
        if os.path.exists(os.path.join(genomeDir,'chrName.txt')):
            with open(os.path.join(genomeDir,'chrName.txt'),'r') as file:
                chrlist =  file.read().splitlines()
                union = list(set(chrmtlist) & set(chrlist))
                if union:
                    chrMT = union[0]
                else:
                    chrMT = 'None'
        else:
            chrMT = 'None'
    else:
        chrMT = chrM

    if chrMT != 'None':
        MTgene = []
        with open(gtf,'r') as gtffile:
            for line in gtffile:
                line = line.strip()
                if line.startswith('#'):
                    pass
                else:
                    lst = line.split('\t')
                    if lst[0] == chrMT:
                        if lst[2] == "gene":
                            aDict = collections.OrderedDict()
                            pattern = re.compile(r'(\S+?)\s*"(.*?)"')
                            for m in re.finditer(pattern, lst[-1]):
                                key = m.group(1)
                                value = m.group(2)
                                aDict[key] = value
                            if "gene_name" in aDict:
                                MTgene.append(aDict['gene_name'])
                            else:
                                MTgene.append(aDict['gene_id'])
        if MTgene:
            result = open(os.path.join(genomeDir,'mtgene.list'),'w')
            for line in MTgene:
                result.write(line+'\n')
            result.close()
    return chrMT

def write_config(genomeDir,species,chrM,fasta,gtf):
    ref_dict = collections.OrderedDict()
    ref_dict['species'] = str(species)
    ref_dict['genome'] = os.path.abspath(fasta)
    ref_dict['gtf'] = os.path.abspath(gtf)
    ref_dict['genomeDir'] = os.path.abspath(genomeDir)
    chrMT = mtgenes_list(chrM,genomeDir,gtf)
    ref_dict['chrmt'] = str(chrMT)
    if os.path.exists(os.path.join(genomeDir,'mtgene.list')):
        ref_dict['mtgenes'] = os.path.abspath("%s/mtgene.list"%genomeDir)
    else:
        ref_dict['mtgenes'] = 'None'
    with open("%s/ref.json"%genomeDir, "w", encoding='utf-8') as jsonfile:
        json.dump(ref_dict, jsonfile, indent=4,  ensure_ascii=False) 
        jsonfile.write("\n")

class Ref:
    def __init__(self, args):
        self.ingtf = os.path.abspath(args.ingtf)
        self.fasta = os.path.abspath(args.fasta)
        self.species = args.species
        self.chrM = args.chrM
        self.genomeDir = os.path.abspath(args.genomeDir)
        self.limitram = args.limitram
        self.threads = args.threads
        self.noindex = args.noindex

    def run(self):
        change_path()
        judgeFilexits(self.ingtf,self.fasta)
        self.genomeDir = os.path.abspath(self.genomeDir)
        if not self.noindex:
            star_index(self.fasta,self.ingtf,self.genomeDir,"%s/software/STAR"%__root_dir__,self.limitram,self.threads)
        write_config(self.genomeDir,self.species,self.chrM,self.fasta,self.ingtf)
        print("\033[0;32;40mAnalysis Complete\033[0m")

def mkref(args):
    Ref(args).run()

def helpInfo_mkref(parser):
    parser.add_argument(
        '--fasta',
        metavar='<FASTA>',
        help='Path to the genome file in FASTA format.'
        )
    parser.add_argument(
        '--ingtf', 
        metavar='<GTF>' ,
        help='Path to the genome annotation file in GTF format.'
        )
    parser.add_argument(
        '--species',
        metavar='<SPECIES>',
        default='undefined',
        help='Species name, [default: undefined].'
        )
    parser.add_argument(
        '--chrM',
        metavar='<MT>',
        help='Mitochondrial chromosome name, [default: auto].',
        default='auto'
        )
    parser.add_argument(
        '--genomeDir',
        metavar='<DATABASE>',
        default=os.getcwd(),
        help='Path to the directory where database files will be stored., [default: current dir].'
        )
    parser.add_argument(
        '--limitram',
        metavar='<MEMORY>',
        help='Maximum available RAM (bytes) for genome index generation.',
        default=125000000000
        )
    parser.add_argument(
        '--threads',
        metavar='<CORENUM>', 
        default=4,
        help='Number of threads used for analysis, [default: 4].'
        )
    parser.add_argument(
        '--noindex',
        action='store_true',
        help='Only generate ref.json without constructing the genome index.'
        )
    return parser



# ####cell-sketcher
# import subprocess
# from loguru import logger
# from pathlib import Path
# from cell_sketcher.utils import get_package_dir

# source_dir = get_package_dir()

# def build_rna_star_index(args):
#     genome = Path(args.genome)
#     if genome.exists():
#         raise ValueError(f"目录 {genome} 已存在，删除文件夹后重新运行")
#     genome.mkdir(parents=True,exist_ok=False)
#     star_bin = source_dir / "bin" / "STAR"
#     cmd = (
#         f'{star_bin} --runMode genomeGenerate '
#         f'--runThreadN {args.thread} '
#         f'--genomeDir {genome} '
#         f'--genomeFastaFiles {args.fasta} '
#         f'--sjdbGTFfile {args.gtf} '
#     )
#     logger.info(f"构建索引命令：{cmd}")
#     subprocess.check_call(cmd, shell=True)

#     if args.mt_genes is not None:
#         subprocess.check_call(['cp -f', args.mt_genes, str(genome / "mtGene.txt")])
#     return None
