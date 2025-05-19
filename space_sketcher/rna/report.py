import os,shutil
from pathlib import Path

class Report:
    def __init__(self,args):
        self.outdir = os.path.abspath(os.path.join(args.outdir,args.name))
        self.name = args.name
        self.kit = args.chemistry
        self.reference = args.reference
        self.dev = args.dev

    def run(self):
        ### import lib
        from space_sketcher.tools.utils import str_mkdir
        from space_sketcher.rna.src.generate_report import generate_report

        ###generate report
        reportdir = os.path.join(self.outdir, "04.report")
        str_mkdir(reportdir) 
        ### run
        generate_report(
            self.outdir,
            self.name,
            self.kit,
            self.reference,
            self.dev
        )
        
        ###copy files
        os.system(f"cp -r {self.outdir}/02.oligo/SCST {reportdir}")
        (Path(self.outdir) / ".report.done").touch()


def report(args):
    Report(args).run()

def helpInfo_report(parser):
    parser.add_argument(
        '-n', '--name',
        type=str,
        metavar='STR',
        help='Sample name.'
        )
    parser.add_argument(
        '-o', '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '-rf', '--reference',
        type=str, 
        metavar='STR',
        help='Reference.'
        )
    parser.add_argument(
        '-c', '--chemistry',
        metavar='STR',
        choices=["10X","leader_v1","other"],
        help='Chemistry version, can be "10X", "leader_v1", "other". If set to other, needs to be used with --mapparams. [default: leader_v1].',
        default='leader_v1'
        )
    parser.add_argument(
        '-d', '--dev',
        type=lambda x: x.lower() == 'true',
        choices=[True, False],
        default=True,
        help='Development mode, selected from [False, True],(default: True)'
    )
    return parser
