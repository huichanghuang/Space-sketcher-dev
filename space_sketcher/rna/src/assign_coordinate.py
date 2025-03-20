import os,gzip, argparse
from collections import defaultdict
import dnaio
from itertools import product
import pandas as pd
from space_sketcher.tools.utils import add_log
import editdistance
"""
This function assign coordinate to each spatial barcode by puckfile
"""

