import argparse
from genemasker.__version__ import version
import genemasker.logging as logging

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version', version="%(prog)s v" + version)
parser.add_argument('--chunk-size', type=int, default=100000, help='number of annot rows per chunk')
parser.add_argument('--n-partitions', type=int, default=10, help='number of partitions per chunk')
parser.add_argument('--training-frac', type=float, default=0.25, help='fraction of data to use for training models')
parser.add_argument('--annot', required=True, help='VEP annotation files')
parser.add_argument('--stat', required=True, help='MAF file')
parser.add_argument('--stat-id-col', required=True, help='Column name in MAF matching annotation file')
parser.add_argument('--stat-maf-col', required=True, help='MAF column name')
parser.add_argument('--out', required=True, help='Output file prefix')
args = parser.parse_args()

logger, logger_handler = logging.setup_logger(f"{args.out}.genemasker.log")
