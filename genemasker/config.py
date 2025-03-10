import argparse
from genemasker.__version__ import version
import genemasker.logging as logging

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version', version="%(prog)s v" + version)
parser.add_argument('--chunk-size', type=int, default=None, help='number of annot rows per chunk')
parser.add_argument('--n-partitions', type=int, default=1, help='number of partitions per chunk')
parser.add_argument('--pca-fit-method', choices=['incremental','standard'], default='standard', help='PCA model fit method')
parser.add_argument('--pca-training-frac', type=float, help='fraction of data to use for training pca models')
parser.add_argument('--impute-training-frac', type=float, help='fraction of data to use for training imputer models')
parser.add_argument('--ica-training-frac', type=float, help='fraction of data to use for training ica models')
parser.add_argument('--impute-center-scale', action='store_true', help='center and scale data prior to training impute model')
parser.add_argument('--pca-center-scale', action='store_true', help='center and scale data prior to training pca model')
parser.add_argument('--ica-center-scale', action='store_true', help='center and scale data prior to training ica model')
parser.add_argument('--annot', required=True, help='VEP annotation files')
parser.add_argument('--stat', required=True, help='MAF file')
parser.add_argument('--stat-id-col', required=True, help='Column name in MAF matching annotation file')
parser.add_argument('--stat-maf-col', required=True, help='MAF column name')
parser.add_argument('--out', required=True, help='Output file prefix')
args = parser.parse_args()

logger, logger_handler = logging.setup_logger(f"{args.out}.genemasker.log")

logger.info(f"genemasker v{version}")
logger.info("user-supplied arguments:")
for key, value in vars(args).items():
    logger.info(f"  {key}: {value}")
