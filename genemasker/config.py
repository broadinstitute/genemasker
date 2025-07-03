import argparse
from genemasker.__version__ import version
import genemasker.logging as logging

parser = argparse.ArgumentParser()
parser.add_argument('--version', action='version', version="%(prog)s v" + version)
parser.add_argument('--chunk-size', type=int, default=None, help='number of annot rows per chunk')
parser.add_argument('--pca-fit-method', choices=['incremental','standard'], default='standard', help='PCA model fit method')
parser.add_argument('--pca-training-frac', type=float, help='fraction of data to use for training pca models')
parser.add_argument('--impute-training-frac', type=float, help='fraction of data to use for training imputer models')
parser.add_argument('--ica-training-frac', type=float, help='fraction of data to use for training ica models')
parser.add_argument('--impute-standardized', action='store_true', help='standardize prior to training impute model')
parser.add_argument('--pca-standardized', action='store_true', help='standardize prior to training pca model')
parser.add_argument('--impute-tol', type=float, default=0.001, help='imputation tolerance')
parser.add_argument('--ignore-mask-maf', action='store_true', help='ignore maf filters in masks (allow for downstream maf filtering)')
parser.add_argument('--impute-pipeline', help='an impute pipeline from previous run')
parser.add_argument('--pca-pipeline', help='a pca pipeline from previous run')
parser.add_argument('--ica-pipeline', help='an ica pipeline from previous run')
parser.add_argument('--rankscore-maf-corr', help='rankscore~maf correlations from previous run')
parser.add_argument('--pc-maf-corr', help='pc~maf correlations from previous run')
parser.add_argument('--ic-maf-corr', help='ic~maf correlations from previous run')
parser.add_argument('--annot', required=True, help='VEP annotation files')
parser.add_argument('--stat', help='MAF file')
parser.add_argument('--stat-id-col', help='Column name in MAF matching annotation file')
parser.add_argument('--stat-maf-col', help='MAF column name')
parser.add_argument('--save-all', action='store_true', help='save all temporary output files')
parser.add_argument('--out', required=True, help='Output file prefix')
args = parser.parse_args()

reconstruct_option = {}
for action in parser._actions:
	if action.option_strings:
		reconstruct_option[action.dest] = action.option_strings[0]

args_proj = ['impute_pipeline', 'pca_pipeline', 'ica_pipeline', 'rankscore_maf_corr', 'pc_maf_corr', 'ic_maf_corr']
args_proj_txt = "[--impute-pipeline, --pca-pipeline, --ica-pipeline, --rankscore-maf-corr, --pc-maf-corr, --ic-maf-corr]"
args_noproj = ['stat', 'stat_id_col', 'stat_maf_col']
args_noproj_txt = "[--stat, --stat-id-col, --stat-maf-col]"

args_proj_provided = [name for name in args_proj if getattr(args, name) is not None]
args_noproj_provided = [name for name in args_noproj if getattr(args, name) is not None]

if args_proj_provided and args_noproj_provided:
	parser.error(f"Provide either all of {args_proj_txt} or all of {args_noproj_txt}, not both.")
elif args_proj_provided:
	if len(args_proj_provided) < len(args_proj):
		missing = [f'{reconstruct_option[name]}' for name in args_proj if getattr(args, name) is None]
		parser.error(f"Missing arguments from {args_proj_txt}: {', '.join(missing)}.")
elif args_noproj_provided:
	if len(args_noproj_provided) < len(args_noproj):
		missing = [f'{reconstruct_option[name]}' for name in args_noproj if getattr(args, name) is None]
		parser.error(f"Missing arguments from {args_noproj_txt}: {', '.join(missing)}.")
else:
	parser.error(f"You must provide all of {args_proj_txt} or all of {args_noproj_txt}")

logger, logger_handler = logging.setup_logger(f"{args.out}.genemasker.log")

logger.info(f"genemasker v{version}")
logger.info("user-supplied arguments:")
for key, value in vars(args).items():
    logger.info(f"  {key}: {value}")
