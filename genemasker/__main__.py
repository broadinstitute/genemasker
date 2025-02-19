#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from sklearn.decomposition import FastICA, PCA
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from genemasker.definitions import annot_cols
from genemasker.tracking import resource_tracker
import genemasker.fxns as fxns
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import math
import random
import shutil

def count_notna(row):
	return len([x for x in row if not math.isnan(x)])

def main(args=None):

	args = fxns.args
	logger = fxns.logger
	logger_handler = fxns.logger_handler

	with ProgressBar():

		@resource_tracker(logger)
		def gene_masker(logger, logger_handler):
	
			logger.info("making temporary directory for storing checkpoints")
			tmpdir = f"{args.out}_tmp"
			os.makedirs(tmpdir, exist_ok=False)
		
			logger.info("Reading MAF file")
			maf_df = pd.read_csv(args.stat, usecols=[args.stat_id_col, args.stat_maf_col], sep="\t")
			maf_df.rename(columns={args.stat_id_col: "#Uploaded_variation", args.stat_maf_col: "MAF"}, inplace=True)
		
			logger.info("Reading and processing annotation file")
			chunk_paths_orig, chunk_count_orig = fxns.process_annotation_file(annot = args.annot, cols = annot_cols, maf = maf_df, out_dir = tmpdir, chunk_size = args.chunk_size, n_partitions = args.n_partitions)
			ddf = dd.read_parquet(chunk_paths_orig)
		
			rankscore_cols = [col for col in ddf.columns if col.endswith('rankscore')]
		
			logger.info("Dropping rows with all missing rankscores")
			ddf['__non_miss_rankscore__'] = ddf[rankscore_cols].apply(count_notna, axis=1, meta=('x','int'))
			ddf = ddf.loc[ddf['__non_miss_rankscore__'] > 0]
			
			logger.info("Calculating missing proportions")
			na_count = ddf[[col for col in ddf.columns if col.endswith('rankscore')]].isna().mean().reset_index().compute()
			na_count.columns = ['Algorithm', 'MissingProportion']
			na_count.to_csv(f"{args.out}.rankscore_miss.tsv", sep='\t', index=False)
			logger.info(f"Missing proportions written to {args.out}.rankscore_miss.tsv")
			
			logger.info("Dropping columns with >20% missing data")
			columns_to_drop = list(na_count.loc[na_count['MissingProportion'] > 0.2, 'Algorithm'])
			ddf = ddf.drop(columns=columns_to_drop)
			
			rankscore_cols_keep = [x for x in rankscore_cols if x not in columns_to_drop]
		
			logger.info("Calculating rankscore means on full data")
			rankscore_means = ddf[rankscore_cols_keep].mean().compute()
			logger.info("Calculating rankscore stddevs on full data")
			rankscore_stddevs = ddf[rankscore_cols_keep].std().compute()
		
			logger.info("Selecting random variants for use in IterativeImputer training data")
			total_rows = ddf.shape[0].compute()
			logger.info(f"Found {total_rows} variants with valid rankscores")
			variant_ids_nonmiss = list(ddf[ddf['__non_miss_rankscore__'] == len(rankscore_cols_keep)]['#Uploaded_variation'].compute())
			logger.info(f"Found {len(variant_ids_nonmiss)} variants with no missing rankscores")
			variant_ids_somemiss = list(ddf[ddf['__non_miss_rankscore__'] != len(rankscore_cols_keep)]['#Uploaded_variation'].compute())
			logger.info(f"Found {len(variant_ids_somemiss)} variants with some missing rankscores")
			random.seed(123)
			random_ids_nonmiss = random.sample(variant_ids_nonmiss, math.ceil(total_rows * args.training_frac * 0.6))
			random_ids_somemiss = random.sample(variant_ids_nonmiss, math.ceil(total_rows * args.training_frac * 0.4))
			logger.info(f"Setting total rows in training dataset to {math.ceil(total_rows * args.training_frac)}")
			logger.info(f"Setting number of variants with no missing rankscores to {math.ceil(total_rows * args.training_frac * 0.6)}")
			logger.info(f"Setting number of variants with some missing rankscores to {math.ceil(total_rows * args.training_frac * 0.4)}")
		
			logger.info(f"Extracting training data for columns {rankscore_cols_keep} and centering on full mean values")
			training_df = fxns.filter_annotation_file(chunk_paths_orig, ["#Uploaded_variation"] + rankscore_cols_keep, random_ids_nonmiss + random_ids_somemiss)
			training_df[rankscore_cols_keep] = (training_df[rankscore_cols_keep] - rankscore_means) / rankscore_stddevs
		
			logger.info("Fitting IterativeImputer for training data")
			iter_imputer = IterativeImputer(min_value=0, max_value=1, random_state=1398, max_iter=30, verbose=2)
			iter_imputer_fit = iter_imputer.fit(training_df[rankscore_cols_keep])
		
			logger.info("Deleting IterativeImputer training data")
			del training_df
		
			logger.info("Imputing annotation file")
			chunk_paths_imp = fxns.impute_annotation_file(chunk_paths = chunk_paths_orig, iter_imp = iter_imputer_fit, rankscore_cols = rankscore_cols_keep, means = rankscore_means, stddevs = rankscore_stddevs)
			ddf = dd.read_parquet(chunk_paths_imp)
		
			logger.info("Selecting random variants for use in PCA and ICA training data")
			random.seed(456)
			variant_ids = list(ddf['#Uploaded_variation'].compute())
			logger.info(f"Found {len(variant_ids)} variants available for training data")
			random_ids = random.sample(variant_ids, math.ceil(total_rows * args.training_frac))
			logger.info(f"Setting number of variants for training data to {math.ceil(total_rows * args.training_frac)}")
		
			logger.info(f"Extracting training data from imputed rankscores for columns {rankscore_cols_keep} and centering on full mean values")
			training_df = fxns.filter_annotation_file(chunk_paths_imp, ["#Uploaded_variation"] + rankscore_cols_keep, random_ids)
			training_df[rankscore_cols_keep] = (training_df[rankscore_cols_keep] - rankscore_means) / rankscore_stddevs
			
			logger.info(f"Fitting model for PCA on columns {rankscore_cols_keep} using centered training data")
			n = min(len(training_df), len(rankscore_cols_keep))
			pca = PCA(n_components = n)
			pca_fit = pca.fit(training_df[rankscore_cols_keep])
			pca_explained_variance = pd.DataFrame({'VarianceExplained': pca.explained_variance_ratio_.cumsum()})
		
			ica_cols_keep = [col for col in rankscore_cols_keep if col not in ['Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore', 'MAF']]
			logger.info(f"Fitting model for ICA on columns {ica_cols_keep} using centered training data")
			fast_ica = FastICA(n_components=None, random_state=0, whiten='unit-variance', max_iter=5000, whiten_solver='eigh')
			fast_ica_fit = fast_ica.fit(training_df[ica_cols_keep])
		
			logger.info("Deleting pca and ica training data")
			del training_df
		
			logger.info("Calculating pca and ica scores")
			chunk_paths_scored = fxns.calculate_pca_ica_scores(chunk_paths = chunk_paths_imp, pca_fit = pca_fit, ica_fit = fast_ica_fit, pca_n = n, rankscore_cols = rankscore_cols_keep, ica_cols = ica_cols_keep)
			ddf = dd.read_parquet(chunk_paths_scored)
		
			logger.info("Calculating Rankscore ~ MAF correlations")
			correlations = fxns.calculate_correlations(ddf, rankscore_cols_keep, 'MAF')
			correlations.to_csv(f"{args.out}.rankscore_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"Rankscore ~ MAF correlations written to {args.out}.rankscore_maf_corr.tsv")
		
			logger.info("Calculating PC ~ MAF correlations")
			pc_correlations = fxns.calculate_correlations(ddf, [f"pc{i+1}" for i in range(n)], 'MAF')
			pc_correlations.to_csv(f"{args.out}.pc_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"PC ~ MAF correlations written to {args.out}.pc_maf_corr.tsv")
		
			logger.info("Calculating IC ~ MAF correlations")
			ic_correlations = fxns.calculate_correlations(ddf, [f"ic{i+1}" for i in range(len(ica_cols_keep))], 'MAF')
			ic_correlations.to_csv(f"{args.out}.ic_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"IC ~ MAF correlations written to {args.out}.ic_maf_corr.tsv")
		
			logger.info("Calculating proportion of damage for each prediction algorithm in the annot file")
			pred_columns = [x for x in list(ddf.columns) if x.endswith("_pred")]
			perc_damage = fxns.calculate_percent_damaging(ddf, pred_cols = pred_columns)
			perc_damage.to_csv(f"{args.out}.damaging_prop.tsv", sep='\t', index=False)
		
			logger.info("Calculating combined scores")
			chunk_paths_predicted = fxns.calculate_damage_predictions(chunk_paths = chunk_paths_scored, rankscore_cols = rankscore_cols_keep, ica_cols = ica_cols_keep, pca_n = n, pc_corr = pc_correlations, pca_exp_var = pca_explained_variance, ic_corr = ic_correlations)
			ddf = dd.read_parquet(chunk_paths_predicted)
		
			logger.info("Calculating combined score correlations")
			score_corr = fxns.get_damaging_pred_all(ddf)
			score_corr.to_csv(f"{args.out}.combined_score_corr.tsv", sep='\t', index=False)
			logger.info(f"Combined score correlations written to {args.out}.combined_score_corr.tsv")
		
			logger.info("Calculating gene mask filters")
			chunk_paths_filters = fxns.calculate_mask_filters(chunk_paths_predicted)
			ddf = dd.read_parquet(chunk_paths_filters)
	
			logger.info("Generating Regenie group files")
			fxns.generate_regenie_groupfiles(ddf, out = args.out)
		
			## Write all to file
			#ddf.compute().to_csv(f"{args.out}.results.tsv.gz", sep='\t', index=False, compression='gzip')
		
			logger.info("Execution complete. Cleaning up temporary directory")
			shutil.rmtree(tmpdir)
	
		gene_masker(logger, logger_handler)
	
		# Close logging
		logger_handler.close()

if __name__ == "__main__":
	main()
	os._exit(0)
