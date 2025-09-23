#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA, FastICA
from sklearn.experimental import enable_iterative_imputer
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.impute import IterativeImputer
from genemasker.tracking import resource_tracker
import genemasker.fxns as fxns
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import math
import random
import glob
import shutil
import joblib

def avg_chunk_size(chunk_list):
	total_size = sum(os.path.getsize(chunk) for chunk in chunk_list if os.path.isfile(chunk))
	result = total_size / len(chunk_list) if chunk_list else 0
	for u in ['b', 'kb', 'mb', 'gb', 'tb']:
		if result < 1024:
			return f"{result:.2f} {u}"
		result = result / 1024
	return f"{result:.2f} tb"

def main(args=None):

	args = fxns.args
	logger = fxns.logger
	logger_handler = fxns.logger_handler

	with ProgressBar():

		if args.user_definitions:
			logger.info(f"Loading user definitions: {args.user_definitions}")
			fxns.load_user_defined_definitions(args.user_definitions)
			
		if args.user_defined_filters:
			logger.info(f"Loading user defined filters: {args.user_defined_filters}")
			fxns.load_user_defined_filters(args.user_defined_filters)

		run_masks = []
		if args.run_masks_file:
			with open(args.run_masks_file, 'r') as file:
				run_masks = [item.strip() for item in file.readlines()]
		elif args.run_masks:
			run_masks = args.run_masks.split(',')

		@resource_tracker(logger)
		def gene_masker(logger, logger_handler):
			
			if not args.generate_from_filtered:
				
				if not args.generate_from_scored:

					logger.info("making temporary directory for storing checkpoints")
					tmpdir = f"{args.out}_tmp"
					os.makedirs(tmpdir, exist_ok=False)
		
					stat_df = None
					if args.stat:
						stat_cols = [args.stat_id_col] + [x for x in [args.stat_maf_col, args.stat_mac_col] if x is not None]
						stat_cols_dtype = ['str'] + ['float32' for x in [args.stat_maf_col, args.stat_mac_col] if x is not None]
						stat_cols_rename = ["#Uploaded_variation"] + ["MAF" for x in [args.stat_maf_col] if x is not None] + ["MAC" for x in [args.stat_mac_col] if x is not None]
						logger.info("Reading stat file")
						stat_df = pd.read_csv(args.stat, usecols=stat_cols, dtype=dict(zip(stat_cols, stat_cols_dtype)), sep="\t")
						stat_df.rename(columns=dict(zip(stat_cols, stat_cols_rename)), inplace=True)
						logger.info(f"loaded stat file with columns: {stat_cols} -> {list(stat_df.columns)}")
		
					logger.info("Reading and processing annotation file")
					chunk_paths_orig, chunk_missense_paths_orig, chunk_count_orig, raw_variant_count, stored_variant_count, stored_missense_variant_count, rankscore_cols, var_id, uid = fxns.process_annotation_file(annot = args.annot, stat = stat_df, out_dir = tmpdir, chunk_size = args.chunk_size, n_partitions = 1, conserved_domains_only = args.conserved_domains_only, include_transcripts = args.include_transcripts)
					ddf = dd.read_parquet(chunk_missense_paths_orig)
					logger.info(f"""Average partition size of raw annot file: {avg_chunk_size(chunk_paths_orig)}""")
					logger.info(f"Found {raw_variant_count} variants in raw annot file and stored {stored_variant_count}")
					total_rows = ddf.shape[0].compute()
					logger.info(f"Found {total_rows} missense variants in annot file")
		
					if args.rankscore_miss:
						logger.info(f"Importing rankscore missingness from {args.rankscore_miss}")
						na_count = pd.read_csv(args.rankscore_miss, dtype={'Algorithm': 'str', 'MissingProportion': 'float32'}, sep="\t")
					else:
						logger.info("Calculating missing proportions")
						na_count = ddf[[col for col in ddf.columns if col.endswith('rankscore')]].isna().mean().reset_index().compute()
						na_count.columns = ['Algorithm', 'MissingProportion']
						na_count.sort_values('MissingProportion',ascending=False).to_csv(f"{args.out}.rankscore_miss.tsv", sep='\t', index=False)
						logger.info(f"Missing proportions written to {args.out}.rankscore_miss.tsv")
					
					columns_to_drop = list(na_count.loc[na_count['MissingProportion'] > 0.2, 'Algorithm'])
					logger.info(f"Dropping ranskcores with >20% missing data: {columns_to_drop}")
					ddf = ddf.drop(columns=columns_to_drop)
					rankscore_cols_keep = [x for x in rankscore_cols if x not in columns_to_drop]
					logger.info(f"Keeping ranskcores: {rankscore_cols_keep}")
		
					logger.info("Update non-missing rankscores count")
					ddf['__non_miss_rankscore_keep__'] = ddf[rankscore_cols_keep].apply(fxns.count_notna, axis=1, meta=('x','int'))
		
					if args.impute_pipeline:
						logger.info(f"Importing imputation pipeline from {args.impute_pipeline}")
						impute_pipeline = joblib.load(args.impute_pipeline)
					else:
						impute_tasks = []
						if args.impute_standardized:
							logger.info("Adding standard scaler to imputation pipeline")
							impute_tasks.append(('scaler', StandardScaler(with_mean=True, with_std=False)))
			
						logger.info("Selecting variants for imputer training data")
						valid_rows = ddf.shape[0].compute()
						logger.info(f"Found {valid_rows} variants with valid rankscores")
						variant_ids_nonmiss = list(ddf[ddf['__non_miss_rankscore_keep__'] == len(rankscore_cols_keep)][uid].compute())
						logger.info(f"Found {len(variant_ids_nonmiss)} variants with no missing rankscores")
						variant_ids_somemiss = list(ddf[ddf['__non_miss_rankscore_keep__'] != len(rankscore_cols_keep)][uid].compute())
						logger.info(f"Found {len(variant_ids_somemiss)} variants with some missing rankscores")
			
						impute_tasks.append(('iterative_imputer', IterativeImputer(min_value=0, max_value=1, random_state=1398, max_iter=30, verbose=2, tol=args.impute_tol)))
						impute_pipeline = Pipeline(impute_tasks)
			
						if args.impute_training_frac:
							logger.info("Selecting random variants for imputer training data")
							random.seed(123)
							random_ids_nonmiss = random.sample(variant_ids_nonmiss, math.ceil(valid_rows * args.impute_training_frac * 0.6))
							random_ids_somemiss = random.sample(variant_ids_somemiss, math.ceil(valid_rows * args.impute_training_frac * 0.4))
							logger.info(f"Setting total rows in training dataset to {math.ceil(valid_rows * args.impute_training_frac)} ({len(random_ids_nonmiss) + len(random_ids_somemiss)} found)")
							logger.info(f"Setting number of variants with no missing rankscores to {math.ceil(valid_rows * args.impute_training_frac * 0.6)} ({len(random_ids_nonmiss)} found)")
							logger.info(f"Setting number of variants with some missing rankscores to {math.ceil(valid_rows * args.impute_training_frac * 0.4)} ({len(random_ids_somemiss)} found)")
			
							logger.info(f"Extracting imputer training data for columns {rankscore_cols_keep}")
							training_df = fxns.filter_annotation_file(chunk_missense_paths_orig, [uid] + rankscore_cols_keep, random_ids_nonmiss + random_ids_somemiss, uid)
							logger.info(f"Imputer training data contains {training_df.shape[0]} variants")
			
							logger.info("Fitting iterative imputer for training data")
							impute_pipeline.fit(training_df[rankscore_cols_keep])
			
							logger.info("Deleting iterative imputer training data")
							del training_df
			
						else:
							logger.info(f"Extracting imputer training data for columns {rankscore_cols_keep}")
							full_df = fxns.filter_annotation_file(chunk_missense_paths_orig, [uid] + rankscore_cols_keep, variant_ids_nonmiss + variant_ids_somemiss, uid)
							logger.info(f"Imputer full data contains {full_df.shape[0]} variants")
			
							logger.info("Fitting iterative imputer for full dataset")
							impute_pipeline.fit(full_df[rankscore_cols_keep])
			
							logger.info("Deleting iterative imputer full data")
							del full_df
		
					logger.info("Imputing annotation file")
					if args.save_all:
						joblib.dump(impute_pipeline, f"{args.out}.impute_pipeline.pkl")
					chunk_paths_imp, rankscore_cols_keep_imputed = fxns.impute_annotation_file(chunk_paths = chunk_missense_paths_orig, impute_pipeline = impute_pipeline, rankscore_cols = rankscore_cols_keep, save_all = args.save_all)
					ddf = dd.read_parquet(chunk_paths_imp)
					logger.info(f"""Average partition size of imputed annot file: {avg_chunk_size(chunk_paths_imp)}""")
		
					if args.pca_pipeline:
						logger.info(f"Importing pca pipeline from {args.pca_pipeline}")
						pca_pipeline = joblib.load(args.pca_pipeline)
					else:
						pca_tasks = []
						pca_scaler = None
						if args.pca_standardized:
							logger.info("Adding standard scaler on imputed rankscores to pca pipeline")
							pca_scaler = StandardScaler(with_mean=True, with_std=False)
			
						if args.pca_fit_method == 'incremental':
							logger.info(f"Fitting incremental PCA on columns {rankscore_cols_keep_imputed}")
							ipca = IncrementalPCA(n_components = len(rankscore_cols))
							ipca, pca_scaler = fxns.fit_incremental_pca(chunk_paths_imp, rankscore_cols_keep_imputed, ipca = ipca, pca_scaler = pca_scaler)
							pca_pipeline = Pipeline([pca_scaler, ipca])
			
						else:
							pca_tasks.append(('scaler', pca_scaler))
							pca_tasks.append(('pca', PCA(n_components = len(rankscore_cols_keep_imputed))))
							pca_pipeline = Pipeline(pca_tasks)
							if args.pca_training_frac:
								logger.info("Selecting random variants for use in PCA training data")
								random.seed(456)
								variant_ids = list(ddf[uid].compute())
								logger.info(f"Found {len(variant_ids)} variants available for training data")
								random_ids = random.sample(variant_ids, math.ceil(valid_rows * args.pca_training_frac))
								logger.info(f"Setting number of variants for PCA training data to {math.ceil(valid_rows * args.pca_training_frac)}")
			
								logger.info(f"Extracting PCA training data from imputed rankscores for columns {rankscore_cols_keep_imputed}")
								training_df = fxns.filter_annotation_file(chunk_paths_imp, [uid] + rankscore_cols_keep_imputed, random_ids, uid)
								logger.info(f"ICA training data contains {training_df.shape[0]} variants")
			
								logger.info(f"Fitting model for PCA on columns {rankscore_cols_keep_imputed}")
								pca_pipeline.fit(training_df[rankscore_cols_keep_imputed])
			
								logger.info("Deleting PCA training data")
								del training_df
			
							else:
								logger.info(f"Extracting PCA full data from imputed rankscores for columns {rankscore_cols_keep_imputed}")
								full_df = ddf[[uid] + rankscore_cols_keep_imputed].compute()
								logger.info(f"PCA full data contains {full_df.shape[0]} variants")
			
								logger.info(f"Fitting model for PCA on columns {rankscore_cols_keep_imputed}")
								pca_pipeline.fit(full_df[rankscore_cols_keep_imputed])
			
								logger.info("Deleting PCA full data")
								del full_df
		
					if args.save_all:
						joblib.dump(pca_pipeline, f"{args.out}.pca_pipeline.pkl")
					pca_explained_variance = pd.DataFrame({'VarianceExplained': pca_pipeline.named_steps['pca'].explained_variance_ratio_.cumsum()})
					pca_explained_variance.to_csv(f"{args.out}.pca_explained_variance.tsv", sep='\t', index=False)
					logger.info(f"PCA explained variance written to {args.out}.pca_explained_variance.tsv")
		
					logger.info("Calculating pca scores")
					chunk_paths_pca = fxns.calculate_pca_scores(chunk_paths = chunk_paths_imp, pca_pipeline = pca_pipeline, pca_n = len(rankscore_cols_keep_imputed), rankscore_cols = rankscore_cols_keep_imputed, save_all = args.save_all)
					ddf = dd.read_parquet(chunk_paths_pca)
					logger.info(f"""Average partition size of imputed missense annot file with pca scores: {avg_chunk_size(chunk_paths_pca)}""")
		
					ica_cols_keep = [col for col in rankscore_cols_keep if col not in ['Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore']]
					ica_cols_keep_imputed = [x + "_imputed" for x in ica_cols_keep]
					
					if args.ica_pipeline:
						logger.info(f"Importing ica pipeline from {args.ica_pipeline}")
						ica_pipeline = joblib.load(args.ica_pipeline)
					else:
						ica_pipeline = Pipeline([('ica', FastICA(n_components=None, random_state=0, max_iter=5000, whiten='unit-variance', whiten_solver='eigh'))])
						if args.ica_training_frac:
							logger.info("Selecting random variants for use in ICA training data")
							random.seed(789)
							variant_ids = list(ddf[uid].compute())
							logger.info(f"Found {len(variant_ids)} variants available for training data")
							random_ids = random.sample(variant_ids, math.ceil(valid_rows * args.ica_training_frac))
							logger.info(f"Setting number of variants for ICA training data to {math.ceil(valid_rows * args.ica_training_frac)}")
			
							logger.info(f"Extracting ICA training data from imputed rankscores for columns {ica_cols_keep_imputed}")
							training_df = fxns.filter_annotation_file(chunk_paths_imp, [uid] + ica_cols_keep_imputed, random_ids, uid)
							logger.info(f"ICA training data contains {training_df.shape[0]} variants")
			
							logger.info(f"Fitting model for ICA on columns {ica_cols_keep_imputed}")
							ica_pipeline.fit(training_df[ica_cols_keep_imputed])
			
							logger.info("Deleting ICA training data")
							del training_df
			
						else:
							logger.info(f"Extracting ICA full data from imputed rankscores for columns {ica_cols_keep_imputed}")
							full_df = ddf[[uid] + ica_cols_keep_imputed].compute()
							logger.info(f"ICA full data contains {full_df.shape[0]} variants")
			
							logger.info(f"Fitting model for ICA on columns {ica_cols_keep_imputed} and transforming full data")
							ica_pipeline.fit(full_df[ica_cols_keep_imputed])
			
							logger.info("Deleting ICA full data")
							del full_df
		
					if args.save_all:
						joblib.dump(ica_pipeline, f"{args.out}.ica_pipeline.pkl")
					logger.info("Calculating ica scores")
					chunk_paths_ica = fxns.calculate_ica_scores(chunk_paths = chunk_paths_pca, ica_pipeline = ica_pipeline, ica_cols = ica_cols_keep_imputed, save_all = args.save_all)
					ddf = dd.read_parquet(chunk_paths_ica)
					logger.info(f"""Average partition size of imputed missense annot file with ica scores: {avg_chunk_size(chunk_paths_ica)}""")
		
					if args.rankscore_maf_corr:
						logger.info("Importing Rankscore ~ MAF correlations from previous run")
						correlations = pd.read_csv(args.rankscore_maf_corr, sep="\t")
					else:
						logger.info("Calculating Rankscore ~ MAF correlations")
						correlations = fxns.calculate_correlations(ddf, rankscore_cols_keep_imputed, 'MAF')
						correlations.to_csv(f"{args.out}.rankscore_maf_corr.tsv", sep='\t', index=False)
						logger.info(f"Rankscore ~ MAF correlations written to {args.out}.rankscore_maf_corr.tsv")
		
					if args.pc_maf_corr:
						logger.info("Importing PC ~ MAF correlations from previous run")
						pc_correlations = pd.read_csv(args.pc_maf_corr, sep="\t")
					else:
						logger.info("Calculating PC ~ MAF correlations")
						pc_correlations = fxns.calculate_correlations(ddf, [f"pc{i+1}" for i in range(len(rankscore_cols_keep_imputed))], 'MAF')
						pc_correlations.to_csv(f"{args.out}.pc_maf_corr.tsv", sep='\t', index=False)
						logger.info(f"PC ~ MAF correlations written to {args.out}.pc_maf_corr.tsv")
		
					if args.ic_maf_corr:
						logger.info("Importing IC ~ MAF correlations from previous run")
						ic_correlations = pd.read_csv(args.ic_maf_corr, sep="\t")
					else:
						logger.info("Calculating IC ~ MAF correlations")
						ic_correlations = fxns.calculate_correlations(ddf, [f"ic{i+1}" for i in range(len(ica_cols_keep_imputed))], 'MAF')
						ic_correlations.to_csv(f"{args.out}.ic_maf_corr.tsv", sep='\t', index=False)
						logger.info(f"IC ~ MAF correlations written to {args.out}.ic_maf_corr.tsv")
		
					logger.info("Calculating combined scores")
					chunk_paths_predicted = fxns.calculate_damage_predictions(chunk_paths = chunk_paths_ica, rankscore_cols = rankscore_cols_keep, ica_cols = ica_cols_keep, pca_n = len(rankscore_cols_keep), pc_corr = pc_correlations, pca_exp_var = pca_explained_variance, ic_corr = ic_correlations, uid = uid, save_all = args.save_all)
					logger.info(f"""Average partition size of imputed missense annot file with pca scores, ica scores, and damage prediction scores: {avg_chunk_size(chunk_paths_predicted)}""")
		
					logger.info("Calculating combined score correlations")
					ddf = dd.read_parquet(chunk_paths_predicted)
					score_corr = fxns.get_damaging_pred_all(ddf)
					score_corr.to_csv(f"{args.out}.combined_score_corr.tsv", sep='\t', index=False)
					logger.info(f"Combined score correlations written to {args.out}.combined_score_corr.tsv")
		
					logger.info("Merge combined scores into annot file")
					chunk_paths_orig_scored = fxns.merge_annot_scores(chunk_paths_orig, uid)
				
				else:
					chunk_paths_orig_scored = glob.glob(args.generate_from_scored)
					
				logger.info("Calculating gene mask filters")
				chunk_paths_filters = fxns.calculate_mask_filters(chunk_paths_orig_scored, run_masks, save_all = args.save_all)
			
			else:
				chunk_paths_filters = glob.glob(args.generate_from_filtered)
			
			ddf = dd.read_parquet(chunk_paths_filters)
			logger.info(f"""Average partition size of imputed missense annot file with pca scores, ica scores, damage prediction scores, and filters: {avg_chunk_size(chunk_paths_filters)}""")

			if not args.skip_calc_perc_damaging:
				logger.info("Calculating proportion of damage for each prediction algorithm in the annot file")
				pred_columns = [x for x in list(ddf.columns) if x.endswith("_pred")]
				perc_damage = fxns.calculate_percent_damaging(ddf, pred_cols = pred_columns)
				perc_damage.to_csv(f"{args.out}.damaging_prop.tsv", sep='\t', index=False)

			logger.info("Generating Regenie group files")
			fxns.generate_regenie_groupfiles(ddf, run_masks, var_id = var_id, uid = uid, out = args.out)

			# Write all to file
			ddf.compute().to_csv(f"{args.out}.results.tsv.gz", sep='\t', index=False, compression='gzip')

			logger.info("Execution complete. Cleaning up temporary directory")
			#shutil.rmtree(tmpdir)

		gene_masker(logger, logger_handler)

		# Close logging
		logger_handler.close()

if __name__ == "__main__":
	main()
	os._exit(0)
