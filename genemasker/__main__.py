#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA, FastICA
from sklearn.experimental import enable_iterative_imputer
from sklearn.preprocessing import StandardScaler
from sklearn.impute import IterativeImputer
from genemasker.definitions import annot_cols
from genemasker.tracking import resource_tracker
import genemasker.fxns as fxns
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import math
import random
import shutil

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

		@resource_tracker(logger)
		def gene_masker(logger, logger_handler):

			logger.info("making temporary directory for storing checkpoints")
			tmpdir = f"{args.out}_tmp"
			os.makedirs(tmpdir, exist_ok=False)

			logger.info("Reading MAF file")
			maf_df = pd.read_csv(args.stat, usecols=[args.stat_id_col, args.stat_maf_col], dtype={args.stat_id_col: "str", args.stat_maf_col: "float32"}, sep="\t")
			maf_df.rename(columns={args.stat_id_col: "#Uploaded_variation", args.stat_maf_col: "MAF"}, inplace=True)

			logger.info("Reading and processing annotation file")
			chunk_paths_orig, chunk_missense_paths_orig, chunk_count_orig, raw_variant_count, stored_variant_count, stored_missense_variant_count, rankscore_cols = fxns.process_annotation_file(annot = args.annot, cols = annot_cols, maf = maf_df, out_dir = tmpdir, chunk_size = args.chunk_size, n_partitions = 1)
			print(chunk_missense_paths_orig)
			ddf = dd.read_parquet(chunk_missense_paths_orig)
			logger.info(f"""Average partition size of raw annot file: {avg_chunk_size(chunk_paths_orig)}""")
			logger.info(f"Found {raw_variant_count} variants in raw annot file and stored {stored_variant_count}")
			total_rows = ddf.shape[0].compute()
			logger.info(f"Found {total_rows} missense variants in annot file")

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

			scaler_fit = None
			if args.impute_standardized:
				logger.info("Fitting standard scaler on full rankscores for imputation")
				scaler = StandardScaler()
				scaler_fit = scaler.fit(ddf[rankscore_cols_keep].compute())

			logger.info("Selecting variants for imputer training data")
			valid_rows = ddf.shape[0].compute()
			logger.info(f"Found {valid_rows} variants with valid rankscores")
			variant_ids_nonmiss = list(ddf[ddf['__non_miss_rankscore_keep__'] == len(rankscore_cols_keep)]['#Uploaded_variation'].compute())
			logger.info(f"Found {len(variant_ids_nonmiss)} variants with no missing rankscores")
			variant_ids_somemiss = list(ddf[ddf['__non_miss_rankscore_keep__'] != len(rankscore_cols_keep)]['#Uploaded_variation'].compute())
			logger.info(f"Found {len(variant_ids_somemiss)} variants with some missing rankscores")

			iter_imputer = IterativeImputer(min_value=0, max_value=1, random_state=1398, max_iter=30, verbose=2, tol=args.impute_tol)

			if args.impute_training_frac:
				logger.info("Selecting random variants for imputer training data")
				random.seed(123)
				random_ids_nonmiss = random.sample(variant_ids_nonmiss, math.ceil(valid_rows * args.impute_training_frac * 0.6))
				random_ids_somemiss = random.sample(variant_ids_somemiss, math.ceil(valid_rows * args.impute_training_frac * 0.4))
				logger.info(f"Setting total rows in training dataset to {math.ceil(valid_rows * args.impute_training_frac)} ({len(random_ids_nonmiss) + len(random_ids_somemiss)} found)")
				logger.info(f"Setting number of variants with no missing rankscores to {math.ceil(valid_rows * args.impute_training_frac * 0.6)} ({len(random_ids_nonmiss)} found)")
				logger.info(f"Setting number of variants with some missing rankscores to {math.ceil(valid_rows * args.impute_training_frac * 0.4)} ({len(random_ids_somemiss)} found)")

				logger.info(f"Extracting imputer training data for columns {rankscore_cols_keep}")
				training_df = fxns.filter_annotation_file(chunk_missense_paths_orig, ["#Uploaded_variation"] + rankscore_cols_keep, random_ids_nonmiss + random_ids_somemiss)
				logger.info(f"Imputer training data contains {training_df.shape[0]} variants")
				if args.impute_standardized:
					logger.info(f"Centering and scaling training data for imputation")
					training_df[rankscore_cols_keep] = scaler_fit.transform(training_df[rankscore_cols_keep])

				logger.info("Fitting iterative imputer for training data")
				iter_imputer_fit = iter_imputer.fit(training_df[rankscore_cols_keep])

				logger.info("Deleting iterative imputer training data")
				del training_df

			else:
				logger.info(f"Extracting imputer training data for columns {rankscore_cols_keep}")
				full_df = fxns.filter_annotation_file(chunk_missense_paths_orig, ["#Uploaded_variation"] + rankscore_cols_keep, variant_ids_nonmiss + variant_ids_somemiss)
				logger.info(f"Imputer full data contains {full_df.shape[0]} variants")
				if args.impute_standardized:
					logger.info(f"Centering and scaling training data for imputation")
					full_df[rankscore_cols_keep] = scaler_fit.transform(full_df[rankscore_cols_keep])

				logger.info("Fitting iterative imputer for full dataset")
				iter_imputer_fit = iter_imputer.fit(full_df[rankscore_cols_keep])

				logger.info("Deleting iterative imputer full data")
				del full_df

			logger.info("Imputing annotation file")
			chunk_paths_imp, rankscore_cols_keep_imputed = fxns.impute_annotation_file(chunk_paths = chunk_missense_paths_orig, iter_imp = iter_imputer_fit, rankscore_cols = rankscore_cols_keep, scaler_fit = scaler_fit, save_all = args.save_all)
			ddf = dd.read_parquet(chunk_paths_imp)
			logger.info(f"""Average partition size of imputed annot file: {avg_chunk_size(chunk_paths_imp)}""")

			imputed_scaler_fit = None
			if args.pca_standardized:
				logger.info("Fitting standard scaler on imputed rankscores for pca")
				imputed_scaler = StandardScaler()
				imputed_scaler_fit = imputed_scaler.fit(ddf[rankscore_cols_keep_imputed].compute())

			pca_fit = None
			full_df_pca = None
			if args.pca_fit_method == 'incremental':
				logger.info(f"Fitting model for incremental PCA on columns {rankscore_cols_keep_imputed}")
				pca_fit = fxns.fit_incremental_pca(chunk_paths_imp, rankscore_cols_keep_imputed, scaler_fit = imputed_scaler_fit)

			else:
				pca = PCA(n_components = len(rankscore_cols_keep_imputed))
				if args.pca_training_frac:
					logger.info("Selecting random variants for use in PCA training data")
					random.seed(456)
					variant_ids = list(ddf['#Uploaded_variation'].compute())
					logger.info(f"Found {len(variant_ids)} variants available for training data")
					random_ids = random.sample(variant_ids, math.ceil(valid_rows * args.pca_training_frac))
					logger.info(f"Setting number of variants for PCA training data to {math.ceil(valid_rows * args.pca_training_frac)}")

					logger.info(f"Extracting PCA training data from imputed rankscores for columns {rankscore_cols_keep_imputed}")
					training_df = fxns.filter_annotation_file(chunk_paths_imp, ["#Uploaded_variation"] + rankscore_cols_keep_imputed, random_ids)
					logger.info(f"ICA training data contains {training_df.shape[0]} variants")
					if args.pca_standardized:
						logger.info(f"Centering and scaling training data for pca")
						training_df[rankscore_cols_keep_imputed] = imputed_scaler_fit.transform(training_df[rankscore_cols_keep_imputed])

					logger.info(f"Fitting model for PCA on columns {rankscore_cols_keep_imputed}")
					pca_fit = pca.fit(training_df[rankscore_cols_keep_imputed])

					logger.info("Deleting PCA training data")
					del training_df

				else:
					logger.info(f"Extracting PCA full data from imputed rankscores for columns {rankscore_cols_keep_imputed}")
					full_df = ddf[["#Uploaded_variation"] + rankscore_cols_keep_imputed].compute()
					logger.info(f"PCA full data contains {full_df.shape[0]} variants")
					if args.pca_standardized:
						logger.info(f"Centering and scaling training data for pca")
						full_df[rankscore_cols_keep_imputed] = imputed_scaler_fit.transform(full_df[rankscore_cols_keep_imputed])

					logger.info(f"Fitting model for PCA on columns {rankscore_cols_keep_imputed}")
					pca_fit = pca.fit(full_df[rankscore_cols_keep_imputed])
					full_df_pca = pd.DataFrame(data = pca_fit.transform(full_df[rankscore_cols_keep_imputed]), columns = [f"pc{i+1}" for i in range(len(rankscore_cols_keep_imputed))])
					full_df_pca['#Uploaded_variation'] = full_df['#Uploaded_variation'].to_list()

					logger.info("Deleting PCA full data")
					del full_df

			pca_explained_variance = pd.DataFrame({'VarianceExplained': pca_fit.explained_variance_ratio_.cumsum()})
			pca_explained_variance.to_csv(f"{args.out}.pca_explained_variance.tsv", sep='\t', index=False)
			logger.info(f"PCA explained variance written to {args.out}.pca_explained_variance.tsv")

			logger.info("Calculating pca scores")
			chunk_paths_pca = fxns.calculate_pca_scores(chunk_paths = chunk_paths_imp, pca_fit = pca_fit, full_df_pca = full_df_pca, pca_n = len(rankscore_cols_keep_imputed), rankscore_cols = rankscore_cols_keep_imputed, pca_scaler_fit = imputed_scaler_fit, save_all = args.save_all)
			ddf = dd.read_parquet(chunk_paths_pca)
			logger.info(f"""Average partition size of imputed missense annot file with pca scores: {avg_chunk_size(chunk_paths_pca)}""")

			fast_ica = FastICA(n_components=None, random_state=0, max_iter=5000, whiten='unit-variance', whiten_solver='eigh')
			ica_cols_keep = [col for col in rankscore_cols_keep if col not in ['Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore']]
			ica_cols_keep_imputed = [x + "_imputed" for x in ica_cols_keep]

			fast_ica_fit = None
			full_df_ica = None
			if args.ica_training_frac:
				logger.info("Selecting random variants for use in ICA training data")
				random.seed(789)
				variant_ids = list(ddf['#Uploaded_variation'].compute())
				logger.info(f"Found {len(variant_ids)} variants available for training data")
				random_ids = random.sample(variant_ids, math.ceil(valid_rows * args.ica_training_frac))
				logger.info(f"Setting number of variants for ICA training data to {math.ceil(valid_rows * args.ica_training_frac)}")

				logger.info(f"Extracting ICA training data from imputed rankscores for columns {ica_cols_keep_imputed}")
				training_df = fxns.filter_annotation_file(chunk_paths_imp, ["#Uploaded_variation"] + ica_cols_keep_imputed, random_ids)
				logger.info(f"ICA training data contains {training_df.shape[0]} variants")

				logger.info(f"Fitting model for ICA on columns {ica_cols_keep_imputed}")
				fast_ica_fit = fast_ica.fit(training_df[ica_cols_keep_imputed])

				logger.info("Deleting ICA training data")
				del training_df

			else:
				logger.info(f"Extracting ICA full data from imputed rankscores for columns {ica_cols_keep_imputed}")
				full_df = ddf[["#Uploaded_variation"] + ica_cols_keep_imputed].compute()
				logger.info(f"ICA full data contains {full_df.shape[0]} variants")

				logger.info(f"Fitting model for ICA on columns {ica_cols_keep_imputed} and transforming full data")
				full_df_ica = pd.DataFrame(data = fast_ica.fit_transform(full_df[ica_cols_keep_imputed]), columns = [f"ic{i+1}" for i in range(len(ica_cols_keep_imputed))])
				full_df_ica['#Uploaded_variation'] = full_df['#Uploaded_variation'].to_list()

				logger.info("Deleting ICA full data")
				del full_df

			logger.info("Calculating ica scores")
			chunk_paths_ica = fxns.calculate_ica_scores(chunk_paths = chunk_paths_pca, ica_fit = fast_ica_fit, ica_cols = ica_cols_keep_imputed, full_df_ica = full_df_ica, save_all = args.save_all)
			ddf = dd.read_parquet(chunk_paths_ica)
			logger.info(f"""Average partition size of imputed missense annot file with ica scores: {avg_chunk_size(chunk_paths_ica)}""")

			logger.info("Calculating Rankscore ~ MAF correlations")
			correlations = fxns.calculate_correlations(ddf, rankscore_cols_keep_imputed, 'MAF')
			correlations.to_csv(f"{args.out}.rankscore_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"Rankscore ~ MAF correlations written to {args.out}.rankscore_maf_corr.tsv")

			logger.info("Calculating PC ~ MAF correlations")
			print(list(ddf.columns))
			pc_correlations = fxns.calculate_correlations(ddf, [f"pc{i+1}" for i in range(len(rankscore_cols_keep_imputed))], 'MAF')
			pc_correlations.to_csv(f"{args.out}.pc_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"PC ~ MAF correlations written to {args.out}.pc_maf_corr.tsv")

			logger.info("Calculating IC ~ MAF correlations")
			ic_correlations = fxns.calculate_correlations(ddf, [f"ic{i+1}" for i in range(len(ica_cols_keep_imputed))], 'MAF')
			ic_correlations.to_csv(f"{args.out}.ic_maf_corr.tsv", sep='\t', index=False)
			logger.info(f"IC ~ MAF correlations written to {args.out}.ic_maf_corr.tsv")

			logger.info("Calculating combined scores")
			chunk_paths_predicted = fxns.calculate_damage_predictions(chunk_paths = chunk_paths_ica, rankscore_cols = rankscore_cols_keep, ica_cols = ica_cols_keep, pca_n = len(rankscore_cols_keep), pc_corr = pc_correlations, pca_exp_var = pca_explained_variance, ic_corr = ic_correlations, save_all = args.save_all)
			logger.info(f"""Average partition size of imputed missense annot file with pca scores, ica scores, and damage prediction scores: {avg_chunk_size(chunk_paths_predicted)}""")

			logger.info("Calculating combined score correlations")
			ddf = dd.read_parquet(chunk_paths_predicted)
			score_corr = fxns.get_damaging_pred_all(ddf)
			score_corr.to_csv(f"{args.out}.combined_score_corr.tsv", sep='\t', index=False)
			logger.info(f"Combined score correlations written to {args.out}.combined_score_corr.tsv")

			logger.info("Merge combined scores into annot file")
			chunk_paths_orig_scored = fxns.merge_annot_scores(chunk_paths_orig)

			logger.info("Calculating gene mask filters")
			chunk_paths_filters = fxns.calculate_mask_filters(chunk_paths_orig_scored, ignore_mask_maf = args.ignore_mask_maf, save_all = args.save_all)
			ddf = dd.read_parquet(chunk_paths_filters)
			logger.info(f"""Average partition size of imputed missense annot file with pca scores, ica scores, damage prediction scores, and filters: {avg_chunk_size(chunk_paths_filters)}""")

			logger.info("Calculating proportion of damage for each prediction algorithm in the annot file")
			pred_columns = [x for x in list(ddf.columns) if x.endswith("_pred")]
			perc_damage = fxns.calculate_percent_damaging(ddf, pred_cols = pred_columns)
			perc_damage.to_csv(f"{args.out}.damaging_prop.tsv", sep='\t', index=False)

			logger.info("Generating Regenie group files")
			fxns.generate_regenie_groupfiles(ddf, out = args.out)

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
