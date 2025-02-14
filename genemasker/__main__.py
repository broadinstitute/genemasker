#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import os
from sklearn.decomposition import FastICA, PCA
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from genemasker.definitions import annot_cols
from genemasker.tracking import resource_tracker
import genemasker.fxns as fxns
from genemasker.__version__ import version
import dask.dataframe as dd
import math

def count_notna(row):
	return len([x for x in row if not math.isnan(x)])

@resource_tracker
def main(args=None):

	parser = argparse.ArgumentParser()
	parser.add_argument('--version', action='version', version="%(prog)s v" + version)
	parser.add_argument('--chunk-size', type=int, default=100000, help='number of annot rows per chunk')
	parser.add_argument('--n-partitions', type=int, default=10, help='number of partitions per chunk')
	parser.add_argument('--annot', required=True, help='VEP annotation files')
	parser.add_argument('--stat', required=True, help='MAF file')
	parser.add_argument('--stat-id-col', required=True, help='Column name in MAF matching annotation file')
	parser.add_argument('--stat-maf-col', required=True, help='MAF column name')
	parser.add_argument('--out', required=True, help='Output file prefix')
	args = parser.parse_args()

	# Read and process annotation files
	print("Reading and processing annotation files")
	chunk_paths = fxns.process_annotation_file(annot = args.annot, cols = annot_cols, out_dir = f"{args.out}_tmp", chunk_size = args.chunk_size, n_partitions = args.n_partitions)
	ddf = dd.read_parquet(chunk_paths)

	# Drop rows with all NaN rankscore columns
	print("Dropping rows with all missing rankscores")
	ddf['__non_miss_rankscore__'] = ddf[[col for col in ddf.columns if col.endswith('rankscore')]].apply(count_notna, axis=1, meta=('x','int'))
	ddf = ddf.loc[ddf['__non_miss_rankscore__'] > 0]
	#ddf.dropna(axis=0, how='all', subset=[col for col in ddf.columns if col.endswith('rankscore')], inplace=True)
	#ddf.reset_index(drop=True, inplace=True)
	
	# Calculate missing proportions
	print("Calculating missing proportions")
	na_count = ddf[[col for col in ddf.columns if col.endswith('rankscore')]].isna().mean().reset_index().compute()
	na_count.columns = ['Algorithm', 'MissingProportion']
	na_count.to_csv(f"{args.out}.rankscore_miss.tsv", sep='\t', index=False)
	print(f"Missing proportions written to {args.out}.rankscore_miss.tsv")
	
	# Drop columns with >20% missing data
	print("Dropping columns with >20% missing data")
	columns_to_drop = na_count.loc[na_count['MissingProportion'] > 0.2, 'Algorithm']
	ddf = ddf.drop(columns=columns_to_drop)

	#rankscore_cols = [col for col in ddf.columns if col.endswith('rankscore')]
	#ddf[rankscore_cols].compute().to_csv(f"{args.out}.TEST.tsv", sep='\t', index=False)
	#return

	# Impute missing data
	print("Imputing missing data")
	imputer = IterativeImputer(min_value=0, max_value=1, random_state=1398, max_iter=30, verbose=2)
	rankscore_cols = [col for col in ddf.columns if col.endswith('rankscore')]
	df = ddf[rankscore_cols + ["#Uploaded_variation"]].compute().reset_index(drop=True)
	df[rankscore_cols] = imputer.fit_transform(df[rankscore_cols])
	print("Imputation complete")

	# Merge with MAF file
	print("Merging with MAF file")
	maf = pd.read_csv(args.stat, usecols=[args.stat_id_col, args.stat_maf_col], sep="\t")
	maf.rename(columns={args.stat_id_col: "#Uploaded_variation", args.stat_maf_col: "MAF"}, inplace=True)
	df = df.merge(maf, on="#Uploaded_variation", how="left")
	print("Merge complete")

	# Calculate rankscore correlations
	correlations = fxns.calculate_correlations(df, rankscore_cols, 'MAF')
	correlations.to_csv(f"{args.out}.rankscore_maf_corr.tsv", sep='\t', index=False)
	print(f"Rankscore correlations written to {args.out}.rankscore_maf_corr.tsv")

	# Write imputed rankscores to file
	df.to_csv(f"{args.out}.imputed_rankscores.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"Imputed rankscores written to {args.out}.imputed_rankscores.tsv.gz")

	# PCA
	print("Performing PCA")
	n = min(len(df), len(rankscore_cols))
	pca = PCA(n_components = n)
	principal_components = pca.fit_transform(df[rankscore_cols])
	explained_variance = pd.DataFrame({'VarianceExplained': pca.explained_variance_ratio_.cumsum()})
	#explained_variance.to_csv(f"{args.out}.pca_var_expl.tsv", sep='\t', index=False)
	print(f"PCA variance explained written to {args.out}.pca_var_expl.tsv")
	#pd.DataFrame(list(pca.components_), columns = rankscore_cols).to_csv(f"{args.out}.pca_components.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"PCA components written to {args.out}.pca_components.tsv.gz")

	df = pd.concat([df, pd.DataFrame(data = principal_components, columns = [f"pc{i+1}" for i in range(n)])], axis = 1)
	#df[["#Uploaded_variation"] + [f"pc{i+1}" for i in range(n)]].to_csv(f"{args.out}.pca_scores.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"PCA scores written to {args.out}.pca_scores.tsv.gz")

	# Calculate PC correlations
	pc_correlations = fxns.calculate_correlations(df, [f"pc{i+1}" for i in range(n)], 'MAF')
	pc_correlations.to_csv(f"{args.out}.pc_maf_corr.tsv", sep='\t', index=False)
	print(f"PC correlations written to {args.out}.pc_maf_corr.tsv")

	# ICA
	print("Performing ICA")
	ica_cols = [col for col in rankscore_cols if col not in ['Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore', 'MAF']]
	transformer = FastICA(n_components=None, random_state=0, whiten='unit-variance', max_iter=5000, whiten_solver='eigh')
	X_transformed = transformer.fit_transform(df[ica_cols])
	df = pd.concat([df, pd.DataFrame(data = X_transformed, columns = [f"ic{i+1}" for i in range(len(ica_cols))])], axis = 1)
	#df[["#Uploaded_variation"] + [f"ic{i+1}" for i in range(len(ica_cols))]].to_csv(f"{args.out}.ica_scores.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"ICA scores written to {args.out}.ica_scores.tsv.gz")
	#pd.DataFrame(list(transformer.components_), columns = ica_cols).to_csv(f"{args.out}.ica_components.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"ICA components written to {args.out}.ica_components.tsv.gz")
	#pd.DataFrame(list(transformer.mixing_), columns = [f"ic{i+1}" for i in range(len(ica_cols))]).to_csv(f"{args.out}.pca_components.tsv.gz", sep='\t', index=False, compression='gzip')
	print(f"ICA mixing written to {args.out}.ica_mixing.tsv.gz")

	# Calculate IC correlations
	ic_correlations = fxns.calculate_correlations(df, [f"ic{i+1}" for i in range(len(ica_cols))], 'MAF')
	ic_correlations.to_csv(f"{args.out}.ic_maf_corr.tsv", sep='\t', index=False)
	print(f"IC correlations written to {args.out}.ic_maf_corr.tsv")

	# Calculate percent damaging
	perc_damage = fxns.calculate_percent_damaging(df)
	perc_damage.to_csv(f"{args.out}.damaging_prop.tsv", sep='\t', index=False)

	# Calculate Combo_mean_score_og
	print("Calculating Combo_mean_score_og")
	df['Combo_mean_score_og'] = fxns.get_damaging_pred_og(df[rankscore_cols])
	#df.to_csv(f"{args.out}.damaging_pred.og.tsv.gz", sep='\t', index=False, compression='gzip')

	# Calculate Combo_mean_score_ic
	print("Calculating Combo_mean_score_ic")
	df['Combo_mean_score_ic'] = fxns.get_damaging_pred_ica(df[[f"ic{i+1}" for i in range(len(ica_cols))]], ic_correlations)
	#df.to_csv(f"{args.out}.damaging_pred.ica.tsv.gz", sep='\t', index=False, compression='gzip')

	# Calculate Combo_mean_score_pc
	print("Calculating Combo_mean_score_pc")
	df['Combo_mean_score_pc'] = fxns.get_damaging_pred_pca(df[[f"pc{i+1}" for i in range(n)]], pc_correlations, explained_variance)
	#df.to_csv(f"{args.out}.damaging_pred.pca.tsv.gz", sep='\t', index=False, compression='gzip')

	# Calculate Combo_mean_score correlations
	print("Calculating Combo_mean_score correlations")
	score_corr = fxns.get_damaging_pred_all(df)
	score_corr.to_csv(f"{args.out}.missense_pred.corr.tsv", sep='\t', index=False)

	# Write Combo_mean_score's to file
	df.to_csv(f"{args.out}.missense_pred.all.tsv.gz", sep='\t', index=False, compression='gzip')

	print("Execution complete")

if __name__ == "__main__":
	main()
	os._exit(0)
