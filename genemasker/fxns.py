import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from genemasker.tracking import resource_tracker
import math
import os
import dask.dataframe as dd

def get_max(x):
	if x == "NaN":
		result = np.nan
	else:
		y=[a for a in str(x).split(",") if a not in ["","."]]
		if len(y) == 0:
			result = np.nan
		else:
			result = max([float(b) for b in y])
	return result

@resource_tracker
def process_annotation_file(annot, cols, out_dir, chunk_size=100000, n_partitions = 10):

	def partition_name(n):
		return f"{os.path.basename(annot)}-{n}.parquet"

	compr = 'gzip' if annot.endswith(".bgz") else 'infer'
	print(f"Processing annotation file: {annot} (compression={compr})")
	iter_csv = pd.read_csv(
		annot,
		iterator=True,
		chunksize=chunk_size,
		compression=compr,
		sep="\t",
		na_values='-',
		usecols=list(cols.keys()),
		dtype = cols
	)
	chunk_count = 0
	chunk_paths = []
	for chunk in iter_csv:
		chunk_count += 1
		print(f"  Reading chunk {chunk_count} from {annot}")
		chunk = chunk[(chunk['Consequence'].str.contains("missense_variant", na=False)) & (chunk['PICK'] == "1")]
		chunk['REVEL_score_max'] = chunk['REVEL_score'].apply(lambda x: get_max(x))
		chunk = dd.from_pandas(chunk, npartitions=n_partitions)
		dd.to_parquet(chunk, f"{out_dir}/chunk{chunk_count}", write_metadata_file=True, name_function = partition_name)
		chunk_paths.append(f"{out_dir}/chunk{chunk_count}/*.parquet")
	print(f"Finished processing annotation file: {annot}")
	return chunk_paths

@resource_tracker
def calculate_correlations(df, cols, maf_col):
	print("Calculating correlations")
	correlations = []
	for c in cols:
		print(f"  Calculating correlation for {c}")
		valid_data = df.dropna(subset=[c, maf_col])
		if len(valid_data) > 0:
			r, p = pearsonr(valid_data[c], valid_data[maf_col])
			correlations.append({'Feature': c, 'Pearsonr': r, 'P': p})
	print("Finished calculating correlations")
	return pd.DataFrame(correlations)

@resource_tracker
def calculate_percent_damaging(df):
	print("Calculating proportion of damage for each prediction algorithm in the annot file")
	pred_columns = []
	for col in list(df.columns):
		if "_pred" in col:
			pred_columns.append(col)
	results = pd.DataFrame({"Algorithm": [], "DamagingProp": []})
	for alg in pred_columns:
		print(f"  Calculating proportion of damage for {alg}")
		df_alg = df[alg].dropna().reset_index()
		if alg == "MutationTaster_pred":
			prop = len(df_alg[(df_alg[alg].str.contains("D")) | (df_alg[alg].str.contains("A"))]) / len(df_alg)
		elif alg == "MutationAssessor_pred":
			prop = len(df_alg[(df_alg[alg].str.contains("H"))]) / len(df_alg)
		elif alg == "Aloft_pred":
			prop = len(df_alg[(df_alg[alg].str.contains("Recessive")) | (df_alg[alg].str.contains("Dominant"))]) / len(df_alg)
		else:
			prop = len(df_alg[(df_alg[alg].str.contains("D"))]) / len(df_alg)
		results = pd.concat([results, pd.DataFrame({"Algorithm": [alg], "DamagingProp": [prop]})], ignore_index=True)
	results = results.sort_values("DamagingProp")
	results = pd.concat([results, pd.DataFrame({"Algorithm": ["Mean"], "DamagingProp": [results["DamagingProp"].mean()]})], ignore_index=True)
	return results

def damaging_pred(row):
	if (math.isnan(row)):
		return np.nan 
	elif row > 0.67:
		return 1
	else:
		return 0

@resource_tracker
def get_damaging_pred_og(df):
	cols_drop=['#Uploaded_variation','Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore', 'MutPred_rankscore', 'LINSIGHT_rankscore']
	score = (df[[c for c in df.columns if c not in cols_drop]] > 0.67).mean(axis=1, numeric_only=True)
	return score

@resource_tracker
def get_damaging_pred_ica(df, ic_corr):
	ic_corr['Dir'] = np.nan
	ic_corr.loc[ic_corr['Pearsonr']>0, 'Dir'] = -1.0
	ic_corr.loc[ic_corr['Pearsonr']<0, 'Dir'] = 1.0
	ic_corr.loc[ic_corr['P']>=0.01, 'Dir'] = 0.0
	for col in df.columns:
		if ic_corr.loc[ic_corr['Feature']==col]['Dir'].iloc[0] == 0:
			df = df.drop(columns=col)
		else:
			df.loc[:,col] = df[col]*ic_corr.loc[ic_corr['Feature']==col]['Dir'].iloc[0]
	df_rankscore = df.rank(method="min", ascending=True, pct=True)
	df_rankscore.columns = [f'{col}_rankscore' for col in list(df_rankscore.columns)]
	score = (df_rankscore > 0.67).mean(axis=1, numeric_only=True)
	return score

def weighted_pc_score(row, weights):
    weights['PC'] = weights['PC'] + "_rankscore"
    row = row.reset_index()
    row.columns = ['PC', 'PC_pred']
    df = pd.merge(row, weights, how='left', on='PC')
    df['Weightted_PC_pred'] = df['PC_pred']*df['VarianceExplained']
    return df['Weightted_PC_pred'].sum()

@resource_tracker
def get_damaging_pred_pca(df, pc_corr, pc_weight):
	df = df.copy()
	pc_corr['Dir'] = np.nan
	pc_corr.loc[pc_corr['Pearsonr']>0, 'Dir'] = -1
	pc_corr.loc[pc_corr['Pearsonr']<0, 'Dir'] = 1
	pc_corr.loc[pc_corr['P']>=0.01, 'Dir'] = 0   
	pc_to_drop = []
	for col in df.columns:
		if pc_corr.loc[pc_corr['Feature']==col]['Dir'].iloc[0] == 0:
			pc_to_drop.append(col)
			df = df.drop(columns=col)
		else:
			df.loc[:,col] = df[col]*pc_corr.loc[pc_corr['Feature']==col]['Dir'].iloc[0]
	df_rankscore = df.rank(method="min", ascending=True, pct=True)
	df_rankscore.columns = [f'{col}_rankscore' for col in list(df_rankscore.columns)]
	for col in list(df_rankscore.columns):
		df_rankscore[col] = df_rankscore[col].apply(damaging_pred)
	pc_weight = pd.concat([pc_weight, pd.DataFrame({'VarianceExplained': [0]})], ignore_index=True).sort_values('VarianceExplained')
	pc_weight['VarianceExplained'] = pc_weight['VarianceExplained'].diff()
	pc_weight = pc_weight.dropna()
	pc_weight['PC'] = [f'pc{i+1}' for i in range(len(pc_weight))]
	pc_weight.loc[pc_weight['PC'].isin(pc_to_drop),'VarianceExplained'] = np.nan
	pc_weight = pc_weight.dropna()
	sum_weight = pc_weight['VarianceExplained'].sum()
	pc_weight['VarianceExplained'] = pc_weight['VarianceExplained']/sum_weight
	score = df_rankscore.apply(weighted_pc_score, args=[pc_weight], axis=1)
	return score

@resource_tracker
def get_damaging_pred_all(df):
	corr1, p1 = pearsonr(df['Combo_mean_score_og'], df['Combo_mean_score_ic'])
	corr2, p2 = pearsonr(df['Combo_mean_score_ic'], df['Combo_mean_score_pc'])
	corr3, p3 = pearsonr(df['Combo_mean_score_pc'], df['Combo_mean_score_og'])
	results = pd.DataFrame({'Mean_score1': ['Combo_mean_score_og', 'Combo_mean_score_ic', 'Combo_mean_score_pc'], 'Mean_score2': ['Combo_mean_score_ic', 'Combo_mean_score_pc', 'Combo_mean_score_og'], 'Pearsonr': [corr1, corr2, corr3], 'Pvalue': [p1, p2, p3]})
	return results
