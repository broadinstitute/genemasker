import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.decomposition import IncrementalPCA
from genemasker.tracking import resource_tracker
from genemasker.definitions import annot_na_values
import genemasker.filters as filters
import math
import os
import dask.dataframe as dd
import glob
import argparse
import genemasker.config as config

args = config.args
logger = config.logger
logger_handler = config.logger_handler

def count_notna(row):
	return len([x for x in row if not math.isnan(x)])

def get_max(x):
	if x == "-":
		result = np.nan
	else:
		y=[a for a in str(x).split(",") if a not in ["","."]]
		if len(y) == 0:
			result = np.nan
		else:
			result = max([float(b) for b in y])
	return result

@resource_tracker(logger)
def process_annotation_file(annot, cols, maf, out_dir, chunk_size=None, n_partitions = 1):
	compr = 'gzip' if annot.endswith(".bgz") else 'infer'
	rankscore_cols = [col for col in cols.keys() if col.endswith('rankscore') or col.endswith('rankscore_hg19')]
	logger.info(f"Processing annotation file: {annot} (compression={compr})")
	iter_csv = pd.read_csv(
		annot,
		iterator=True,
		chunksize=chunk_size,
		compression=compr,
		sep="\t",
		na_values=annot_na_values,
		usecols=list(cols.keys()),
		dtype = cols
	)
	raw_variant_count = 0
	stored_variant_count = 0
	chunk_count = 0
	chunk_paths = []
	for chunk in iter_csv:
		chunk_count += 1
		logger.info(f"  Reading chunk {chunk_count} from {annot}")
		raw_variant_count = raw_variant_count + chunk.shape[0]
		chunk = chunk[(chunk['Consequence'].str.contains("missense_variant", na=False)) & (chunk['PICK'] == "1")]
		chunk = chunk.dropna(how='all', subset=rankscore_cols)
		chunk['REVEL_score_max'] = chunk['REVEL_score'].apply(lambda x: get_max(x))
		chunk = chunk.merge(maf.loc[maf['#Uploaded_variation'].isin(chunk['#Uploaded_variation'])], on="#Uploaded_variation", how="left")
		stored_variant_count = stored_variant_count + chunk.shape[0]
		chunk = dd.from_pandas(chunk, npartitions=n_partitions)
		dd.to_parquet(chunk, f"{out_dir}/chunk{chunk_count}", write_metadata_file=True, name_function = lambda n: f"{os.path.basename(annot)}-{chunk_count}-{n}.parquet")
		chunk_paths.extend(glob.glob(f"{out_dir}/chunk{chunk_count}/*.parquet"))
	logger.info(f"Finished processing annotation file: {annot}")
	return chunk_paths, chunk_count, raw_variant_count, stored_variant_count, rankscore_cols

@resource_tracker(logger)
def filter_annotation_file(chunk_paths, cols, ids):
	i = 0
	fdfs = []
	n = 0
	m = 0
	for p in chunk_paths:
		i = i + 1
		df = pd.read_parquet(p, columns = cols)
		m = m + df.shape[0]
		df = df[df["#Uploaded_variation"].isin(ids)]
		fdfs = fdfs + [df]
		n = n + df.shape[0]
		logger.info(f"  ({i}/{len(chunk_paths)}) filter {os.path.basename(p)}: {df.shape[0]}[{n}/{m}]")
	return pd.concat(fdfs, sort=False, ignore_index=True) if fdfs else pd.DataFrame()

@resource_tracker(logger)
def impute_annotation_file(chunk_paths, iter_imp, rankscore_cols, means = None, stddevs = None):
	if (means is not None and stddevs is None) or (means is None and stddevs is not None):
		raise ValueError("impute_annotation_file(): means and stddevs must both be passed together")
	chunk_paths_out = []
	i = 0
	for p in chunk_paths:
		i = i + 1
		o = p.replace(".parquet",".imputed.parquet")
		df = pd.read_parquet(p)
		if means is not None:
			# mean and center with respect to full dataset
			df_center_scale = df.copy()
			df_center_scale[rankscore_cols] = (df_center_scale[rankscore_cols] - means) / stddevs
			df[rankscore_cols] = iter_imp.transform(df_center_scale[rankscore_cols])
		else:
			df[rankscore_cols] = iter_imp.transform(df[rankscore_cols])
		df.to_parquet(o, engine="pyarrow", index=False)
		#os.remove(p)
		logger.info(f"  ({i}/{len(chunk_paths)}) {os.path.basename(p)} -> {os.path.basename(o)}")
		chunk_paths_out.append(o)
	return chunk_paths_out

@resource_tracker(logger)
def fit_incremental_pca(chunk_paths, rankscore_cols, means, stddevs):
	pca = IncrementalPCA(n_components = len(rankscore_cols))
	i = 0
	for p in chunk_paths:
		i = i + 1
		df = pd.read_parquet(p)
		# mean and center with respect to full dataset
		df[rankscore_cols] = (df[rankscore_cols] - means) / stddevs
		pca.partial_fit(df[rankscore_cols])
		logger.info(f"  ({i}/{len(chunk_paths)}) incremental pca partial fit on {os.path.basename(p)} complete")
	return pca

@resource_tracker(logger)
def calculate_pca_ica_scores(chunk_paths, pca_fit, pca_n, ica_fit, rankscore_cols, ica_cols, pca_means = None, pca_stddevs = None, ica_means = None, ica_stddevs = None):
	if (pca_means is not None and pca_stddevs is None) or (pca_means is None and pca_stddevs is not None):
		raise ValueError("calculate_pca_ica_scores(): pca_means and pca_stddevs must both be passed together")
	if (ica_means is not None and ica_stddevs is None) or (ica_means is None and ica_stddevs is not None):
		raise ValueError("calculate_pca_ica_scores(): ica_means and ica_stddevs must both be passed together")
	chunk_paths_out = []
	i = 0
	for p in chunk_paths:
		i = i + 1
		o = p.replace(".parquet",".scores.parquet")
		df = pd.read_parquet(p)
		if pca_means is not None:
			# mean and center with respect to full dataset
			df_pca = df.copy()
			df_pca[rankscore_cols] = (df_pca[rankscore_cols] - pca_means) / pca_stddevs
			df_trans = pca_fit.transform(df_pca[rankscore_cols])
		else:
			df_trans = pca_fit.transform(df[rankscore_cols])
		df = pd.concat([df, pd.DataFrame(data = df_trans, columns = [f"pc{i+1}" for i in range(pca_n)])], axis = 1)
		if ica_means is not None:
			# mean and center with respect to full dataset
			df_ica = df.copy()
			df_ica[ica_cols] = (df_ica[ica_cols] - pca_means) / pca_stddevs
			df_trans = ica_fit.transform(df_ica[ica_cols])
		else:
			df_trans = ica_fit.transform(df[ica_cols])
		df = pd.concat([df, pd.DataFrame(data = df_trans, columns = [f"ic{i+1}" for i in range(len(ica_cols))])], axis = 1)
		df.to_parquet(o, engine="pyarrow", index=False)
		#os.remove(p)
		logger.info(f"  ({i}/{len(chunk_paths)}) {os.path.basename(p)} -> {os.path.basename(o)}")
		chunk_paths_out.append(o)
	return chunk_paths_out

@resource_tracker(logger)
def calculate_correlations(ddf, cols, maf_col):
	correlations = []
	i = 0
	for c in cols:
		i = i + 1
		logger.info(f"  ({i}/{len(cols)}) corr({c} ~ MAF)")
		df = ddf[[c, maf_col]].compute()
		df = df.dropna(subset=[c, maf_col])
		if len(df) > 0:
			r, p = pearsonr(df[c], df[maf_col])
			correlations.append({'Feature': c, 'Pearsonr': r, 'P': p})
	return pd.DataFrame(correlations)

@resource_tracker(logger)
def calculate_percent_damaging(ddf, pred_cols):
	results = pd.DataFrame({"Algorithm": [], "DamagingProp": []})
	i = 0
	for alg in pred_cols:
		i = i + 1
		logger.info(f"  ({i}/{len(pred_cols)}) %_damage({alg})")
		df_alg = ddf[[alg]].compute()
		df_alg = df_alg[alg].dropna().reset_index()
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

@resource_tracker(logger)
def get_damaging_pred_og(df):
	cols_drop=['#Uploaded_variation','Eigen-PC-raw_coding_rankscore', 'Eigen-raw_coding_rankscore', 'Polyphen2_HVAR_rankscore', 'CADD_raw_rankscore_hg19', 'BayesDel_noAF_rankscore', 'MutPred_rankscore', 'LINSIGHT_rankscore']
	score = (df[[c for c in df.columns if c not in cols_drop]] > 0.67).mean(axis=1, numeric_only=True)
	return score

@resource_tracker(logger)
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

@resource_tracker(logger)
def get_damaging_pred_pca(df, pc_corr, pc_weight):
	#df = df.copy()
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
	pc_weight['PC'] = pc_weight['PC'] + "_rankscore"
	sum_weight = pc_weight['VarianceExplained'].sum()
	pc_weight['VarianceExplained'] = pc_weight['VarianceExplained']/sum_weight
	pc_weight = pc_weight.set_index('PC',drop=True)
	weighted_scores = df_rankscore * pc_weight['VarianceExplained']
	return weighted_scores.sum(axis=1)

@resource_tracker(logger)
def get_damaging_pred_all(ddf):
	df = ddf[['Combo_mean_score_og','Combo_mean_score_pc','Combo_mean_score_ic']].compute()
	corr1, p1 = pearsonr(df['Combo_mean_score_og'], df['Combo_mean_score_ic'])
	corr2, p2 = pearsonr(df['Combo_mean_score_ic'], df['Combo_mean_score_pc'])
	corr3, p3 = pearsonr(df['Combo_mean_score_pc'], df['Combo_mean_score_og'])
	results = pd.DataFrame({'Mean_score1': ['Combo_mean_score_og', 'Combo_mean_score_ic', 'Combo_mean_score_pc'], 'Mean_score2': ['Combo_mean_score_ic', 'Combo_mean_score_pc', 'Combo_mean_score_og'], 'Pearsonr': [corr1, corr2, corr3], 'Pvalue': [p1, p2, p3]})
	return results

@resource_tracker(logger)
def calculate_damage_predictions(chunk_paths, rankscore_cols, ica_cols, pca_n, pc_corr, pca_exp_var, ic_corr):
	chunk_paths_out = []
	i = 0
	for p in chunk_paths:
		i = i + 1
		o = p.replace(".parquet",".predictions.parquet")
		df = pd.read_parquet(p)
		df['Combo_mean_score_og'] = get_damaging_pred_og(df[rankscore_cols])
		df['Combo_mean_score_pc'] = get_damaging_pred_pca(df[[f"pc{i+1}" for i in range(pca_n)]], pc_corr, pca_exp_var)
		df['Combo_mean_score_ic'] = get_damaging_pred_ica(df[[f"ic{i+1}" for i in range(len(ica_cols))]], ic_corr)
		df.to_parquet(o, engine="pyarrow", index=False)
		#os.remove(p)
		logger.info(f"  ({i}/{len(chunk_paths)}) {os.path.basename(p)} -> {os.path.basename(o)}")
		chunk_paths_out.append(o)
	return chunk_paths_out

@resource_tracker(logger)
def calculate_mask_filters(chunk_paths):
	chunk_paths_out = []
	i = 0
	for p in chunk_paths:
		i = i + 1
		o = p.replace(".parquet",".filters.parquet")
		df = pd.read_parquet(p)
		df['new_damaging_ic2'] = filters.new_damaging_ic2(df)
		df['new_damaging_og25'] = filters.new_damaging_og25(df)
		df['new_damaging_og25_0_01'] = filters.new_damaging_og25_0_01(df)
		df['new_damaging_og50'] = filters.new_damaging_og50(df)
		df['new_damaging_og50_0_01'] = filters.new_damaging_og50_0_01(df)
		df['x23633568_m1'] = filters.x23633568_m1(df)
		df['x24507775_m6_0_01'] = filters.x24507775_m6_0_01(df)
		df['x29177435_m1'] = filters.x29177435_m1(df)
		df['x29378355_m1_0_01'] = filters.x29378355_m1_0_01(df)
		df['x30269813_m4'] = filters.x30269813_m4(df)
		df['x30828346_m1'] = filters.x30828346_m1(df)
		df['x31118516_m5_0_001'] = filters.x31118516_m5_0_001(df)
		df['x31383942_m10'] = filters.x31383942_m10(df)
		df['x31383942_m4'] = filters.x31383942_m4(df)
		df['x32141622_m4'] = filters.x32141622_m4(df)
		df['x32141622_m7'] = filters.x32141622_m7(df)
		df['x32141622_m7_0_01'] = filters.x32141622_m7_0_01(df)
		df['x32853339_m1'] = filters.x32853339_m1(df)
		df['x34183866_m1'] = filters.x34183866_m1(df)
		df['x34216101_m3_0_001'] = filters.x34216101_m3_0_001(df)
		df['x36327219_m3'] = filters.x36327219_m3(df)
		df['x36411364_m4_0_001'] = filters.x36411364_m4_0_001(df)
		df['x37348876_m8'] = filters.x37348876_m8(df)
		df.to_parquet(o, engine="pyarrow", index=False)
		#os.remove(p)
		logger.info(f"  ({i}/{len(chunk_paths)}) {os.path.basename(p)} -> {os.path.basename(o)}")
		chunk_paths_out.append(o)
	return chunk_paths_out

@resource_tracker(logger)
def generate_regenie_groupfiles(ddf, out):
	masks = ['new_damaging_ic2','new_damaging_og25','new_damaging_og25_0_01','new_damaging_og50','new_damaging_og50_0_01','x23633568_m1','x24507775_m6_0_01','x29177435_m1','x29378355_m1_0_01','x30269813_m4','x30828346_m1','x31118516_m5_0_001','x31383942_m10','x31383942_m4','x32141622_m4','x32141622_m7','x32141622_m7_0_01','x32853339_m1','x34183866_m1','x34216101_m3_0_001','x36327219_m3','x36411364_m4_0_001','x37348876_m8']
	df = ddf[["#Uploaded_variation",'Gene'] + masks].compute()
	df[['chr','pos','ref','alt']]=df["#Uploaded_variation"].str.split(":",expand=True)
	df['chr_num']=df['chr'].str.replace("chr","").replace({'X': '23', 'Y': '24', 'XY': '25', 'MT': '26', 'M': '26'})
	df.chr_num=df.chr_num.astype(int)
	df.pos=df.pos.astype(int)
	df.sort_values(by=['chr_num','pos'],inplace=True)
	df.reset_index(drop=True, inplace=True)
	genes=df['Gene'].unique()
	logger.info("found " + str(len(genes)) + " genes")
	logger.info("extracting minimum positions for genes")
	df_first_pos=df[['Gene','chr','chr_num','pos']].drop_duplicates(subset=['Gene'], keep='first')
	with open(f"{out}.regenie.setlist.tsv", 'w') as setlist:
		logger.info("grouping variants into genes")
		setlist_df=df[['Gene','#Uploaded_variation']]
		setlist_df=setlist_df.groupby('Gene', as_index=False, sort=False).agg(','.join)
		setlist_df=df_first_pos.merge(setlist_df)
		setlist_df.sort_values(by=['chr_num','pos'],inplace=True)
		setlist_df[['Gene','chr','pos','#Uploaded_variation']].to_csv(setlist, header=False, index=False, sep="\t", na_rep="NA")
	for mask in masks:
		logger.info("adding annotations for mask " + mask)
		mask_df=df[df[mask] == True][['#Uploaded_variation','Gene']]
		mask_df['mask']=mask
		with open(f"{out}.regenie.annotations.{mask}.tsv", 'w') as annots:
			logger.info(f"writing annotations to file for mask {mask}")
			mask_df.to_csv(annots, header=False, index=False, sep="\t", na_rep="NA")
		with open(f"{out}.regenie.mask.{mask}.tsv", 'w') as mask_file:
			logger.info(f"writing mask annotations used for mask {mask}")
			mask_file.write(mask + "\t" + mask)
