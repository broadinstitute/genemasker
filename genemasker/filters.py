import pandas as pd

# numeric variant filters
def CADD_phred_20(row):
	return(row['CADD_phred_hg19'] >= 20)

def REVEL_score_0_55(row):
	return(row['REVEL_score_max'] >= 0.55)

def REVEL_score_0_6(row):
	return(row['REVEL_score_max'] >= 0.6)

def REVEL_score_0_75(row):
	return(row['REVEL_score_max'] >= 0.75)

def combo_ic25(row):
	return(row['Combo_mean_score_ic'] >= 0.25)

def combo_og25(row):
	return(row['Combo_mean_score_og'] >= 0.25)

def combo_og50(row):
	return(row['Combo_mean_score_og'] >= 0.5)

def maf0_01(row):
	return(row['MAF'] < 0.0001)

def maf0_1(row):
	return(row['MAF'] < 0.001)

def maf1(row):
	return(row['MAF'] < 0.01)

def maf10(row):
	return(row['MAF'] < 0.1)

def maf_gt_0_5(row):
	return(row['MAF'] < 0.005)

# categorical variant filters
def ClinVar_pred_P_or_LP(row):
	return(any(x in row['clinvar_clnsig'] for x in ['Pathogenic', 'Likely_pathogenic']))

def LoF_HC(row):
	return(any(x in row['LoF'] for x in ['HC']))

def LRT_pred_D(row):
	return(any(x in row['LRT_pred'] for x in ['D']))

def MutationTaster_pred_D_or_A(row):
	return(any(x in row['MutationTaster_pred'] for x in ['D','A']))

def Polyphen2_HDIV_pred_D(row):
	return(any(x in row['Polyphen2_HDIV_pred'] for x in ['D']))

def Polyphen2_HDIV_pred_D_or_P(row):
	return(any(x in row['Polyphen2_HDIV_pred'] for x in ['D','P']))

def Polyphen2_HDIV_pred_P(row):
	return(any(x in row['Polyphen2_HDIV_pred'] for x in ['P']))

def Polyphen2_HVAR_pred_D(row):
	return(any(x in row['Polyphen2_HVAR_pred'] for x in ['D']))

def Polyphen2_HVAR_pred_D_or_P(row):
	return(any(x in row['Polyphen2_HVAR_pred'] for x in ['D','P']))
 
def Polyphen2_HVAR_pred_P(row):
	return(any(x in row['Polyphen2_HVAR_pred'] for x in ['P']))

def SIFT_pred_D(row):
	return(any(x in row['SIFT_pred'] for x in ['D']))

def essential_splice(row):
	return(any(x in row['Consequence'] for x in ['splice_acceptor_variant','splice_donor_variant']))

def frameshift(row):
	return(any(x in row['Consequence'] for x in ['frameshift_variant']))

def incomplete_terminal_codon(row):
	return(any(x in row['Consequence'] for x in ['incomplete_terminal_codon_variant']))

def indels(row):
	return(any(x in row['Consequence'] for x in ['inframe_insertion','inframe_deletion']))

def missense(row):
	return(any(x in row['Consequence'] for x in ['missense_variant']))

def splice_region_variant(row):
	return(any(x in row['Consequence'] for x in ['splice_region_variant']))

def start_retained(row):
	return(any(x in row['Consequence'] for x in ['start_retained_variant']))

def stop_gained(row):
	return(any(x in row['Consequence'] for x in ['stop_gained']))

def stop_lost(row):
	return(any(x in row['Consequence'] for x in ['stop_lost']))

def stop_retained(row):
	return(any(x in row['Consequence'] for x in ['stop_retained_variant']))

def synonymous(row):
	return(any(x in row['Consequence'] for x in ['synonymous_variant']))

def IMPACT_HIGH(row):
	return(any(x in row['IMPACT'] for x in ['HIGH']))

def IMPACT_MODERATE(row):
	return(any(x in row['IMPACT'] for x in ['MODERATE']))

# masks
def new_damaging_ic25(df):
	return(df.apply(lambda x: LoF_HC(x) or (missense(x) and combo_ic25(x)), axis=1))

def new_damaging_og25(df):
	return(df.apply(lambda x: LoF_HC(x) or (missense(x) and combo_og25(x)), axis=1))

def new_damaging_og25_0_01(df):
	return(df.apply(lambda x: (LoF_HC(x) or (missense(x) and combo_og25(x))) and maf1(x), axis=1))

def new_damaging_og50(df):
	return(df.apply(lambda x: LoF_HC(x) or (missense(x) and combo_og50(x)), axis=1))

def new_damaging_og50_0_01(df):
	return(df.apply(lambda x: (LoF_HC(x) or (missense(x) and combo_og50(x))) and maf1(x), axis=1))

def x23633568_m1(df):
	return(df.apply(lambda x: stop_gained(x) or frameshift(x) or (missense(x) and (SIFT_pred_D(x) or Polyphen2_HVAR_pred_D_or_P(x) or Polyphen2_HDIV_pred_D_or_P(x))), axis=1))

def x24507775_m6_0_01(df):
	return(df.apply(lambda x: (stop_gained(x) or stop_lost(x) or essential_splice(x) or Polyphen2_HDIV_pred_D(x) or Polyphen2_HVAR_pred_D(x)) & maf1(x), axis=1))

def x29177435_m1(df):
	return(df.apply(lambda x: (incomplete_terminal_codon(x) or start_retained(x) or stop_retained(x) or synonymous(x) or IMPACT_MODERATE(x) or IMPACT_HIGH(x)) & maf1(x), axis=1))

def x29378355_m1_0_01(df):
	return(df.apply(lambda x: (stop_gained(x) | missense(x) | essential_splice(x)) & maf1(x) & maf_gt_0_5(x), axis=1))

def x30269813_m4(df):
	return(df.apply(lambda x: (stop_gained(x) | essential_splice(x) | frameshift(x) | indels(x) | (missense(x) & (Polyphen2_HDIV_pred_D(x) | Polyphen2_HVAR_pred_D(x)))) & maf0_1(x), axis=1))

def x30828346_m1(df):
	return(df.apply(lambda x: (missense(x) | synonymous(x) | stop_gained(x) | stop_lost(x)) & maf10(x), axis=1))

def x31118516_m5_0_001(df):
	return(df.apply(lambda x: (LoF_HC(x) | (Polyphen2_HDIV_pred_D(x) & Polyphen2_HVAR_pred_D(x) & SIFT_pred_D(x) & LRT_pred_D(x) & MutationTaster_pred_D_or_A(x))) & maf0_1(x), axis=1))

def x31383942_m10(df):
	return(df.apply(lambda x: (stop_gained(x) | stop_lost(x) | frameshift(x) | essential_splice(x) | ClinVar_pred_P_or_LP(x)) & maf0_1(x), axis=1))

def x31383942_m4(df):
	return(df.apply(lambda x: (stop_gained(x) | stop_lost(x) | frameshift(x) | essential_splice(x) | (missense(x) & REVEL_score_0_55(x))) & maf0_1(x), axis=1))

def x32141622_m4(df):
	return(df.apply(lambda x: indels(x), axis=1))

def x32141622_m7(df):
	return(df.apply(lambda x: stop_gained(x) | essential_splice(x) | frameshift(x) | splice_region_variant(x), axis=1))

def x32141622_m7_0_01(df):
	return(df.apply(lambda x: (stop_gained(x) | essential_splice(x) | frameshift(x) | splice_region_variant(x)) & maf1(x), axis=1))

def x32853339_m1(df):
	return(df.apply(lambda x: (IMPACT_HIGH(x) | ClinVar_pred_P_or_LP(x)) & maf1(x), axis=1))

def x34183866_m1(df):
	return(df.apply(lambda x: (stop_gained(x) | essential_splice(x) | frameshift(x) | indels(x) | (missense(x) & Polyphen2_HDIV_pred_D(x))) & maf0_01(x), axis=1))

def x34216101_m3_0_001(df):
	return(df.apply(lambda x: (Polyphen2_HDIV_pred_P(x) | Polyphen2_HVAR_pred_P(x)) & maf0_1(x), axis=1))

def x36327219_m3(df):
	return(df.apply(lambda x: (IMPACT_HIGH(x) | (indels(x) & IMPACT_MODERATE(x)) | (missense(x) & Polyphen2_HDIV_pred_D(x) & Polyphen2_HVAR_pred_D(x) & LRT_pred_D(x) & MutationTaster_pred_D_or_A(x) & SIFT_pred_D(x))) & maf1(x), axis=1))

def x36411364_m4_0_001(df):
	return(df.apply(lambda x: (LoF_HC(x) | (missense(x) & REVEL_score_0_75(x))) & maf0_1(x), axis=1))

def x37348876_m8(df):
	return(df.apply(lambda x: (CADD_phred_20(x) | LoF_HC(x)) & maf0_1(x), axis=1))
