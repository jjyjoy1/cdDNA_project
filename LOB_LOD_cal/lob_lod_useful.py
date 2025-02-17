import glob
import pandas as pd
import numpy as np
import re
from scipy.stats import norm
import statsmodels.api as sm


def merged_table(filenames):
    merged_data = pd.DataFrame()
    for filename in filenames:
        data = pd.read_csv(filename, sep='\t')
        merged_data = pd.concat([merged_data, data])
    return merged_data


gene_to_chr = {
    'AKT1': '14',
    'TP53': '17',
    'ERBB2' : '17',
    'ESR1' : '6',
    'KRAS' : '12',
    'PIK3CA' : '3'
}

#Convert mutant calling files to matrix table
def files2mdata_calltab(fls):
    mdata = merged_table(fls)
    mdata['Chr'] = mdata['Gene'].apply(lambda x: gene_to_chr.get(x, 'Unknown'))
    mdata = mdata[mdata['MutantFraction']<0.2].fillna(0)
    mdata['Chr:hg19Pos:CDSChange'] = mdata['Chr'].astype(str) + ':' + mdata['hg19Pos'].astype(str) + ':' + mdata['CDSChange'].astype(str)
    mdata['MutantMolecules'] = mdata['MutantMolecules'].astype(float)
    mdata['#UIDs/Amplicon'] = mdata['#UIDs/Amplicon'].astype(int)
    mdata['GE'] = mdata['GE'].astype(int)
    mdata['VD'] = (mdata['MutantMolecules'] * mdata['#UIDs/Amplicon']) / mdata['GE']
    return mdata


def pivottable1(merged_data, field):
    pivot_table = merged_data.pivot_table(index='SampleId', columns='Chr:hg19Pos:CDSChange', values=field, fill_value=0)
    return pivot_table

def pivottable2(merged_data, field):
    pivot_table = merged_data.pivot_table(index='SampleId', columns='CDSChange', values=field, fill_value=0)
    return pivot_table

#Remove all zero columns, which means using all potential SafeSeq calling variants LOB calculation
pdata_lob = pdata_lob.loc[:, (pdata_lob != 0).any(axis=0)]
pdata_lob.shape
#(40, 1214), global lob99  will be in position 12~13
lob99 = pdata_lob.apply(lambda x: np.percentile(x, 99))
lob99.sort_values(ascending=False).head(20)

flat_series = lob99.values.flatten()
lob99_global = np.percentile(flat_series, 99)

###Calculate lob function
def lob_cal(pdata):
    lob99 = pdata_lob.apply(lambda x: np.percentile(x, 99))
    print (lob99.sort_values(ascending=False).head(20))
    flat_series = lob99.values.flatten()
    lob99_global = np.percentile(flat_series, 99)
    print (lob99_global)
    return lob99, lob99_global

###Prepare binary table for probit lod assay
df_reversed = pdata_lod_mc
#df_reversed = pdata_lod_mc['c.2324insATACGTGATGGC']    #can select different table
df_binary_reversed = pd.DataFrame().reindex_like(df_reversed)
for col in df_reversed.columns:
    df_binary_reversed[col] = (df_reversed[col] > Spick_lob99[col]).astype(int)


df_binary_reversed_reset = df_binary_reversed.reset_index()
long_binary_df = df_binary_reversed_reset.melt(id_vars=['MM_level'], var_name='Gene', value_name='Outcome')

LOD95_with_CI_results_corrected = {}

unique_genes = long_binary_df['Gene'].unique()
# Loop through each unique gene to fit the Probit model and calculate LOD95 along with its 95% CI
for gene in unique_genes:
    subset_df = long_binary_df[long_binary_df['Gene'] == gene]
    endog = subset_df['Outcome'].values
    exog = sm.add_constant(subset_df['MM_level'])  # Adding a constant term for intercept
    probit_model = sm.Probit(endog, exog)
    try:
        probit_result = probit_model.fit(disp=0)  # disp=0 suppresses output
        # Calculate the LOD at 95% positive outcome
        beta_0 = probit_result.params[0]
        beta_1 = probit_result.params[1]
        LOD95 = (norm.ppf(0.95) - beta_0) / beta_1
        # Extract the covariance matrix
        cov_matrix = probit_result.cov_params()
        var_beta_0 = cov_matrix.loc['const', 'const']
        var_beta_1 = cov_matrix.loc['MM_level', 'MM_level']
        cov_beta_0_1 = cov_matrix.loc['const', 'MM_level']
        # Calculate the variance for LOD95
        var_LOD95 = ((1 / beta_1)**2 * var_beta_0) + ((-beta_0 / beta_1**2)**2 * var_beta_1) - \
                    (2 * (1 / beta_1) * (-beta_0 / beta_1**2) * cov_beta_0_1)
        # Calculate the 95% confidence interval for LOD95
        lower_bound_LOD95 = LOD95 - 1.96 * (var_LOD95**0.5)
        upper_bound_LOD95 = LOD95 + 1.96 * (var_LOD95**0.5)
        LOD95_with_CI_results_corrected[gene] = {
            'LOD95': LOD95,
            'Lower_Bound_95CI': lower_bound_LOD95,
            'Upper_Bound_95CI': upper_bound_LOD95
        }
    except Exception as e:
        LOD95_with_CI_results_corrected[gene] = {'Error': str(e)}

###Lod dictionary
LOD95_with_CI_results_corrected


###Convert to table
LOD95_with_CI_results_df = pd.DataFrame.from_dict(LOD95_with_CI_results_corrected, orient='index')
LOD95_with_CI_results_csv_path = '/scratch/raw_data/LOD95_with_95CI_results.csv'
LOD95_with_CI_results_df.to_csv(LOD95_with_CI_results_csv_path)

#data = pdata_lod_mc
def pdata2lod(data):
    df_reversed = data
    df_binary_reversed = pd.DataFrame().reindex_like(df_reversed)
    for col in df_reversed.columns:
        df_binary_reversed[col] = (df_reversed[col] > Spick_lob99[col]).astype(int)
    df_binary_reversed_reset = df_binary_reversed.reset_index()
    long_binary_df = df_binary_reversed_reset.melt(id_vars=['MM_level'], var_name='Gene', value_name='Outcome')
    LOD95_with_CI_results_corrected = {}
    unique_genes = long_binary_df['Gene'].unique()
    # Loop through each unique gene to fit the Probit model and calculate LOD95 along with its 95% CI
    for gene in unique_genes:
        subset_df = long_binary_df[long_binary_df['Gene'] == gene]
        endog = subset_df['Outcome'].values
        exog = sm.add_constant(subset_df['MM_level'])  # Adding a constant term for intercept
        probit_model = sm.Probit(endog, exog)
        try:
            probit_result = probit_model.fit(disp=0)  # disp=0 suppresses output
        # Calculate the LOD at 95% positive outcome
            beta_0 = probit_result.params[0]
            beta_1 = probit_result.params[1]
            LOD95 = (norm.ppf(0.95) - beta_0) / beta_1
        # Extract the covariance matrix
            cov_matrix = probit_result.cov_params()
            var_beta_0 = cov_matrix.loc['const', 'const']
            var_beta_1 = cov_matrix.loc['MM_level', 'MM_level']
            cov_beta_0_1 = cov_matrix.loc['const', 'MM_level']
        # Calculate the variance for LOD95
            var_LOD95 = ((1 / beta_1)**2 * var_beta_0) + ((-beta_0 / beta_1**2)**2 * var_beta_1) - \
                    (2 * (1 / beta_1) * (-beta_0 / beta_1**2) * cov_beta_0_1)
        # Calculate the 95% confidence interval for LOD95
            lower_bound_LOD95 = LOD95 - 1.96 * (var_LOD95**0.5)
            upper_bound_LOD95 = LOD95 + 1.96 * (var_LOD95**0.5)
            LOD95_with_CI_results_corrected[gene] = {
                'LOD95': LOD95,
                'Lower_Bound_95CI': lower_bound_LOD95,
                'Upper_Bound_95CI': upper_bound_LOD95
            }
        except Exception as e:
            LOD95_with_CI_results_corrected[gene] = {'Error': str(e)}
    results_df = pd.DataFrame.from_dict(LOD95_with_CI_results_corrected, orient='index')
    return results_df




