import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kendalltau, linregress, rankdata
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import numpy as np
import piflib 
from piflib import compute_cigs
import warnings



# Suppress SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=pd.core.common.SettingWithCopyWarning)

# Activate pandas <-> R DataFrame conversion
pandas2ri.activate()
from rpy2.robjects.conversion import localconverter


# Import the sdcMicro package
sdcMicro = importr('sdcMicro')



def stats(suda, pif, suda_att, k_combined_field):
    
    
    print('\n' + '='*40)
    print('        ROW LEVEL CORRELATION      ')
    print('='*40 + '\n')
    
    
    # Spearman Rank Correlation
    correlation_spearman, p_value = spearmanr(suda['dis-score'], pif['RIG'])
    print(f"Spearman Rank Correlation between suda & pif: {correlation_spearman:.2f}, P-value: {p_value:.4f}")
    
    # Kendall's Tau Correlation
    tau, p_value = kendalltau(suda['dis-score'], pif['RIG'])
    print(f"Kendall's Tau between suda & pif: {tau:.2f}, p-value: {p_value:.4f}")
    
    # Pearson Correlation
    correlation_pearson = suda['dis-score'].corr(pif['RIG'])
    
    #field relationship
    
    sum_df = pif.drop(columns=['RIG']).sum()
    sum_df = pd.DataFrame(list(sum_df.items()), columns=['variable', 'sum'])
    
    
    # Merge the two DataFrames on 'variable'
    merged_field_values = pd.merge(suda_att, sum_df, on='variable')
    merged_field_values = pd.merge(merged_field_values, k_combined_field, on='variable') 
    


    # Calculate the correlation between 'contribution' and 'Sum'
    correlation_pif_suda = merged_field_values['contribution'].corr(merged_field_values['sum'])
    correlation_pif_k = merged_field_values['Normalized Difference'].corr(merged_field_values['sum'])
    correlation_suda_k = merged_field_values['contribution'].corr(merged_field_values['Normalized Difference'])

    
    
    
    
    # Plot with trend line
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.scatter(suda['dis-score'], pif['RIG'], alpha=0.3, label='Data Points')
    
    # Fit and plot the trend line
    slope, intercept, _, _, _ = linregress(suda['dis-score'], pif['RIG'])
    ax.plot(suda['dis-score'], slope * suda['dis-score'] + intercept, color='red', label='Trend Line')
    
    ax.set_title(f'Pearson Correlation: {correlation_pearson:.2f}')
    ax.set_xlabel('suda["dis-score"]')
    ax.set_ylabel('pif["RIG"]')
    ax.grid(False)
    ax.legend()
    
    plt.show()    
    plt.close(fig)
    print('____________________________________________________________________________________')
    
    print('\n' + '='*40)
    print('      FIELD LEVEL CORRELATIONS     ')
    print('='*40 + '\n')
    
    result_text_field_s_p = f"Pearson Correlation between SUDA & PIF field level: {correlation_pif_suda:.2f}"
    result_text_field_p_k = f"Pearson Correlation between K-combined & PIF field level: {correlation_pif_k:.2f}"
    result_text_field_s_k = f"Pearson Correlation between SUDA & K-combined field level: {correlation_suda_k:.2f}"

    print(result_text_field_s_p)
    print(result_text_field_p_k)
    print(result_text_field_s_k)
    
    print('\n')
    
    
    # Spearman correlation between PIF & SUDA
    correlation_spearman_pif_suda, p_value = spearmanr(merged_field_values['contribution'], merged_field_values['sum'])
    print(f"Spearman Rank Correlation between PIF & SUDA: {correlation_spearman_pif_suda:.2f}, P-value: {p_value:.4f}")

    # Spearman correlation between PIF & K
    correlation_spearman_pif_k, p_value = spearmanr(merged_field_values['Normalized Difference'], merged_field_values['sum'])
    print(f"Spearman Rank Correlation between PIF & K: {correlation_spearman_pif_k:.2f}, P-value: {p_value:.4f}")

    # Spearman correlation between SUDA & K
    correlation_spearman_suda_k, p_value = spearmanr(merged_field_values['contribution'], merged_field_values['Normalized Difference'])
    print(f"Spearman Rank Correlation between SUDA & K: {correlation_spearman_suda_k:.2f}, P-value: {p_value:.4f}")

    return 




def convert_to_numeric(df):
    for col in df.select_dtypes(include=['object']).columns:
        df[col] = df[col].astype('category').cat.codes
    return df



def compute_suda2(dataframe, sample_fraction=0.2, missing_value=np.NaN):
    with localconverter(pandas2ri.converter):
        r_df = robjects.conversion.py2rpy(dataframe)
        
    
    convert_to_numeric(dataframe) 
    
    suda_result = sdcMicro.suda2(r_df, missing=missing_value, DisFraction=sample_fraction)
    contribution_percent = list(suda_result.rx2('contributionPercent'))
    score = list(suda_result.rx2('score'))
    dis_score = list(suda_result.rx2('disScore'))

    return (contribution_percent, score, dis_score,suda_result)





import piflib.pif_calculator as pif
def calculate_summed_dis_scores(k_combined_all, AOMIC, sample_fraction=0.3, missing_value= -999):

    dis_scores = []
    pif_score = []
    
    k_combined_all = k_combined_all[k_combined_all['Combination'].str.count(',') >= 2]
    
    for index, row in k_combined_all.iterrows():
        combination_fields = row['Combination'].split(', ')
       
        df_subset = AOMIC[combination_fields]
        df_subset = df_subset.fillna(-999)
        
        contribution_percent, score, dis_score, attribute_contributions = compute_suda2(
            df_subset, sample_fraction=sample_fraction, missing_value=missing_value)  
        
        summed_dis_score = sum(dis_score)
        dis_scores.append(summed_dis_score)
        
        cigs = piflib.compute_cigs(df_subset)
        pif_95 = pif.compute_pif(cigs, 0.95)
        pif_score.append(pif_95)
        

    k_combined_all['pif_score'] = pif_score
    k_combined_all['suda_score'] = dis_scores
    k_combined_all.fillna(0, inplace=True)
    
    return k_combined_all





def plot_calc(k_combined_all): 
    k_combined_all = k_combined_all.replace([np.inf, -np.inf], 0)

    # SUDA vs. K-Combined
    spearman_corr, _ = spearmanr(k_combined_all['Score'], k_combined_all['suda_score'])
    pearson_corr = k_combined_all['Score'].corr(k_combined_all['suda_score'])
    print(f"Spearman Correlation between Suda sum and K-combined: {spearman_corr:.2f}")
    print(f"Pearson Correlation between Suda sum and K-combined: {pearson_corr:.2f}")
    
    plt.figure(figsize=(8, 4))
    plt.scatter(rankdata(k_combined_all['Score']), rankdata(k_combined_all['suda_score']), alpha=0.6)
    plt.title(f'Scatter Plot of Ranked k-combined Score vs. SUDA Score\n(Spearman Correlation: {spearman_corr:.2f})')
    plt.xlabel('Ranked K-combined Score')
    plt.ylabel('Ranked SUDA Score')
    plt.grid(False)
    plt.show()
    
    print('________________________________________________________________________')
    print('\n')
    
    # PIF vs. K-Combined
    spearman_corr, _ = spearmanr(k_combined_all['Score'], k_combined_all['pif_score'])
    pearson_corr = k_combined_all['Score'].corr(k_combined_all['pif_score'])
    print(f"Spearman Correlation between PIF 95% and K-combined: {spearman_corr:.2f}")
    print(f"Pearson Correlation between PIF 95% and K-combined: {pearson_corr:.2f}")
    
    plt.figure(figsize=(8, 4))
    plt.scatter(rankdata(k_combined_all['Score']), rankdata(k_combined_all['pif_score']), alpha=0.6)
    plt.title(f'Scatter Plot of Ranked k-combined Score vs. PIF Score\n(Spearman Correlation: {spearman_corr:.2f})')
    plt.xlabel('Ranked K-combined Score')
    plt.ylabel('Ranked PIF Score')
    plt.grid(False)
    plt.show()
    
    print('________________________________________________________________________')
    print('\n')
    
    # PIF vs. SUDA
    spearman_corr, _ = spearmanr(k_combined_all['suda_score'], k_combined_all['pif_score'])
    pearson_corr = k_combined_all['suda_score'].corr(k_combined_all['pif_score'])
    print(f"Spearman Correlation between PIF 95% and SUDA: {spearman_corr:.2f}")
    print(f"Pearson Correlation between PIF 95% and SUDA: {pearson_corr:.2f}")
    
    plt.figure(figsize=(8, 4))
    plt.scatter(rankdata(k_combined_all['suda_score']), rankdata(k_combined_all['pif_score']), alpha=0.6)
    plt.title(f'Scatter Plot of Ranked SUDA Score vs. PIF Score\n(Spearman Correlation: {spearman_corr:.2f})')
    plt.xlabel('Ranked SUDA Score')
    plt.ylabel('Ranked PIF Score')
    plt.grid(False)
    plt.show()
    
    return





def rst_outlier_case2(data, column, k=2.2414):
    column_data = data[column]
    column_data_array = np.array(column_data)
    median = np.nanmedian(column_data_array)
    mad = np.nanmedian(np.abs(column_data_array - median))
    madn = mad / 0.6745

    z_scores = (column_data_array - median) / madn
    class_outliers = (np.abs(z_scores) > k).astype(int)
    outlier_indices = data[class_outliers == 1].index.tolist()

    above_outliers = (z_scores > k).astype(int)
    above_outlier_indices = data[above_outliers == 1].index.tolist()

    

    return class_outliers, madn, mad, outlier_indices, above_outlier_indices



def plot_boxplot_with_outliers(data, class_outliers, median, madn, k=2.2414):
    data = np.array(data)
    lower_whisker = median - k * madn
    upper_whisker = median + k * madn
    fig, ax = plt.subplots()
    ax.bxp([{
        'med': median,
        'q1': np.percentile(data, 25),
        'q3': np.percentile(data, 75),
        'whislo': lower_whisker,
        'whishi': upper_whisker
    }], vert=False, showfliers=False, patch_artist=True, boxprops=dict(facecolor='lightblue'))
    outliers = data[class_outliers == 1]
    ax.scatter(outliers, [1] * len(outliers), color='red', zorder=3, label='Outliers')
    ax.set_title('Boxplot with MAD-based Whiskers and Outliers')
    ax.set_xlabel('Value')
    ax.legend()
    plt.show()


