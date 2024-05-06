import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from typing import Dict, List, Literal
from IPython.display import display

def ssgsea_score(ranks: pd.DataFrame, genes: List, use_old_formula: bool = False) -> pd.Series:
    """
    Calculates single sample GSEA score based on vector of gene expression ranks.
    Only overlapping genes will be analyzed.
    The original article describing the ssGSEA formula: https://doi.org/10.1038/nature08460.
    We use adapted fast function. Result is the same as in analogous packages (like GSVA).

    Note: formula was updated in November 2023.
    Please visit Wiki-page for more details: https://bostongene.atlassian.net/wiki/spaces/BIT/pages/3427991553/ssGSEA

    :param ranks: DataFrame with gene expression ranks; samples in columns and genes in rows
    :param genes: list or set, genes of interest
    :param use_old_formula: (Default: False) if true calculates old ssGSEA-BostonGene Formula
    :return: Series with ssGSEA scores for samples
    """

    # Finding common_genes
    # Note: List is needed here because pandas can not do .loc with sets
    common_genes = list(set(genes).intersection(set(ranks.index)))

    # If not intersections were found
    if not common_genes:
        return pd.Series([0.0] * len(ranks.columns), index=ranks.columns)

    # Ranks of genes inside signature
    sranks = ranks.loc[common_genes]

    if use_old_formula:
        return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (len(ranks.index) - len(common_genes) + 1) / 2

    return (sranks ** 1.25).sum() / (sranks ** 0.25).sum() - (len(ranks.index) + 1) / 2

def ssgsea_formula(expressions: pd.DataFrame, gmt: Dict,
                   rank_method: Literal['average', 'min', 'max', 'first', 'dense'] = 'max',
                   use_old_formula: bool = False) -> pd.DataFrame:
    """
    Return DataFrame with ssGSEA scores for gene signatures from gmt (dict of GeneSets)
    Only overlapping genes will be analyzed
    The original article describing the ssGSEA formula: https://doi.org/10.1038/nature08460

    :param expressions: DataFrame with gene expressions; samples in columns and genes in rows
    :param gmt: keys - signature names, values - bioreactor.gsea.GeneSet
    :param rank_method: {'min', 'max', etc} how to rank genes that have the same expression value
    :param use_old_formula: (Default: False) if true calculates old ssGSEA-BostonGene Formula
    :return: DataFrame with ssGSEA scores, index - signatures, columns - samples
    """

    # Calculate ranks
    ranks = expressions.rank(method=rank_method, na_option='bottom')

    # Calculate ssGSEA scores for all signatures - at the same moment merging it in DataFrame
    return pd.DataFrame(
        {gs_name: ssgsea_score(ranks, geneset, use_old_formula) for gs_name, geneset in gmt.items()}
         ).T

def calculate_signatures(df_coh, df_cms_list):
    """
    This function counts signatures for all patients for all genes under study.

    :param df_coh: Cohort loaded in the read_data function
    :param df_cms_list: Examined genes loaded in the read_data function
    :return: Signatures for all patients for all genes investigated. Rows - patients
    """
    # Удаление столбца 'Unnamed: 0'
    expressions = df_coh.drop('Unnamed: 0', axis=1)

    # Вызов функции ssgsea_formula() с транспонированным датафреймом и списком генов
    signatures = ssgsea_formula(expressions.T, {'CMS4': df_cms_list}).T

    return signatures

def calculate_correlations(df_kin, df_coh, signatures):
    """
    This function takes all tested kinases, finds their expression in a cohort, and then correlates them with patient signatures
    for the tested genes.

    :param df_kin: Kinases loaded in analyze_data functions
    :param df_coh: Cohort loaded in the read_data function
    :param signatures: Signatures for all patients obtained in the calculate_signatures function
    :return: A dataset with gene names in the "Gene" column and their correlation to the signature in the "Correlation" column
    """
    correlations = []

    for gene_name in df_kin['Name']:
        if gene_name in df_coh.columns:
            # Извлечение столбца с геном из df_coh
            gene_expression = df_coh[gene_name]

            # Вычисление корреляции с столбцом CMS4 из signatures
            correlation = np.corrcoef(gene_expression, signatures['CMS4'])[0, 1]

            # Добавление корреляции в список
            correlations.append(correlation)
        else:
            # Если ген не найден в df_coh, добавляем NaN в список
            correlations.append(np.nan)

    # Создание DataFrame с результатами
    results = pd.DataFrame({'Gene': df_kin['Name'], 'Correlation': correlations})

    return results

def filter_and_print_correlations(results):
    """
    This function shows which kinases have the most obvious correlation (>0.85 or <-0.85), as well as the kinases with the lowest
    correlation (lower to -1 on the segment -1 - 1).

    :param results: Correlations of kinases with signatures obtained in the calculate_correlations function
    :return: Two datasets with values filtered by condition
    """
    # Фильтрация DataFrame по условию
    filtered_results = results[(results['Correlation'] > 0.85) | (results['Correlation'] < -0.85)]

    # Вывод имен и значений, удовлетворяющих условию
    print("Гены с корреляцией > 0.85 или < -0.85:")
    display(filtered_results[['Gene', 'Correlation']])

    # Находим индексы 5 строк с самыми маленькими значениями в столбце 'Correlation'
    smallest_correlations = results.nsmallest(5, 'Correlation')

    # Выводим названия генов и их значения корреляции
    print("\nГены с наименьшей корреляцией:")
    display(smallest_correlations[['Gene', 'Correlation']])

def normalize_correlation(results):
    """
    This function normalizes correlations and converts their values from -1 - 1 to 0 - 1 (add 1 to the value and then divide by 2).

    :param results: Correlations of kinases with signatures obtained in the calculate_correlations function
    :return: Dataset with gene names in the "Gene" column, their signature correlation in the "Correlation" column, and their
    normalized correlation in the "Normalized_Correlation" column
    """
    # Нормализуем корреляции и добавляем в новый столбец
    results['Normalized_Correlation'] = (results['Correlation'] + 1) / 2

    return results
