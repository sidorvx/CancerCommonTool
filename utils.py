import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from typing import Dict, List, Literal
from IPython.display import display

def read_data(df_coh_path, df_cms_list):
    """
    Reads the data you feed into the input: expression cohort and genes under study

    :param df_coh_path: The path to the .tsv file, where the rows are patients and the columns are genes. Expressions in log2(TMP+1)
    format
    :param df_cms_list: Python list with the genes you're researching
    :return: Data read out
    """
    data = {}

    data['df_coh'] = pd.read_csv(df_coh_path, sep='\t')
    data['df_cms_list'] = df_cms_list

    return data

def analyze_data(data, df_kin_path, df_aff_path, df_cms_list):
    """
    Subloads human kinome and drug affinity data, and analyzes whether all of the genes and kinases under investigation are in
    your cohort

    :param data: Your downloaded data from the previous function (cohort and genes)
    :param df_kin_path: The path to the human kinome file. Names of kinases in the "Name" column
    :param df_aff_path: Path to the file with drugs with their affinities to targets. Rows - drugs, columns - targets
    :param df_cms_list: Python list with the genes you're researching
    :return: All data required for the work (cohort, genes, kinome, preparations) and a warning if kinases or genes under
    study are missing from the cohort
    """
    data['df_kin'] = pd.read_csv(df_kin_path)
    data['df_aff'] = pd.read_excel(df_aff_path)

    df_coh_genes = set(data['df_coh'].columns)
    df_cms_genes = df_cms_list

    missing_kinases = [kinase for kinase in set(data['df_kin']['Name']) if kinase not in df_coh_genes]
    missing_genes = [gene for gene in df_cms_genes if gene not in df_coh_genes]

    if missing_kinases:
        print("Предупреждение: Следующие киназы из df_kin отсутствуют в df_coh:")
        display(pd.DataFrame(missing_kinases, columns=['Missing Kinases']))

    if missing_genes:
        print("Предупреждение! Следующие гены из df_cms_list отсутствуют в df_coh:")
        display(pd.DataFrame(missing_genes, columns=['Missing Genes']))

    return data

def plot_score_distribution(scores):
    """
    This function constructs the distribution of efficiency metrics.

    :param scores: Performance metric values obtained in the calculate_drug_scores function
    :return: A graph of the distribution of the performance metric data. The X axis is the value of the metric, and the Y axis
    is the frequency of occurrence
    """
    # Получение значений score
    score_values = list(scores['Score'])  # Получаем значения оценок из столбца 'Score'

    # Построение гистограммы
    plt.hist(score_values, bins=10, color='skyblue', edgecolor='black')

    # Настройка заголовка и меток осей
    plt.title('Распределение данных метрик эффективности')
    plt.xlabel('Значение метрики')
    plt.ylabel('Частота встречаемости')

    # Отображение гистограммы
    plt.show()
