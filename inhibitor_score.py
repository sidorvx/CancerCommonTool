import pandas as pd
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
from typing import Dict, List, Literal
from IPython.display import display

def preprocess_affinity_data(df_aff):
    """
    This function preprocesses the read dataset with drug affinities. It replaces all affinity values with inverse values
    (this step is performed based on the formula for calculating the efficacy metric) and replaces all NaNs with 0.

    :param df_aff: Affinities of preparations loaded in read_data functions
    :return: Conditionally processed dataset with drug affinities
    """
    # Заменяем все NaN на 0
    df_aff.replace(np.nan, 0, inplace=True)

    # Заменяем числовые значения на их обратные значения
    df_aff = df_aff.applymap(lambda x: 1/x if isinstance(x, (int, float)) and x != 0 else x)

    return df_aff

def calculate_drug_scores(df_aff_red, results):
    """
    This function calculates the efficacy metric for each drug using the following methodology. For each drug target, the correlation
    value is found and then multiplied by the inverse affinity of that drug to that target. The results for all targets are then
    summarized within a single drug.

    :param df_aff_red: Edited dataset with drug affinities obtained in the preprocess_affinity_data function
    :param results: Normalized correlations of kinases with signatures obtained in the normalize_correlation function
    :return: A dataset in which the Drugs column contains the drugs and the Score column contains the efficacy metrics for the drugs
    """
    genes1 = set(results['Gene'])
    genes2 = set(df_aff_red.columns)
    genes = list(genes1.intersection(genes2))
    df_aff_red = df_aff_red.set_index('Unnamed: 0')
    df_aff_red.index.name = 'Drugs'  # Переименуем индекс в 'Drugs'
    results = results.set_index('Gene')
    relevant_results = results.loc[genes]
    ordered_aff = df_aff_red.loc[:,genes]

    # Извлекаем имена препаратов из df_aff_red
    drug_names = ordered_aff.index

    # Умножаем значения столбцов df_aff_red на значения из relevant_results и суммируем по строкам
    scores = ordered_aff.dot(relevant_results['Normalized_Correlation'])

    # Создаем DataFrame с результатами
    scores_df = pd.DataFrame({'Score': scores})

    return scores_df  # Возвращаем только столбец с оценками

def print_top_drugs(scores, top_n=5):
    """
    This function outputs the top of the best drugs (with the highest metric).

    :param scores: Performance metric values obtained in the calculate_drug_scores function
    :param top_n=n: Number is the number of drugs with the best metric to be derived
    :return: Dataset with the conditionally filtered value of the drug efficacy metric
    """
    # Сортировка DataFrame по столбцу "Score" в порядке убывания
    sorted_scores = scores.sort_values(by='Score', ascending=False)

    # Выводим указанное количество лучших препаратов
    print(f"{top_n} лучших препаратов:")
    display(sorted_scores.head(top_n))