from utils import read_data, analyze_data, plot_score_distribution
from correlation_calc import ssgsea_score, ssgsea_formula, calculate_signatures, calculate_correlations, filter_and_print_correlations, normalize_correlation
from inhibitor_score import preprocess_affinity_data, calculate_drug_scores, print_top_drugs

def master_function(df_kin_path, df_coh_path, df_aff_path, df_cms_list, top_n=5):
    """
    This function is a master function and is a combination of all the functions presented above. It outputs all key results obtained
    in all steps.

    :param df_kin_path: The path to the human kinome file. Names of kinases in the "Name" column
    :param df_coh_path: The path to the .tsv file, where the rows are patients and the columns are genes. Expressions in log2(TMP+1)
    format
    :param df_aff_path: Path to the file with drugs with their affinities to targets. Rows - drugs, columns - targets
    :param df_csm_list: Python list with the genes you're researching
    :param top_n=n: Number is the number of drugs with the best metric to be derived
    :return: Kinases missing in the cohort, genes needed for the study, genes missing in the cohort; 2 filtered dataframes.
    The first contains genes with correlation > 0.85 or < -0.85, and the second contains genes with the lowest correlation; Dataframe
    with 2 columns - Drugs and Score. They contain a pair of drug name - efficacy metric respectively; A graph on the abscissa axis
    with a range of efficacy metric values and on the ordinate axis - frequency of occurrence of these values; A dataframe with
    2 columns and the number of rows you specify - Drugs and Score. They contain a pair of drug name - efficacy metric respectively.
    Only drugs are arranged in descending order of effectiveness metric.
    """
    # Чтение данных пользователя
    data = read_data(df_coh_path, df_cms_list)
    df_coh = data['df_coh']
    df_cms_list = data['df_cms_list']

    # Подгрузка и анализ всех данных
    data = analyze_data(data, df_kin_path, df_aff_path, df_cms_list)
    df_kin = data['df_kin']
    df_aff = data['df_aff']

    # Вычисление сигнатур
    signatures = calculate_signatures(df_coh, df_cms_list)

    # Расчет корреляций
    results = calculate_correlations(df_kin, df_coh, signatures)

    # Фильтрация и вывод корреляций
    filter_and_print_correlations(results)

    # Нормализация корреляций
    normalized_results = normalize_correlation(results)

    # Предобработка данных о взаимодействии
    df_aff_red = preprocess_affinity_data(df_aff)

    # Расчет оценок препаратов
    scores = calculate_drug_scores(df_aff_red, normalized_results)
    print(scores)

    # Построение распределения оценок
    plot_score_distribution(scores)

    # Вывод топ-N препаратов
    print_top_drugs(scores, top_n=top_n)

# Пример использования мастер-функции
df_kin_path = r'C:\Users\sidor\Desktop\Daniil\BG\Algoritm\HumanKinome.csv'
df_coh_path = r'C:\Users\sidor\Desktop\Daniil\BG\Algoritm\transposed_expressions.tsv'
df_aff_path = r'C:\Users\sidor\Desktop\Daniil\BG\Algoritm\Affinity.xlsx'
df_cms_list = ['ADAM8', 'AGR2', 'ANXA1', 'ANXA3', 'ANXA8L2', 'AREG', 'B3GNT3', 'BIK', 'C15orf48', 'C19orf33', 'CDH1', 'CDH3', 'CEACAM6', 'COL17A1', 'CST6', 'CTSL2', 'DEFB1', 'DHRS9', 'DKK1', 'EPCAM', 'EPS8L1', 'FA2H', 'FGFBP1', 'FXYD3', 'GALNT3', 'GPR87', 'HS3ST1', 'IFI27', 'IL20RB', 'ITGB4', 'KCNK1', 'KIAA1522', 'KLK10', 'KLK8', 'KRT15', 'KRT16', 'KRT17', 'KRT17P3', 'KRT18P30', 'KRT19', 'KRT6A', 'KRT6B', 'KRT6C', 'KRT7', 'LAD1', 'LAMB3', 'LAMC2', 'LCN2', 'LEMD1', 'LOC440335', 'LY6D', 'MAL2', 'MIA', 'MMP7', 'NMU', 'OCIAD2', 'PERP', 'PHACTR3', 'PKP3', 'PLAT', 'PLAU', 'PLEK2', 'PMEPA1', 'PRSS22', 'PSCA', 'S100A14', 'S100A2', 'S100A4', 'S100A6', 'S100P', 'SCEL', 'SERPINB3', 'SERPINB4', 'SLC2A1', 'SLPI', 'SPINT2', 'SPRR1B', 'SPRR3', 'STYK1', 'TACSTD2', 'TFF2', 'TMPRSS4', 'TNS4', 'TSPAN1', 'UCA1', 'WFDC2']

master_function(df_kin_path, df_coh_path, df_aff_path, df_cms_list, top_n=5)
