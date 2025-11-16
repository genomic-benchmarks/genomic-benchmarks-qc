import numpy as np
import pandas as pd
import logging
from genbenchQC.utils.bias_model import model


def flag_significant_differences(stats1, stats2):

    results = model(stats1, stats2)
    results['Unique bases'] = {}
    results['Unique bases']['Flag'] = flag_unique_bases(
        stats1, stats2
    )

    results['Duplicate sequences'] = {}
    results['Duplicate sequences']['Flag'] = flag_duplicate_sequences(
        stats1, stats2
    )

    results['Duplication between labels'] = {}
    results['Duplication between labels']['Flag'] = flag_duplication_between_datasets(
        stats1.sequences, stats2.sequences
    )

    results = pd.DataFrame.from_dict(results, orient='index')
    results.index.name = 'Statistic'

    return results

def flag_unique_bases(stats1, stats2):
    if set(stats1.stats['Unique bases']) == set(stats2.stats['Unique bases']):
        return 'Pass'
    else:
        return 'Fail'
    
def flag_duplicate_sequences(stats1, stats2):
    if stats1.stats['Number of sequences'] != stats1.stats['Number of sequences left after deduplication']:
        return 'Warning'
    if stats2.stats['Number of sequences'] != stats2.stats['Number of sequences left after deduplication']:
        return 'Warning'
    return 'Pass'

def flag_duplication_between_datasets(sequences1, sequences2):
    return "Fail" if bool(set(sequences1) & set(sequences2)) else "Pass"

# def flag_per_sequence_content(stats1, stats2, column, threshold):
    
#     df1 = stats1[column]
#     df2 = stats2[column]
    
#     # get columns names
#     bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

#     distances = {}
#     for base in bases:
#         if base not in df1 or base not in df2:
#             distances[base] = np.inf
#         else:
#             distances[base] = wasserstein_distance(df1[base], df2[base])
#             logging.debug(f"Distance for {base}: {distances[base]} (threshold: {threshold})")
#             max_value = max(
#                 max(df1[base]),
#                 max(df2[base])
#             )
#             if max_value > 0:
#                 distances[base] /= max_value
#                 logging.debug(f"Max value for {base}: {max_value}")
#                 logging.debug(f"Distance after normalization for {base}: {distances[base]}")

#     passed = np.all(np.array(list(distances.values())) < threshold)
    
#     return (distances, passed)

# def flag_per_position_nucleotide_content(stats1, stats2, column, threshold, end_position):
    
#     df1 = stats1[column]
#     df2 = stats2[column]

#     # get columns names
#     bases = list(set(list(df1.columns.values) + list(df2.columns.values)))

#     p_values = {}
#     passed = True
#     for base in bases:

#         p_values[base] = []
#         for i in range(end_position):
#             if base not in df1 or base not in df2:
#                 p_values[base].append(np.inf)
#             else:
#                 df1_base = df1[base][i]
#                 df2_base = df2[base][i]
#                 table=[[df1_base * 100, (1 - df1_base) * 100],
#                     [df2_base * 100, (1 - df2_base) * 100]]

#                 _, p_value = fisher_exact(table=table) 
#                 p_values[base].append(p_value)

#         # Correcting for FDR per base
#         _, p_values[base] = fdrcorrection(p_values[base])

#         passed = passed and np.all(np.array(p_values[base]) > threshold)
 
#     return (p_values, passed)

    
# def flag_per_sequence_one_stat(stats1, stats2, column, threshold):

#     distance = wasserstein_distance(
#         stats1[column].values.flatten(),
#         stats2[column].values.flatten()
#     )
#     logging.debug(f"Distance for {column}: {distance} (threshold: {threshold})")
#     max_value = max(
#         max(stats1[column].values.flatten()), 
#         max(stats2[column].values.flatten())
#     )
#     if max_value > 0:
#         distance /= max_value
#         logging.debug(f"Max value for {column}: {max_value}")
#         logging.debug(f"Distance after normalization: {distance}")

#     passed = distance < threshold

#     return (distance, passed)
