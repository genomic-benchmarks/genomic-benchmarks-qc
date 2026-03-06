import pandas as pd
from genbenchQC.utils.bias_model import model


def flag_significant_differences(stats1, stats2, max_class_size=None):

    results = {}

    results['Unique bases'] = {}
    results['Unique bases']['Flag'] = flag_unique_bases(
        stats1, stats2
    )

    results['Sequence Duplications within Labels'] = {}
    results['Sequence Duplications within Labels']['Flag'] = flag_duplicate_sequences(
        stats1, stats2
    )

    results['Duplicate Sequences between Labels'] = {}
    results['Duplicate Sequences between Labels']['Flag'] = flag_duplication_between_datasets(
        stats1.sequences, stats2.sequences
    )

    model_results = model(stats1, stats2, max_class_size=max_class_size)

    results.update(model_results)

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
