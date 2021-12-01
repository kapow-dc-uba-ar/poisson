# encoding: utf-8
from itertools import product
from math import floor, ceil

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

def non_aligned_count(a_sequence, k, lambda_value=1):
    '''
    
    :param a_sequence: a poisson.Sequence object
    :param k:
    :param lambda_value:
    '''
    return non_aligned_count_base(a_sequence.sequence, k, a_sequence.base, lambda_value)


def non_aligned_count_base(x, k, base=10, lam=1):
    '''
    Con superposicion (overlapping)
    Non aligned
    :param x: a sequence
    :param k:
    :param base: a sequence base
    :param lam: lambda value
    '''
    r = 0
    window = 0
    mx = 0
    cnt = [0] * (base ** k)
    lookup_limit = floor(lam * (base ** k)) + (k - 1)
    for l in range(lookup_limit):
        while r - l < k:
            window = base * window + (ord(x[r]) - ord('0'))
            r += 1
        
        window = window % base ** k
        
        cnt[window] += 1
        mx = max(mx, cnt[window])
    return cnt, mx


def aligned_count(a_sequence, k, lambda_value=1):
    base = a_sequence.base
    x = a_sequence.sequence
    r = 0
    window = 0
    mx = 0
    cnt = [0] * (base ** k)
    # lookup_limit = floor(lam*(base**k)) + (k - 1)
    lookup_limit = k * floor(lambda_value * (base ** k)) + (k - 1)
    for l in range(0, lookup_limit, k):
        try:
            window = int(x[l:l + k], base)
        except:
            print('warning! seq length:', len(x), f' ; index: {l}:{l+k}')
        cnt[window] += 1
        mx = max(mx, cnt[window])
    return cnt, mx


def get_frequencies(a_sequence, k=8, lam=1, count_function=non_aligned_count):
    '''
    
    :param a_sequence:
    :param k:
    :param lam:
    :param count_function:
    '''

    cnt, mx = count_function(a_sequence, k, lam)
    cnt_j = [0] * (mx + 1)
    for j in cnt:
        cnt_j[j] += 1
    freq = [val / (a_sequence.base ** k) for val in cnt_j]
    return freq


# generate words of length k from alphabet
def words_from(alphabet, k):
    return [''.join(x) for x in product(alphabet, repeat=k)]


def k_words_generator(sequence, k, increment=1):
    index = 0
    while index < len(sequence) - k:
        yield sequence[index: index + k]
        index += increment


def filter_dict(dictObj, callback):
    newDict = dict()
    # Iterate over all the items in dictionary
    for (key, value) in dictObj.items():
        # Check if item satisfies the given condition then add to new dict
        if callback((key, value)):
            newDict[key] = value
    return newDict


def get_words_occurrences(x, lam, k, inc=1):
    '''
    Returns the number of occurrences for each word of length k
    that's present on x
    @parameter inc: step to find the next word (default: 1)
    '''
    results = {}
    gen = k_words_generator(x, k, increment=inc)
    
    for word in gen:
        if word in results:
            results[word] += 1
        else:
            results[word] = 1
    
    return results


def get_non_aligned_words_occurrences(x, lam, k):
    '''
    Returns word occurrences according to R set criteria
    '''
    lookup_limit = floor(lam * (2 ** k)) + (k - 1)
    if lookup_limit > len(x):
        print(f'WARNING: sequence x is too short (<{lookup_limit})')
    
    return get_words_occurrences(x[:lookup_limit], lam, k, 1)


def get_aligned_words_occurrences(x, lam, k):
    '''
    Returns word occurrences according to Q set criteria
    '''
    lookup_limit = k * floor(lam * (2 ** k))
    if lookup_limit > len(x):
        print(f'WARNING: sequence x is too short (<{lookup_limit})')
    
    return get_words_occurrences(x[:lookup_limit], lam, k, k)


def fill_occurrences(word_occurrences, alphabet, word_length):
    '''
    Returns the dictionary word_occurrences filled with non appearing
    words of length word_length from alphabet
    '''
    all_words_dict = {k:0 for k in words_from(alphabet, word_length)}
    all_words_dict.update(word_occurrences)
    return all_words_dict
  

def plot_j_distribution(sequence_list, ks, lam, wo_function=get_non_aligned_words_occurrences):
    '''
    For each sequence and k combination plots
    the distribution of the J value (number of words repeating
    exactly j times in the sequence).
    :param sequence_list:
    :param ks: a list of k values
    :param lam: a lambda value
    '''
    n_rows = len(sequence_list)
    n_cols = len(ks)

    subplot_titles = [f'{seq.name}; k: {k}' for seq in sequence_list for k in ks ]
    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left", subplot_titles=subplot_titles)

    index = 1
    for sequence in sequence_list:
        alphabet = ''.join([str(n) for n in range(sequence.base)])
        for k in ks:
            wo = fill_occurrences(wo_function(sequence.sequence, lam, k), alphabet, k)  # muy lento

            if len(wo) > 0:
                df = pd.DataFrame.from_dict(wo, orient='index').sort_index()
                fila = ceil(index / n_cols)
                col = ((index - 1) % n_cols) + 1
                fig.add_trace(go.Histogram(x=df[0], name=f'x: {sequence.name}; k: {k}'),
                        row=fila, col=col)
                fig.update_xaxes(type='category')
            index += 1
    fig.update_layout(title_text=f'J values distribution (amount of words appearing exactly j times)  lam:{lam}')

    return fig

