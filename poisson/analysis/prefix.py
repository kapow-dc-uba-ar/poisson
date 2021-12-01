# encoding: utf-8
from math import ceil

import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd

from poisson.analysis.poisson import fill_occurrences, \
    get_non_aligned_words_occurrences, get_aligned_words_occurrences


def non_aligned_set(sequence, lam, k, j):
    '''
    Returns a list of words of length k that appear exactly j times
    inside x[:floor(lam*(2**k))]
    
    '''
    x = sequence.sequence
    alphabet = ''.join([str(n) for n in range(sequence.base)])

    return [word for (word, count) in fill_occurrences(get_non_aligned_words_occurrences(x, lam, k), alphabet, k).items() if count == j]


def aligned_set(sequence, lam, k, j):
    '''
    Returns a list of words of length k that appear exactly j times
    inside x[:k*floor(lam*(2**k))] aligned to k
    
    
    :param sequence: a Sequence object
    :param lam:
    :param k:
    :param j:
    '''
    x = sequence.sequence
    alphabet = ''.join([str(n) for n in range(sequence.base)])
    return [word for (word, count) in fill_occurrences(get_aligned_words_occurrences(x, lam, k), alphabet, k).items() if count == j]


def prefix_distribution(occurrence_dict, prefix_length):
    prefix_occurrences = dict()
    for word in occurrence_dict:
        if len(word) >= prefix_length:
            prefix = word[:prefix_length]
            if prefix in prefix_occurrences:
                prefix_occurrences[prefix] += occurrence_dict[word]
            else:
                prefix_occurrences[prefix] = occurrence_dict[word]
    return prefix_occurrences


def prefix_distribution_from_list(words_list, prefix_length):
    prefix_occurrences = dict()
    for word in words_list:
        if len(word) >= prefix_length:
            prefix = word[:prefix_length]
            if prefix in prefix_occurrences:
                prefix_occurrences[prefix] += 1
            else:
                prefix_occurrences[prefix] = 1
    return prefix_occurrences


def get_multi_figure(word_occurrences, lam, k, set_name, seq_name):

    prefix_length = 6
    n_rows = 4
    n_cols = 2
    n_graphs = n_rows * n_cols
    
    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left")
    index = 1
    for prefix_length in range(1, 1 + n_graphs):
        prefix_dict = prefix_distribution(word_occurrences, prefix_length)
        df = pd.DataFrame.from_dict(prefix_dict, orient='index').sort_index()
        if len(prefix_dict) > 0:
            fila = ceil(index / n_cols)
            col = ((index - 1) % n_cols) + 1
            fig.add_trace(go.Bar(x=df.index, y=df[0], name=f'len(prefix): {prefix_length}'),
                      row=fila, col=col)
            fig.update_xaxes(type='category')
        index += 1
    fig.update_layout(title_text=f'Repetitions in {set_name} (lambda: {lam}, k: {k}, x: {seq_name}) by prefix')
    
    return fig


def display_prefix_hist(data_dict):
    df = pd.DataFrame.from_dict(data_dict, orient='index').sort_index()
    fig = px.bar(df, x=df.index, y=df[0], title='Prefix distribution')
    fig.update_xaxes(type='category')
    fig.show()


def plot_repetitions_lambda_k(sequence, lambdas, ks, prefix_length, set_name, wo_function=get_non_aligned_words_occurrences):

    # prefix_length = 6
    n_rows = len(lambdas)
    n_cols = len(ks)

    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left")
    index = 1
    alphabet = ''.join([str(n) for n in range(sequence.base)])
    for lam in lambdas:
        for k in ks:
            word_occurrences = fill_occurrences(wo_function(sequence.sequence, lam, k), alphabet, k)

            prefix_dict = prefix_distribution(word_occurrences, prefix_length)
            df = pd.DataFrame.from_dict(prefix_dict, orient='index').sort_index()
            fila = ceil(index / n_cols)
            col = ((index - 1) % n_cols) + 1
            fig.add_trace(go.Bar(x=df.index, y=df[0], name=f'lambda: {lam}; k: {k}'),
                    row=fila, col=col)
            fig.update_xaxes(type='category')
            index += 1
    fig.update_layout(title_text=f'Repetitions in {set_name} with prefix length {prefix_length}')

    return fig


def plot_prefix_repetitions_k_j(sequence, ks, lams, js, prefix_length, compute_func=aligned_set):
    n_rows = len(js)
    n_cols = len(ks)

    subplot_titles = [f'j: {j}; k: {k}' for j in js for k in ks ]
    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left", subplot_titles=subplot_titles)

    index = 1
    for j in js:
        for k_idx, k in enumerate(ks):
            # word_occurrences = get_non_aligned_words_occurrences(xs, lam, k)
            lam = lams[k_idx] if k_idx < len(lams) else lams[0]
            wo = compute_func(sequence, lam, k, j)
            if len(wo) > 0: 
                prefix_dict = prefix_distribution_from_list(wo, prefix_length)
                df = pd.DataFrame.from_dict(prefix_dict, orient='index').sort_index()
                fila = ceil(index / n_cols)
                col = ((index - 1) % n_cols) + 1
                fig.add_trace(go.Bar(x=df.index, y=df[0], name=f'j: {j}; k: {k}; lambda: {lam}'),
                        row=fila, col=col)
                fig.update_xaxes(type='category')
            index += 1
    fig.update_layout(title_text=f'Repetitions in {compute_func.__name__} for {sequence.name} with prefix length {prefix_length}')

    return fig


def plot_prefix_repetitions_xs_js_lam(sequence_list, js, k, lam, prefix_length, set_function=aligned_set):
    n_rows = len(sequence_list)
    n_cols = len(js)

    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left")
    index = 1
    for sequence in sequence_list:
        for j in js:
            # word_occurrences = get_non_aligned_words_occurrences(xs, lam, k)
            wo = set_function(sequence, lam, k, j)

            if len(wo) > 0:
                prefix_dict = prefix_distribution_from_list(wo, prefix_length)
                df = pd.DataFrame.from_dict(prefix_dict, orient='index').sort_index()
                fila = ceil(index / n_cols)
                col = ((index - 1) % n_cols) + 1
                fig.add_trace(go.Bar(x=df.index, y=df[0], name=f'x: {sequence.name}; j: {j}'),
                        row=fila, col=col)
                fig.update_xaxes(type='category')
            index += 1
    fig.update_layout(title_text=f'Repetitions in {set_function.__name__} with prefix length {prefix_length}; lam: {lam}; k: {k}')

    return fig


def plot_prefix_repetitions_xs_js(sequence_list, js, k, prefix_length, set_function=aligned_set):
    n_rows = len(sequence_list)
    n_cols = len(js)

    subplot_titles = [f'{sequence.name}; j:{j}' for sequence in sequence_list for j in js ]
    fig = make_subplots(rows=n_rows, cols=n_cols, start_cell="top-left", subplot_titles=subplot_titles)
    index = 1
    for sequence in sequence_list:
        alphabet = ''.join([str(n) for n in range(sequence.base)])
        for j in js:
            # word_occurrences = get_non_aligned_words_occurrences(xs, lam, k)
            lam = 1 / len(alphabet)
            wo = set_function(sequence, lam, k, j)

            if len(wo) > 0:
                prefix_dict = prefix_distribution_from_list(wo, prefix_length)
                df = pd.DataFrame.from_dict(prefix_dict, orient='index').sort_index()
                fila = ceil(index / n_cols)
                col = ((index - 1) % n_cols) + 1
                fig.add_trace(go.Bar(x=df.index, y=df[0], name=f'x: {sequence.name}; j: {j}; lam: {lam}'),
                        row=fila, col=col)
                fig.update_xaxes(type='category')
            index += 1
    fig.update_layout(title_text=f'Repetitions in {set_function.__name__} with prefix length {prefix_length}; k: {k}')

    return fig

