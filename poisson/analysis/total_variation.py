# encoding: utf-8
import plotly.express as px
import plotly.graph_objects as go

import pandas as pd
from scipy.stats import poisson
from poisson.analysis.poisson import non_aligned_count, aligned_count
import math


# Total variation

def total_variation(freq, lam=1):
    remaining = 1
    distance = 0
    for j in range(len(freq)):
        p = poisson.pmf(j,lam)
        distance += abs(freq[j] - p)
        remaining -= p
    distance += remaining
    return distance / 2


def get_variation_limit(a_sequence, max_k=8, lam=1, count_function=non_aligned_count):
    variations = []
    for k in range(1, max_k+1):
        cnt, mx = count_function(a_sequence, k, lam)
        cnt_j = [0] * (mx+1)
        for j in cnt:
            cnt_j[j] += 1
        freq = [val /  (a_sequence.base**k) for val in cnt_j]
        variations.append(total_variation(freq, lam))
    return variations


def plot_variation_limit(a_sequence, max_k=8, lam=1, count_function=non_aligned_count):
    variations = get_variation_limit(a_sequence, max_k, lam, count_function)
    fig = px.line(x=range(len(variations)),y=variations, log_y=True)
    fig.update_layout(
        title=f'Total variation for {a_sequence.name}, lambda: {lam}',
        xaxis_title="K",
        yaxis_title="Total variation",
    )

    fig.show()


def total_variation_dataframe(sequence_list, lambda_value=1):
    df = pd.DataFrame()
    
    for a_sequence in sequence_list:
        max_k = math.floor(math.log(a_sequence.length / lambda_value, a_sequence.base)) - 1
        seq_var = get_variation_limit(a_sequence, max_k=max_k, lam=lambda_value, count_function=non_aligned_count)
        df = pd.concat([df, pd.DataFrame({a_sequence.name: seq_var})], axis=1)
    return df

def total_variation_aligned_dataframe(sequence_list, lambda_value=1):
    df = pd.DataFrame()
    
    for a_sequence in sequence_list:
        # lookup_limit = k*floor(lambda_value*(base**k)) + (k - 1)
        scale_aprox = 5 # ~= log_2(30)
        max_k = math.floor(math.log(a_sequence.length / lambda_value, a_sequence.base)  / scale_aprox)
        seq_var = get_variation_limit(a_sequence, max_k=max_k, lam=lambda_value, count_function=aligned_count)
        df = pd.concat([df, pd.DataFrame({a_sequence.name: seq_var})], axis=1)
    return df


def plot_total_variation_comparison(df, title='Sequences and total variation (Poisson) comparison'):
    fig = go.Figure()

    for col_name in df:
        fig.add_trace(go.Scatter(
            x=df.index,
            y=df[col_name],
            name=col_name,
        ))

    fig.update_yaxes(type="log")
    fig.update_layout(
        title=title,
        xaxis_title="K",
        yaxis_title="total variation (log scale)",
    )
    fig.show()
    
    
def get_max_j_values(sequence, max_k, lam, count_function=non_aligned_count):
    '''
    Returns the min J that doesn't have elements for every K between 0 and max_k
    '''
    values = []
    for k in range(1, max_k):
        cnt, mx = count_function(sequence, k, lambda_value=lam)
        values.append(mx)
    return values


def plot_max_j(sequence, max_k, lambda_set, count_function=non_aligned_count):
    '''
    Plots the min J that doesn't have elements for every K between 0 and max_k
    ie.: plot_max_j(r2_seq, 2, 23, [1/5, 1/3, 1/2, 1, 2, 3, 5], 'Random2')
    
    :param sequence: A Sequence object
    :param base: A sequence base (2, 5, 10, etc)
    :param max_k:
    :param lambda_set: a list or single value of lambda
    :param seq_name: Sequence name
    '''
    
    if not isinstance(lambda_set, list):
        lambda_set = [lambda_set]
    
    fig = go.Figure()
    for lam in lambda_set:
        values = get_max_j_values(sequence, max_k, lam, count_function)
        fig.add_trace(go.Scatter(
            x=list(range(1, max_k)),
            y=values,
            name=f'Lambda: {lam}'

        ))
    fig.update_layout(
        title=f'Min J values with no occurrences. Sequence: {sequence.name}',
        xaxis_title="K",
        yaxis_title="J",
    )
    fig.show()
