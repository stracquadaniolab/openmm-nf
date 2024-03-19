#!/usr/bin/env python3
"""output_data
Calculates ΔΔG relative to the wildtype for using minimum, lowest 3, 5, 10 and 20 outputted folding energies then measures Spearman's Ranking relative to experimental/theoretical benchmarks
Usage:
output_data.py [--bench=<csv>] [--test=<csv>]
Options:
--bench=<csv>    Read original experimental/theoretical data from csv file
--test=<csv>     Read calculated data from openmm-minimise process 
"""
import logging
import csv
import pandas as pd
import numpy as np
from docopt import docopt
import sys
from sys import stdout

def get_energy(test_csv: str):
    df2 = pd.read_csv(test_csv)
    df3 = df2.loc[df2['name'] == 'wildtype']
    wt_min = df3.iloc[0]['min']
    wt_low3 = df3.iloc[0]['lowest3']
    wt_low5 = df3.iloc[0]['lowest5']
    wt_low10 = df3.iloc[0]['lowest10']
    wt_low20 = df3.iloc[0]['lowest20']
    df2['min_ΔΔG'] = df2['min'] - wt_min
    df2['lowest3_ΔΔG'] = df2['lowest3'] - wt_low3
    df2['lowest5_ΔΔG'] = df2['lowest5'] - wt_low5
    df2['lowest10_ΔΔG'] = df2['lowest10'] - wt_low10
    df2['lowest20_ΔΔG'] = df2['lowest20'] - wt_low20
    df2.drop(df2.index[(df2['name'] == 'wildtype')],axis=0,inplace=True)
    df2 = df2[['name', 'min_ΔΔG', 'lowest3_ΔΔG', 'lowest5_ΔΔG', 'lowest10_ΔΔG', 'lowest20_ΔΔG']]
    return df2

def spearman_rank(bench_csv: str, df2):
    df1 = pd.read_csv(bench_csv)
    mergeDf = pd.merge(df1, df2, left_on = ['name'], right_on = ['name'])
    min_rank = mergeDf['ΔΔG'].corr(mergeDf['min_ΔΔG'], method='spearman')
    low3_rank = mergeDf['ΔΔG'].corr(mergeDf['lowest3_ΔΔG'], method='spearman')
    low5_rank = mergeDf['ΔΔG'].corr(mergeDf['lowest5_ΔΔG'], method='spearman')
    low10_rank = mergeDf['ΔΔG'].corr(mergeDf['lowest10_ΔΔG'], method='spearman')
    low20_rank = mergeDf['ΔΔG'].corr(mergeDf['lowest20_ΔΔG'], method='spearman')
    logging.info("Min ΔΔG Spearman's Rank  = %.9f "
            % min_rank)
    logging.info("Lowest 3 ΔΔG Spearman's Rank  = %.9f "
            % low3_rank)
    logging.info("Lowest 5 ΔΔG Spearman's Rank  = %.9f "
            % low5_rank)
    logging.info("Lowest 10 ΔΔG Spearman's Rank  = %.9f "
            % low10_rank)
    logging.info("Lowest 20 ΔΔG Spearman's Rank  = %.9f "
            % low20_rank)
    all_ranks = {'min_rank': [min_rank], 'low3_rank': [low3_rank], 'low5_rank': [low5_rank], 'low10_rank': [low10_rank], 'low20_rank': [low20_rank]}
    df_sr = pd.DataFrame(all_ranks)
    return df_sr

def main():
    arguments = docopt(__doc__, version='output_data.py')
    df2 = get_energy(arguments['--test'])
    df2.to_csv('data_ΔΔG.csv', mode='w', indexr=False)
    df_sr = spearman_rank(arguments['--bench'], df2)
    df_sr.to_csv('data_ΔΔG-spearman.csv', mode='w', indexr=False)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmm-minimise.py')
    logging.getLogger().setLevel(logging.INFO)
    main()
