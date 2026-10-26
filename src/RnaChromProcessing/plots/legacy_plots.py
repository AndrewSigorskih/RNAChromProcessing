from pathlib import Path
from typing import List, Literal

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import seaborn as sns


FIGSIZE = (11.7, 8.27)
TEN_KB = 10_000

#----------------------------------------#
#         UTILITY FUNCTIONS              #
#----------------------------------------#

def _calc_upper_whisker(arr: pd.Series) -> float:
    Q1, Q3 = np.percentile(arr, [25, 75])
    IQR = Q3 - Q1
    hival = Q3 + 1.5 * IQR
    return arr[arr <= hival].max()


def _transform_raw_counts_df(raw_wins: pd.DataFrame,
                             merged_index: List[str],
                             strandness: Literal[0, 1]) -> pd.DataFrame:
    counts = raw_wins.applymap(lambda x: x[strandness])
    counts['id'] = merged_index
    counts = pd.melt(counts, id_vars='id')
    counts['value'] = counts['value'].apply(np.log10)
    return counts


def set_style_white() -> None:
    sns.set_style('white')
    sns.set_palette('husl')
    plt.rc(
        'font',  
        serif = 'Ubuntu',
        monospace = 'Ubuntu Mono',
        size = 10
    )
    plt.rc(
        'axes', 
        labelsize = 16,
        labelweight = 'bold',
        labelpad = 10,
        titlesize = 22,
        titlepad = 10,
        titleweight = 'bold'
    )
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14
    plt.rcParams['legend.fontsize'] = 14


#----------------------------------------#
#           DETECT-STRAND                #
#----------------------------------------#


def rna_strand_barplot(wins: pd.DataFrame,
                       total_genes: int,
                       out_dir: str,
                       prefix: str,
                       plot_format = 'png') -> None:
    # init figure
    fig, ax = plt.subplots()
    fig.set_size_inches(*FIGSIZE)
    ax.set_ylim(-total_genes-2, total_genes+2)
    # variables
    index = wins.index
    width = 0.5  # the width of the bars
    groups = index.get_level_values(0)
    labels = ['_'.join(x) for x in index.to_flat_index()]
    x = np.arange(len(labels))  # the label locations
    # groups coloring
    colors =[]
    patches=[]
    palette = sns.color_palette("husl", len(groups.unique()))
    for i, name in enumerate(groups.unique()):
        colors += [palette[i]] * groups[groups==name].shape[0]
        patches.append(
            mpatches.Patch(
                color=palette[i], 
                label=name
            )
        )
    # plot bars
    rects = []
    negrects = []
    for i in x:
        rects += ax.bar(i, wins.iat[i, 0], width, color=colors[i], alpha=.5)
        negrects += ax.bar(i, -wins.iat[i, 1], width, color=colors[i], alpha=.5)
    # labels and legend
    ax.set_ylim(-total_genes-3, total_genes+3)
    ax.set_ylabel('Numbers of wins', fontsize=20)
    ax.set_title(f'Numbers of wins and losses\nout of {total_genes} genes', fontsize=16)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=85)
    ax.legend(handles=patches, fontsize=14, loc='best')
    ax.set_yticks([])
    
    # exact height values over rects
    def autolabel(rects, offset=1):
        """Attach a text label above each bar in *rects*, displaying its height."""
        for rect in rects:
            height = round(rect.get_height(), 2)
            ax.annotate(f'{np.abs(height)}',
                        xy=(rect.get_x()+width/2, height),
                        xytext=(0, 3*offset),  # 3 points vertical offset
                        textcoords='offset points',
                        ha='center', va='bottom')
    autolabel(rects)
    autolabel(negrects, -3)
    ax.axhline(color='grey')
    # save
    if plot_format == 'png':
        plt.savefig(f'{out_dir}/{prefix}_wins.png', dpi=300, bbox_inches='tight')
    elif plot_format == 'svg':  
        plt.savefig(f'{out_dir}/{prefix}_wins.svg', dpi=300, bbox_inches='tight', format='svg')


def rna_strand_boxplot(raw_wins: pd.DataFrame,
                       total_genes: int,
                       out_dir: str,
                       prefix: str,
                       plot_format='png') -> None:
    # init figure
    fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True)
    fig.set_size_inches(11.7, 8.27)
    index = raw_wins.index
    groups = index.get_level_values(0)
    labels = ['_'.join(x) for x in index.to_flat_index()]
    x = np.arange(len(labels))  # the label locations
    # groups coloring
    colors =[]
    patches=[]
    palette = sns.color_palette("husl", len(groups.unique()))
    for i, name in enumerate(groups.unique()):
        colors += [palette[i]] * groups[groups==name].shape[0]
        patches.append(
            mpatches.Patch(
                color=palette[i], 
                label=name
            )
        )
    # split data to 2 dfs and log transform
    same_counts = _transform_raw_counts_df(raw_wins, labels, 0)
    anti_counts = _transform_raw_counts_df(raw_wins, labels, 1)

    # plot 2 boxplots
    sns.boxplot(same_counts, x='id', y='value', ax=ax0)
    sns.boxplot(anti_counts, x='id', y='value', ax=ax1)
    # customize axes
    ax0.set(xlabel=None)
    ax1.set(xlabel=None)
    ax1.set_xticklabels(labels, rotation=85)
    ax0.set_ylabel('Genes wins distribution', fontsize=12)
    ax1.set_ylabel('Genes losses distribution', fontsize=12)
    ax0.legend(handles=patches, fontsize=12, bbox_to_anchor=(1.01, 1))
    ax0.set_title(f'Coverage distributions of {total_genes} genes,\nlog10-transformed', fontsize=16)
    fig.tight_layout()
    # save
    if plot_format == 'png':
        plt.savefig(f'{out_dir}/{prefix}_cov_boxplot.png', dpi=300, bbox_inches='tight')
    elif plot_format == 'svg':
        plt.savefig(f'{out_dir}/{prefix}_cov_boxplot.svg', dpi=300, bbox_inches='tight', format='svg')


#----------------------------------------#
#                 XRNA                   #
#----------------------------------------#


def plot_length_distribution(tab: pd.DataFrame,
                             out_dir: Path,
                             prefix: str) -> None:
    fig, axes = plt.subplots()
    fig.set_size_inches(*FIGSIZE)
    length = tab['end'] - tab['start'] + 1
    length[length <= TEN_KB].hist(ax=axes, color='grey',
                                  grid=False, bins=20)
    title = 'X-RNAs lengths distribution'
    plt.title(title)
    axes.set_ylabel('Count', fontsize=20)
    axes.set_xlabel('Length, bp', fontsize=20)
    info = (
        f'min: {length.min()}\n'
        f'max: {length.max()}\n'
        f'median: {length.median():.1f}'
    )
    plt.text(0.75*axes.dataLim.bounds[2],
             0.75*axes.dataLim.bounds[3],
             info, ha='left', va='center',
             fontsize=20)
    # save
    plt.savefig(out_dir / f'{prefix}_length_distr.png',
                dpi=300, bbox_inches='tight')
    

def plot_distance_to_closest(tab: pd.DataFrame,
                             out_dir: Path,
                             prefix: str) -> None:
    fig, ax = plt.subplots()
    fig.set_size_inches(*FIGSIZE)
    palette = {'3\'': 'r', '5\'': 'b'}
    # select 5` and 3` distances that are within 10kb window around genes
    tab = tab[['closest_gene_side', 'closest_gene_dist']].copy()
    tab.columns = ['Prime', 'Distance']
    tab['Distance'] = tab['Distance'].abs()
    tab = tab[tab['Distance'] <= TEN_KB]
    # plot histograms
    ax2 = ax.twinx()
    sns.histplot(
        data=tab, x='Distance', hue='Prime', ax=ax,
        palette=palette
    )
    sns.kdeplot(
        data=tab, x='Distance', hue='Prime', ax=ax2,
        clip=(0.0, None), palette=palette, legend=False
    )
    # labels
    title = 'Distance to closest gene by strand'
    plt.title(title)
    ax.set_ylabel('Count', fontsize=20)
    ax.set_xlabel('Distance, bp', fontsize=20)
    ax2.set_ylabel('Density', fontsize=20)
    # save
    plt.savefig(out_dir / f'{prefix}_closest_gene_distances.png',
                dpi=300, bbox_inches='tight')
    

def plot_tpm_expressions(tab: pd.DataFrame,
                         out_dir: Path,
                         prefix: str) -> None:
    # data prep
    colnames = [x for x in tab.columns if x.endswith('_TPM')]
    tab = tab[colnames].copy()
    tab.columns = [x[:-4] for x in colnames]  # remove _TPM suffix
    tab = pd.melt(tab, value_vars=tab.columns)
    # plot
    plt.figure(figsize=(14,10))
    ax = sns.boxplot(
        data=tab, x='variable', y='value', notch=True,
        flierprops=dict(
            marker='.', markersize=5, markerfacecolor='steelblue',
            markeredgecolor='none', alpha=.4
        ),
        boxprops=dict(alpha=.7), width=.7
    )
    # set ylim
    hival = max(
        _calc_upper_whisker(
            tab.loc[tab['variable'] == name, 'value']
        ) 
        for name in tab['variable'].unique()
    )
    plt.ylim(-5, hival+5)
    # add noise to outliers.
    for line in ax.get_lines()[5::6]:
        xoffsets = line.get_xdata()
        line.set_xdata(xoffsets + np.random.uniform(-0.1, 0.1, xoffsets.size))
    # labelling
    plt.xticks(rotation=65)
    plt.title('X-RNA TPM coverage distributions for all biological replicas',
              x=0.5, y=1.0, ha='center', fontsize=20)
    plt.ylabel("TPM values")
    plt.xlabel("Experiment / cell line")
    # save
    plt.savefig(out_dir / f'{prefix}_tpm_distr.png',
                dpi=600, bbox_inches='tight')
