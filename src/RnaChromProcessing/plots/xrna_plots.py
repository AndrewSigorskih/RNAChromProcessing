from contextlib import contextmanager
from pathlib import Path
from typing import Generator, Sequence

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns

FIGSIZE = (11.7, 8.27)
TEN_KB = 10_000


_FIGSIZE = (7.0, 4.5)  # inches; reads well at single/1.5-column width
BAR_COLOR = '#595959'     # ggplot2's default histogram grey
# Paul Tol "bright" blue/red: colourblind-safe, and unlike complementary
# pairs (blue/orange) their translucent overlap reads as a distinct mauve
# rather than muddy grey
PRIME_PALETTE = {"5'": '#4477AA', "3'": '#EE6677'}
HIST_ALPHA = 0.3

_GGPLOT_RC = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.titleweight': 'bold',
    'axes.titlelocation': 'left',
    'axes.titlepad': 10,
    'axes.labelsize': 11,
    'axes.labelweight': 'normal',
    'axes.labelpad': 6,
    'axes.edgecolor': '#333333',
    'axes.linewidth': 0.8,
    'axes.grid': True,
    'axes.grid.axis': 'y',
    'axes.axisbelow': True,
    'grid.color': '#E5E5E5',
    'grid.linewidth': 0.7,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'xtick.color': '#333333',
    'ytick.color': '#333333',
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'legend.fontsize': 9,
    'legend.title_fontsize': 10,
    'legend.frameon': False,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    # keep text editable in vector output (Illustrator/Inkscape)
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
}


@contextmanager
def ggplot_style() -> Generator[None, None, None]:
    """Apply the style only inside the block, leaving global
    rcParams (e.g. those set by set_style_white) untouched afterwards."""
    with sns.axes_style('ticks'), plt.rc_context(_GGPLOT_RC):
        yield


def set_style_publication() -> None:
    """Apply the publication style globally (alternative to set_style_white)."""
    sns.set_style('ticks')
    plt.rcParams.update(_GGPLOT_RC)



def _thousands(axis) -> None:
    axis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))


def _save(fig: plt.Figure, out_dir: Path, stem: str,
          formats: Sequence[str]) -> None:
    for ext in formats:
        fig.savefig(out_dir / f'{stem}.{ext}')
    plt.close(fig)


# --- plots -----------------------------------------------------------------

def plot_length_distribution(tab: pd.DataFrame,
                             out_dir: Path,
                             prefix: str,
                             formats: Sequence[str] = ('png',)) -> None:
    length = tab['end'] - tab['start'] + 1
    shown = length[length <= TEN_KB]

    with ggplot_style():
        fig, ax = plt.subplots(figsize=_FIGSIZE)
        ax.hist(shown, bins=20, color=BAR_COLOR,
                edgecolor='white', linewidth=0.6)

        ax.set_title('X-RNA length distribution')
        ax.set_xlabel('Length (bp)')
        ax.set_ylabel('Count')
        _thousands(ax.xaxis)
        _thousands(ax.yaxis)
        ax.set_ylim(bottom=0)
        ax.margins(x=0.01)
        sns.despine(ax=ax)

        # summary stats are over *all* X-RNAs, not just the ≤10 kb shown
        info = (
            f'All X-RNAs (n = {len(length):,})\n'
            f'Min: {length.min():,.0f} bp\n'
            f'Max: {length.max():,.0f} bp\n'
            f'Median: {length.median():,.1f} bp'
        )
        ax.text(0.98, 0.96, info, transform=ax.transAxes,
                ha='right', va='top', fontsize=9, linespacing=1.6,
                bbox=dict(boxstyle='round,pad=0.6', facecolor='white',
                          edgecolor='#CCCCCC', linewidth=0.6))

        _save(fig, out_dir, f'{prefix}_length_distr', formats)


def plot_distance_to_closest(tab: pd.DataFrame,
                             out_dir: Path,
                             prefix: str,
                             formats: Sequence[str] = ('png',)) -> None:
    # select 5' and 3' distances that are within 10kb window around genes
    data = tab[['closest_gene_side', 'closest_gene_dist']].copy()
    data.columns = ['Side', 'Distance']
    data['Distance'] = data['Distance'].abs()
    data = data[data['Distance'] <= TEN_KB]
    present = set(data['Side'])
    hue_order = [s for s in PRIME_PALETTE if s in present]

    with ggplot_style():
        fig, ax = plt.subplots(figsize=_FIGSIZE)
        ax2 = ax.twinx()

        sns.histplot(
            data=data, x='Distance', hue='Side', ax=ax,
            hue_order=hue_order, palette=PRIME_PALETTE,
            element='step', alpha=HIST_ALPHA, linewidth=1.0,
        )
        sns.kdeplot(
            data=data, x='Distance', hue='Side', ax=ax2,
            hue_order=hue_order, palette=PRIME_PALETTE,
            clip=(0.0, None), linewidth=2.0, legend=False,
        )
        ax2.grid(False)

        # density curves sit on the twin axis (drawn on top), so move the
        # legend there too, otherwise the curves would run over it
        leg = ax.get_legend()
        if leg is not None:
            handles = getattr(leg, 'legend_handles', None) or leg.legendHandles
            labels = [t.get_text() for t in leg.get_texts()]
            leg.remove()
            ax2.legend(handles, labels, title='Closest gene end',
                       loc='upper right', frameon=False)

        ax.set_title('Distance to closest gene by gene end')
        ax.set_xlabel('Distance (bp)')
        ax.set_ylabel('Count')
        ax2.set_ylabel('Density')
        _thousands(ax.xaxis)
        _thousands(ax.yaxis)
        ax2.yaxis.set_major_formatter(
            mticker.ScalarFormatter(useMathText=True))
        ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        # data are filtered to <=10 kb; don't let the KDE tail stretch the axis
        ax.set_xlim(0, TEN_KB)
        ax.set_ylim(bottom=0)
        ax2.set_ylim(bottom=0)
        for a in (ax, ax2):
            a.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        _save(fig, out_dir, f'{prefix}_closest_gene_distances', formats)
