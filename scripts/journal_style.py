"""Journal-wide quantitative-plot grammar at 180 mm final figure width.

The colors identify comparisons; axes, intervals and labels carry the evidence.
Vector schematics keep their separate source-family palette.
"""
INK = '#20272D'
MID = '#5C6871'
GRID = '#E5E9EC'
BLUE = '#2D6998'
TEAL = '#08796D'
AMBER = '#AD762A'
PURPLE = '#806A91'
ROSE = '#A6536D'
MODEL_COLORS = [MID, BLUE, AMBER, PURPLE, '#416B6C', '#55779B', ROSE, TEAL, '#89959D']
LINE_STYLES = ['-', '--', '-.', ':', (0, (5, 1, 1, 1)), (0, (3, 1)), (0, (1, 1)), '-']

def apply(plt):
    plt.rcParams.update({
        'font.family':'Arial', 'font.size':9, 'axes.labelsize':9,
        'axes.titlesize':9.3, 'axes.titleweight':'bold',
        'xtick.labelsize':8.1, 'ytick.labelsize':8.1,
        'mathtext.fontset':'stix', 'axes.formatter.use_mathtext':True,
        'axes.linewidth':.65, 'axes.spines.top':False, 'axes.spines.right':False,
        'axes.edgecolor':'#657078', 'axes.labelcolor':INK, 'text.color':INK,
        'xtick.color':MID, 'ytick.color':MID, 'xtick.major.width':.6,
        'ytick.major.width':.6, 'xtick.major.size':3, 'ytick.major.size':3,
        'xtick.direction':'out', 'ytick.direction':'out',
        'axes.axisbelow':True, 'axes.titlepad':9, 'axes.labelpad':6,
        'grid.color':GRID, 'grid.linewidth':.45, 'lines.solid_capstyle':'round',
        'lines.linewidth':1.1, 'lines.markersize':4,
        'legend.frameon':False, 'legend.fontsize':8.1,
        'legend.handlelength':1.6, 'legend.labelspacing':.35,
        'figure.facecolor':'white', 'axes.facecolor':'white', 'savefig.facecolor':'white',
        'svg.fonttype':'none', 'pdf.fonttype':42,
    })


def polish_figure(fig):
    """Normalize axes without changing observations, limits or error intervals."""
    for ax in fig.axes:
        if not ax.axison:
            continue
        ax.set_axisbelow(True)
        for side in ('left', 'bottom'):
            ax.spines[side].set_color('#657078')
            ax.spines[side].set_linewidth(.65)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis='both', which='major', direction='out',
                       width=.6, length=3, color='#657078', labelcolor=INK)
        for gridline in (*ax.get_xgridlines(), *ax.get_ygridlines()):
            if gridline.get_visible():
                gridline.set_color(GRID)
                gridline.set_linewidth(.45)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)

def normalize_svg(path):
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
