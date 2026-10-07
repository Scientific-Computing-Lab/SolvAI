"""Shared publication typography and color semantics (180 mm full width)."""
INK = '#18303C'
MID = '#647581'
GRID = '#E3E9EC'
BLUE = '#417FA8'
TEAL = '#008C7A'
AMBER = '#C78935'
PURPLE = '#9870AD'
ROSE = '#B05C7C'
MODEL_COLORS = [MID, BLUE, AMBER, PURPLE, '#277476', '#5173B4', ROSE, TEAL, '#8097A5']
LINE_STYLES = ['-', '--', '-.', ':', (0, (5, 1, 1, 1)), (0, (3, 1)), (0, (1, 1)), '-']

def apply(plt):
    plt.rcParams.update({
        'font.family':'Arial', 'font.size':8.5, 'axes.labelsize':8.5,
        'axes.titlesize':10, 'axes.titleweight':'bold',
        'xtick.labelsize':8, 'ytick.labelsize':8,
        'axes.linewidth':.65, 'axes.spines.top':False, 'axes.spines.right':False,
        'axes.edgecolor':MID, 'axes.labelcolor':INK, 'text.color':INK,
        'xtick.color':MID, 'ytick.color':MID, 'xtick.major.width':.65,
        'ytick.major.width':.65, 'xtick.major.size':3, 'ytick.major.size':3,
        'axes.axisbelow':True, 'axes.titlepad':9, 'axes.labelpad':6,
        'figure.facecolor':'white', 'axes.facecolor':'white', 'savefig.facecolor':'white',
        'svg.fonttype':'none', 'pdf.fonttype':42,
    })

def normalize_svg(path):
    path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
