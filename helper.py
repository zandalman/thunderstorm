import matplotlib.pyplot as plt
from datetime import datetime
import os

# coordinates
X, Y, Z = 0, 1, 2

def save_fig(fig_name, filetype="png", dpi=256):
    '''
    Save the current matplotlib figure.

    Args
    name (string): figure name
    filetype (string): file type
    dpi (int): dots per inch
    '''
    datetime_string = datetime.now().strftime("%m%d%Y%H%M")
    filename = "%s-%s.%s" % (fig_name, datetime_string, filetype)
    plt.savefig(os.path.join('figures', filename), bbox_inches="tight", dpi=dpi)
    print("Saved figure as '%s'" % filename)