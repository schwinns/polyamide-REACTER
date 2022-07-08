# Script to quickly plot thermodynamic data from lammps log file

import pandas as pd
import matplotlib.pyplot as plt
from thermo import get_thermo
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-l', '--log', default='log.lammps', help='log file with thermodynamic data')
parser.add_argument('-x','--x', default='step', help='property to plot on x axis')
parser.add_argument('-y','--y', default='pe',nargs='+',
                    help='property to plot on y axis')
parser.add_argument('-xi','--xmin', default=None, type=float, help='minimum x value for plotting')
parser.add_argument('-xa','--xmax', default=None, type=float, help='maximum x value for plotting')
parser.add_argument('-yi','--ymin', default=None, type=float, help='minimum y value for plotting')
parser.add_argument('-ya','--ymax', default=None, type=float, help='maximum y value for plotting')
parser.add_argument('-o', '--options', action='store_true', help='show possible properties to plot')
parser.add_argument('-s', '--save', action='store_true', help='save a png of the figure')
args = parser.parse_args()

df = get_thermo(args.log)

if args.options:
    print('Possible properties: {}'.format(df.columns.to_list()))

df.plot(x=args.x, y=args.y, xlim=(args.xmin,args.xmax), ylim=(args.ymin, args.ymax))

if args.save:
    fig_name = ''
    for y in args.y:
        fig_name += y

    fig_name += '_vs_' + args.x + '.png'
    plt.savefig(fig_name)
plt.show()