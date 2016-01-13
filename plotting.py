'''
Functions for plotting graphs of data

Created on 14/01/2016

@author: iwelsh
'''

import matplotlib.pyplot as plt
import os

def plot_scatters(data, save_path, default_point='.', show_legend=True, xlim=None, ylim=None,
                  x_label='X axis', y_label='Y axis', title='Graph title', file_out='png'):
    """
    Plot data on a scatter plot, with title, axis labels and legend.
    data is a list of dicts containing the data, and any information relating to that data. Each dict must have:
    x the list of x values to plot
    y the list of y values to plot
    and can have:
    label the data set label
    marker what marker colour etc to use
    
    """
    if not save_path.endswith('.'+file_out):
        save_path += '.'+file_out
    if not os.path.isdir('/'.join(save_path.split('/')[:-1])):
        os.makedirs('/'.join(save_path.split('/')[:-1]))
    # plot all the data
    for d_set in data:
        if isinstance(d_set['x'], dict):
            x = [d_set['x'][x] for x in sorted(d_set['x'])]
            y = [d_set['y'][y] for y in sorted(d_set['y'])]
        else:
            x = d_set['x']
            y = d_set['y']
        plt.plot(x, y, default_point if 'marker' not in d_set else d_set['marker'],
                 label=None if 'label' not in d_set else d_set['label'])
    # layout the graph nicely
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    if show_legend:
        plt.legend(loc=8, fontsize='xx-small', numpoints=1, mode='expand',
                   bbox_to_anchor=(0, 0.02, 1, 0.02), ncol=5)
    plt.savefig(save_path, format=file_out, papertype='a4')
    plt.close()
    
        

def main():
    save_path = '/Users/iwelsh/Documents/AdditionRigidTorsions/test.abc'
    x_data = [1,2,3,4,5,6,7,8,9,10]
    y_data = [1,4,9,16,25,36,49,64,81,100]
    data = [{'x':x_data,'y':y_data,'marker':'^--'}]
    x_label = r'$x$'
    y_label = r'$x^2$'
    title = 'Quadratics for 1 to 10'
    plot_scatters(data, save_path, show_legend=False, x_label=x_label, y_label=y_label, title=title, file_out='pdf')

if __name__ == '__main__':
    main()