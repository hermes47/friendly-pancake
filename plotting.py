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
    ax = plt.subplot(111)
    # plot all the data
    for d_set in data:
        if isinstance(d_set['x'], dict):
            x = [d_set['x'][x] for x in sorted(d_set['x'])]
            y = [d_set['y'][y] for y in sorted(d_set['y'])]
        else:
            x = d_set['x']
            y = d_set['y']
        ax.plot(x, y, default_point if 'marker' not in d_set else d_set['marker'],
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
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        plt.legend(loc='upper center', fontsize='xx-small', numpoints=1,
                   bbox_to_anchor=(0.5, -0.12), ncol=5)
    plt.savefig(save_path, format=file_out, papertype='a4')
    plt.close()

def main():
    pass
if __name__ == '__main__':
    main()
