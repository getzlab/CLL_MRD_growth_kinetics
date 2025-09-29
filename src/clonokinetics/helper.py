import pandas as pd


def plot_ccf(df, ax, times_sample):
    import matplotlib.pyplot as plt
    # Keep the necessary columns
    cols = ['Sample_ID', 'Cluster_ID', 'postDP_ccf_mean', 'postDP_ccf_CI_low', 'postDP_ccf_CI_high']
    df = df[cols]
    cluster_list = df.Cluster_ID.unique().tolist()
    number_samples = len(df.Sample_ID.unique())

    tick_list = ['T' + str(i) for i in range(number_samples)]
    x_axis = [i / 365 for i in times_sample]
    ax.set_xticks(x_axis)

    secax = ax.secondary_xaxis('top')
    secax.set_xlabel("Time (years)")
    ax.grid(True)

    for i in cluster_list:
        x = df[df.Cluster_ID == i].Sample_ID
        y = df[df.Cluster_ID == i].postDP_ccf_mean
        ci_low = df[df.Cluster_ID == i].postDP_ccf_CI_low
        ci_high = df[df.Cluster_ID == i].postDP_ccf_CI_high

        #         if x_scale == 'sample':
        #             x_axis = np.arange(0,number_samples)

        #         else:
        #             x_axis = times_sample
        ax.plot(x_axis, y, c=ClusterColors.get_hex_string(i), marker='o', label=i)

        #         ax.plot(x_axis, y,c= ClusterColors.get_hex_string(i), marker ='o')
        ax.fill_between(x_axis, ci_low, ci_high, color=ClusterColors.get_hex_string(i), alpha=0.1)

        ax.set_xlabel('Samples')
        ax.set_xticks(x_axis)

        ax.set_xticklabels(tick_list, fontsize=8)

        ax.set_ylabel('CCF')
        ax.legend()

    cmap = plt.get_cmap("Pastel1")
    xlim = ax.get_xlim()[1]
    for i, row in treatment.iterrows():
        treatment_name = row.tx
        start = row.tx_start / 365
        end = row.tx_end / 365
        if np.isnan(end):
            end = xlim
        length = end - start
        center = (start + end) / 2

        #         ax.axvspan(xmin = start, xmax= end, label = treatment_name, facecolor= cmap(i), alpha = 0.2)
        ax.axvspan(xmin=start, xmax=end, facecolor=cmap(i), alpha=0.2)


#         ax.legend(ncol = treatment.shape[0], loc='upper center', bbox_to_anchor=(0.5, -2), fontsize = 'x-large')


def plot_tree(tree_posterior_file, edge_labels):
    tree_df = pd.read_csv(tree_posterior_file, sep='\t')
    edges = tree_df.loc[0, 'edges'].split(',')
    cluster_list = []
    for i in edges:
        new_list = i.split('-')
        for j in new_list:
            if (j != 'None') & (j not in cluster_list):
                cluster_list.append(j)
    cluster_list = [int(i) for i in cluster_list]

    DG = nx.DiGraph()
    for edge in edges:
        nodes = edge.split('-')
        if nodes[0] != 'None':
            DG.add_edge(int(nodes[0]), int(nodes[1]))

    pos = graphviz_layout(DG, prog='dot')

    edge_color_list = []
    for edge in DG.edges():
        node = edge[1]

        edge_color_list.append(ClusterColors.get_hex_string(node))

    nx.draw(DG, pos, with_labels=True, width=4, node_size=400, arrows=False, font_color='white',
            edge_color=edge_color_list, node_color=[ClusterColors.get_hex_string(i) for i in cluster_list])
    nx.draw_networkx_edge_labels(DG, pos,
                                 edge_labels,
                                 font_color=ClusterColors.get_hex_string(7))

    return DG


def simple_lin_interp(x1, y1, x2, y2, x):
    slope = (y2 - y1) / float(x2 - x1)
    return slope * (x - x1) + y1
def get_neighbor(DG, node):
    try:
        neighbors = DG.successors(DG.predecessors(node)[0])
    except:
        neighbors = []
    return neighbors
def return_smoothed(points, interpolate=False):
    if len(points) < 3: return points
    from scipy.interpolate import splrep
    from scipy.interpolate import splev
    from scipy.interpolate import Akima1DInterpolator

    x, y = zip(*points)

    if interpolate:
        return points
        # return zip(x,Akima1DInterpolator(x,y)(x))

    else:
        from statsmodels.nonparametric.kernel_regression import KernelReg
        kr = KernelReg(y, x, 'c')
        x_ = []  # list(np.linspace(min(x),max(x),3))
        x_.extend(x)
        x_ = sorted(set(x_))
        y_pred, _ = kr.fit(x_)

        return zip(x_, y_pred)


        if len(points) > 2:
            w = [1] + [0.1] * (len(points) - 2) + [1]
        else:
            w = [1] * len(points)

        K = min([len(points) - 1, 3])
        spl = splrep(x, y, s=30., t=[], k=K, w=w)
        # plt.figure()
        # plt.plot(x,splev(x,spl))
        # plt.plot(x,y,'o')
        # plt.show()
        x_ = []  # list(np.linspace(min(x),max(x),3))
        x_.extend(x)
        x_ = sorted(set(x_))
        return zip(x_, splev(x_, spl))
def nx_walk(node,tree,override_list=False):
    """ iterate tree in pre-order depth-first search order """
    if override_list is not False:
        for n in override_list:
            yield n
    else:
        yield node
        for child in tree.successors(node):
            for n in nx_walk(child,tree):
                yield n



class ClusterColors(object):
    # Cluster colors
    color_list = [[166, 17, 129],
                  [39, 140, 24],
                  [103, 200, 243],
                  [248, 139, 16],
                  [16, 49, 41],
                  [93, 119, 254],
                  [152, 22, 26],
                  [104, 236, 172],
                  [249, 142, 135],
                  [55, 18, 48],
                  [83, 82, 22],
                  [247, 36, 36],
                  [0, 79, 114],
                  [243, 65, 132],
                  [60, 185, 179],
                  [185, 177, 243],
                  [139, 34, 67],
                  [178, 41, 186],
                  [58, 146, 231],
                  [130, 159, 21],
                  [161, 91, 243],
                  [131, 61, 17],
                  [248, 75, 81],
                  [32, 75, 32],
                  [45, 109, 116],
                  [255, 169, 199],
                  [55, 179, 113],
                  [34, 42, 3],
                  [56, 121, 166],
                  [172, 60, 15],
                  [115, 76, 204],
                  [21, 61, 73],
                  [67, 21, 74],  # Additional colors, uglier and bad
                  [123, 88, 112],
                  [87, 106, 46],
                  [37, 66, 58],
                  [132, 79, 62],
                  [71, 58, 32],
                  [59, 104, 114],
                  [46, 107, 90],
                  [84, 68, 73],
                  [90, 97, 124],
                  [121, 66, 76],
                  [104, 93, 48],
                  [49, 67, 82],
                  [71, 95, 65],
                  [127, 85, 44],  # even more additional colors, gray
                  [88, 79, 92],
                  [220, 212, 194],
                  [35, 34, 36],
                  [200, 220, 224],
                  [73, 81, 69],
                  [224, 199, 206],
                  [120, 127, 113],
                  [142, 148, 166],
                  [153, 167, 156],
                  [162, 139, 145],
                  [0, 0, 0]]  # black

    @classmethod
    def get_rgb_string(cls, c):
        return 'rgb({},{},{})'.format(*cls.color_list[c])

    @classmethod
    def get_hex_string(cls, c):
        return '#{:02X}{:02X}{:02X}'.format(*cls.color_list[c])


class simple_tree():

    """
    Class to allow slicing points across intervals. Intended for large mafs.
    (works for 1d data, but not arbitrary length segments)
    Similar to a BST

    Built to work for int positions, but could be re-done to work

    """

    def __init__(self, _xvals=None, _data=None, sorted_=True):

        if (_data is not None) and (_xvals is not None):
            if sorted_:
                self._xvals = _xvals
                self._data = _data
            else:  # assume sorted, but if re-sorting is nessary do it.
                self._xvals, self._data = zip(*sorted(zip(_xvals, _data), key=lambda x: x[0]))
                self._xvals=list(self._xvals)
                self._data=list(self._data)
        else:
            self._xvals = []
            self._data = []

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return str(self._data)

    def __getitem__(self, sliced):
        if type(sliced) is int or type(sliced) is float:  # if the slice is a point, return the point
            return [self._xvals[bisect.bisect_left(self._xvals, sliced)], self._data[bisect.bisect_left(self._xvals, sliced)]]
        else:  # otherwise slice from the slice object with bisection
            return zip(self._xvals[bisect.bisect_left(self._xvals, sliced.start):bisect.bisect_left(self._xvals, sliced.stop)],
                       self._data[bisect.bisect_left(self._xvals, sliced.start):bisect.bisect_left(self._xvals, sliced.stop)])

    def add(self, item):
        """
        Add an item, assumes [pos(int),value]
        Append if not sorted, else binary search for position.
        """
        try:  # will throw index error on first attempt.

            _sorted = item[0] >= self._xvals[-1]
        except IndexError:
            _sorted = True  # it is slowish to check every time if the list exists or not, but on the first try it's going to give indexerror since it's empty.

        if _sorted:  # if item is at the end we can do O(1) append.
            self._xvals.append(item[0])
            self._data.append(item[1])

        else:  # if item is not at the end, we need  O(logN) binary search, then  O(N) insert. In all O(NlogN)

            ins_idx = bisect.bisect(self._xvals, item[0])
            self._xvals.insert(ins_idx, item[0])
            self._data.insert(ins_idx, item[1])

    def append(self, item):
        """
        Append an item, assumes [pos(int),value] and item belongs at end of list!
        """
        self._data.append(item[1])
        self._xvals.append(item[0])

    def copy(self):
        return simple_tree(_xvals=self._xvals[:], _data=self._data[:], sorted_=True)  # new tree, parameters should be default but might change.

    def _re_sort(self):
        """
        This function shouldn't need to be called. Re-sort the object.
        """
        self._xvals, self._data = zip(*sorted(zip(self._xvals, self._data), key=lambda x: x[0]))

