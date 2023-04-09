import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import pandas as pd
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from math import sqrt
from torch_geometric.utils import from_networkx


class MakeGraph:

    """
    MakeGraph is used to generate a networkx object from list of 2D coordinates with nodal information

        Typical use:

            graph : MakeGraph(data_dir, file_name, node_weight, max_linking_distance).graph()

        :datadir: string of directory containing nuclear position data
        :filename: string of (.csv) file containing nuclear positional information
    *** :node_weight: list of column name(s) strings in nuclear position dataframe (in 'filename'.csv) that contain node
                     information (i.e. nuclear staining intensity, nuclear optimal density, cellular h-score value)
                     - the column name(s) list must contain less than 2 column names
        :max_linking_distance: (int) value of distance between nuclear centroids (in pixels) to form a link (edge)
                              between nuclei
    """

    def __init__(self, data_dir, file_name, node_weight, max_link_distance):
        # Initializing preprocessing class object
        # automatically generates a pandas dataframe
        self.data_dir = data_dir
        self.file_name = file_name
        self.df = pd.read_csv(f'{data_dir}/{file_name}.csv', index_col=0)
        self.node_weight = node_weight
        self.node_feature_count = len(self.node_weight)
        self.ld = max_link_distance

    def add_edges(self):

        """
            add_edges() generates networkx edge data linking nuclei that are < 'max_linking_distance' apart

            :self: Preprocessing Object described above
            :return: list containing networkx edge data [(node id1, node id2, linear_distance), ...]
        """

        # create a list of all possible x and y pairs
        x_list = zip(self.df.index.values.tolist(), self.df['XCentroid'])
        y_list = zip(self.df.index.values.tolist(), self.df['YCentroid'])
        x_pairs = [x for x in itertools.product(x_list, repeat=2)]
        y_pairs = [y for y in itertools.product(y_list, repeat=2)]
        centroid_pairs = zip(x_pairs, y_pairs)

        edges = []

        for pair in centroid_pairs:
            # retrieve id data from each x, y pair
            id1, id2 = pair[0][0][0], pair[0][1][0]

            if id1 == id2:
                continue
            else:
                # calculate the linear distance between the centroid of two nodes in a 2D plane
                dist = sqrt(((pair[0][1][1] - pair[0][0][1]) ** 2) + ((pair[1][1][1] - pair[1][0][1]) ** 2))
                # generate 'edges' for each pair that are closer than the set 'max_linking_distance'
                if dist < self.ld:
                    edge = (str(id1), str(id2), {'length': dist})
                    edges.append(edge)

        return edges

    def add_nodes(self):

        """
            add_nodes() generates networkx node data from nuclear position and staining intensity values.

            :self: Preprocessing Object described above
            :return: list containing networkx node data [(node id, {weight: [], size: [], pos: ()}), ...]
        """

        # calculate centroid coordinates from [[Xmin, Xmax], [Ymin, Ymax]] values and store in self.df & in 'positions'
        self.df['XCentroid'] = (self.df.XMax + self.df.XMin) / 2
        self.df['YCentroid'] = (self.df.YMax + self.df.YMin) / 2

        positions = self.df[['XCentroid', 'YCentroid']].values.tolist()

        # store nodeweight(s) ('weight' and 'size') and position data for each node in a
        # list of nodes of form: [((ID),{'weight':[], ***'size':[], 'pos':[]}), ...]

        if self.node_feature_count == 2:
            stains = self.df[[self.node_weight[0], self.node_weight[1]]].values.tolist()
        else:
            stains = self.df[self.node_weight[0]].values.tolist()

        node = self.df.index.values.tolist()
        node_data = zip(node, stains, positions)

        nodes = []

        for node_id, stain, position in node_data:
            weight1 = stain[0]
            if self.node_feature_count == 2:
                weight2 = stain[1]
                node = (str(node_id), {'weight': weight1, 'size': weight2, 'pos': (int(position[0]), int(position[1]))})
            else:
                node = (str(node_id), {'weight': weight1, 'pos': (int(position[0]), int(position[1]))})
            nodes.append(node)

        return nodes

    def graph(self, image=False, image_path=None, torch_graph=False):

        """
            graph generates a 'torch_geometric.data.data.Data' object from nuclear position (with staining intensity)
            if torch parameter is set to True, otherwise, graph generates a networkx.Graph() object.

            :self: Preprocessing Object described above
            :image: set True if you want an image of networkx graph
            :image_path: string with path to directory to save graph image
            :torch_graph: set True if you want a torch_geometric.data.data.Data() object instead of networkx.Graph()

            :return: 1. (default): networkx.Graph()
                        (if torch_graph = True): torch_geometric.data.data.Data()

                     2. (if image = True): PNG image file of 2D graph network with node weights represented as color
                        and/or size of node image saved to (default): 'datadir'/'filename'.png
                                                           (if image_path specified): :image_path:/'filename'.png

            """

        # generate networkx.Graph() object and load in our generated nodes and edges data
        graph = nx.Graph()
        graph.add_nodes_from(self.add_nodes())
        graph.add_edges_from(self.add_edges())

        if image:

            # images require both node position and 'weight' data ('size' captured below if len(:node_weights:) > 1
            pos = nx.get_node_attributes(graph, 'pos')
            weight = nx.get_node_attributes(graph, 'weight')

            if image_path:
                out_path = image_path + f'/{self.file_name}.png'
            else:
                out_path = f'{self.data_dir}/{self.file_name}.png'

            # we will capture node feature information as colour (node weight) [and size (node size)]
            if self.node_feature_count == 2:
                size = nx.get_node_attributes(graph, 'size')
                nx.draw(graph, pos, node_color=list(weight.values()), node_size=list(size.values()))
            elif self.node_feature_count == 1:
                nx.draw(graph, pos, node_color=list(weight.values()))

            plt.savefig(out_path)

        if torch_graph:
            # generate a torch_geometric.data.data.Data object if :torch_graph: = True
            graph = from_networkx(graph)

        return graph
