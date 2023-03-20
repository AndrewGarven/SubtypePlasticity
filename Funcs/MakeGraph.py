import pandas as pd
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from math import sqrt
from os import path, mkdir
#from torch_geometric.utils.convert import from_networkx


class MakeGraph:

    """MakeGraph is used to generate a networkx object from list of 2D coordinates with nodal information

        Typical use:

            graph : MakeGraph(datadir, filename, nodeweight, maxlinkingdistance).graph()

        :datadir: string of directory containing nuclear position data
        :filename: string of (.csv) file containing nuclear positional information
    *** :nodeweight: list of column name(s) strings in nuclear position dataframe (in 'filename'.csv) that contain node
                     information (i.e. nuclear staining intensity, nuclear optimal density, cellular h-score value)
                     - the column name(s) list must contain less than 2 column names
        :maxlinkingdistance(int): values of distance between nuclear centroids (in pixels) to form a link (edge) between
        nuclei

        """

    def __init__(self, datadir, filename, nodeweight, maxlinkingdistance):
        # Initializing preprocessing class object
        # automatically generates a pandas dataframe
        self.datadir = datadir
        self.filename = filename
        self.df = pd.read_csv(f'{datadir}/{filename}.csv', index_col=0)
        self.ld = maxlinkingdistance
        self.nodeweight = nodeweight

    def addedges(self):

        """
            addedges() generates networkx edge data linking nuclei that are 'linkingdistance' apart

            :self: Preprocessing Object described above
            :return: list containing networkx edge data [(node id1, node id2, lineardistance), ...]
                            """
        # create a list of all possible x and y pairs
        xlist = zip(self.df.index.values.tolist(), self.df['XCentroid'])
        ylist = zip(self.df.index.values.tolist(), self.df['YCentroid'])
        xpairs = [x for x in itertools.product(xlist, repeat=2)]
        ypairs = [y for y in itertools.product(ylist, repeat=2)]
        centroidpairs = zip(xpairs, ypairs)

        edges = []

        for pair in centroidpairs:
            # retrieve id data from each x, y pair
            ids = pair[0][0][0], pair[0][1][0]
            # calculate the linear distance between the centroid of two nodes in a 2D plane
            dist = sqrt(((pair[0][1][1] - pair[0][0][1])**2) + ((pair[1][1][1] - pair[1][0][1])**2))
            if ids[0] == ids[1]:
                continue
            # generate 'edges' for each pair that are closer than the set 'maxlinkingdistance'
            elif 0 < dist < self.ld:
                edge = (str(ids[0]), str(ids[1]), {'length': dist})
                edges.append(edge)
            else:
                continue

        return edges

    def addnodes(self):

        """
            addnodes() generates networkx node data from nuclear position and staining intensity values.

            :self: Preprocessing Object described above
            :return: list containing networkx node data [(node id, {weight: [], size: [], pos: ()}), ...]
                    """

        # calculate centroid coordinates from [[Xmin, Xmax], [Ymin, Ymax]] values and store in self.df & in 'positions'
        self.df['XCentroid'] = (self.df.XMax + self.df.XMin) / 2
        self.df['YCentroid'] = (self.df.YMax + self.df.YMin) / 2

        positions = self.df[['XCentroid', 'YCentroid']].values.tolist()

        # store nodeweight(s) ('weight' and 'size') and position data for each node in a
        # list of nodes of form: [((ID),{'weight':[], ***'size':[], 'pos':[]}), ...]
        nodes = []
        if len(self.nodeweight) == 2:
            stains = self.df[[self.nodeweight[0], self.nodeweight[1]]].values.tolist()
            nodeids = self.df.index.values.tolist()
            nodedata = zip(nodeids, stains, positions)
            for nodeid, stain, position in nodedata:
                weight1 = stain[0]
                weight2 = stain[1]
                node = (str(nodeid), {'weight': weight1, 'size': weight2, 'pos': (int(position[0]), int(position[1]))})
                nodes.append(node)
        elif len(self.nodeweight) == 1:
            stains = self.df[[self.nodeweight[0]]].values.tolist()
            nodeids = self.df.index.values.tolist()
            nodedata = zip(nodeids, stains, positions)
            for nodeid, stain, position in nodedata:
                weight1 = stain[0]
                node = (str(nodeid), {'weight': weight1, 'pos': (int(position[0]), int(position[1]))})
                nodes.append(node)
        else:
            pass

        return nodes

    def graph(self, image='no'):

        """
            graph generates a networkx.Graph() from nuclear position (with staining intensity) 'Preprocessing' object.

            :self: Preprocessing Object described above
            :image: string of 'yes' or 'no' based on if you want an image of your graph returned

            :return: 1. networkx.Graph() representing nuclear position (with staining intensity)
                     2. PNG image file of 2D graph network with node weights represented as color and/or size of node
                        image saved to 'datadir'/'filename'.png
            """
        # generate networkx.Graph() object and load in our generated nodes and edges data
        graph = nx.Graph()
        graph.add_nodes_from(self.addnodes())
        graph.add_edges_from(self.addedges())
        #geograph = from_networkx(graph)

        if image == 'yes' or image == 'Yes':

            pos = nx.get_node_attributes(graph, 'pos')
            outpath = f'{self.datadir}/GraphImages/{self.filename}.png'

            if len(self.nodeweight) == 2:
                # if there are two variables associated with a single node (i.e. GATA3 & CK5 expression),
                # their data can be stored in two variable types 'size' and 'color' of the node
                size = nx.get_node_attributes(graph, 'size')
                weight = nx.get_node_attributes(graph, 'weight')
                # graph is displayed, saved, and returned as networkx.Graph object
                nx.draw(graph, pos, node_color=list(weight.values()), node_size=list(size.values()))
                if path.exists(f'{self.datadir}/GraphImages'):
                    plt.savefig(outpath, format="PNG")
                    return graph
                else:
                    mkdir(f'{self.datadir}/GraphImages')
                    plt.savefig(outpath, format="PNG")
                    return graph

            elif len(self.nodeweight) == 1:
                # if there is only one variable associated with each node (i.e. p16 expression)
                # the data will be expressed as a variance in 'color' between each node in the graph
                weight = nx.get_node_attributes(graph, 'weight')
                # graph is displayed, saved, and returned as networkx.Graph object
                nx.draw(graph, pos, node_color=list(weight.values()))

                if path.exists(f'{self.datadir}/GraphImages'):
                    plt.savefig(outpath, format="PNG")
                    return graph
                else:
                    mkdir(f'{self.datadir}/GraphImages')
                    plt.savefig(outpath, format="PNG")
                    return graph

            else:
                pass

        else:
            return graph