import pandas as pd
import itertools
import networkx as nx
import matplotlib.pyplot as plt
from math import sqrt


class Preprocess:

    def __init__(self, datadir, filename):
        # Initializing preprocessing class object
        # automatically generates a pandas dataframe
        self.df = pd.read_csv(f"{datadir}" + '/' + f"{filename}" + '.csv', index_col=0)

    def calcdist(self):

        self.df['XCentroid'] = (self.df.XMax + self.df.XMin) / 2
        self.df['YCentroid'] = (self.df.YMax + self.df.YMin) / 2

        # inputs
        # x = list of x centroid coordinates
        # y = list of x centroid coordinates
        xlist = zip(self.df.index.values.tolist(), self.df['XCentroid'])
        ylist = zip(self.df.index.values.tolist(), self.df['YCentroid'])
        xpairs = [x for x in itertools.product(xlist, repeat=2)]
        ypairs = [y for y in itertools.product(ylist, repeat=2)]
        centroidpairs = zip(xpairs, ypairs)

        edges = []

        for pair in centroidpairs:
            ids = pair[0][0][0], pair[0][1][0]
            # calculate the linear distance between the centroid of two nuclei
            dist = sqrt(((pair[0][1][1] - pair[0][0][1])**2) + ((pair[1][1][1] - pair[1][0][1])**2))
            if ids[0] == ids[1]:
                continue
            elif 0 < dist < 75:
                edge = [ids[0], ids[1], dist]
                edges.append(edge)
            else:
                continue

        return edges

    def addnodes(self):

        nodes = self.df[['XCentroid', 'YCentroid']].values.tolist()
        weights = self.df[['Stain 2 Cytoplasm OD']].values.tolist()
        size = self.df[['Stain 1 Nucleus OD']].values.tolist()

        return nodes, weights, size

    def graph(self):
        edges = self.calcdist()
        nodes, weights, size = self.addnodes()
        graph = nx.Graph(nodes=nodes,
                         e=edges,
                         weights=weights,
                         size=size)
        #graph.add_nodes_from(nodes, weights=weights, size=size)
        #graph.add_weighted_edges_from(edges)
        nx.draw(graph)
        plt.show()


Preprocess('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data', 'TEST-tmaCoordinates').graph()
