from Funcs import MakeGraph
import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import torch
from torch_geometric.data import Dataset
from torch_geometric.loader import DataLoader


class NuclearSubtypeDataset(Dataset):

    __metaclass__ = MakeGraph

    def __init__(self, root, transform=None, pre_transform=None, pre_filter=None):
        super(NuclearSubtypeDataset, self).__init__(root, transform, pre_transform, pre_filter)

    # create a generator object to retrieve filenames within a directory that have a certain filetype
    def retrieve_file_path(self, directory, filetype):
        filepaths = []
        for file in os.listdir(directory):
            if os.path.isfile(os.path.join(directory, file)):
                if file.endswith(f'.{filetype}'):
                    filepaths.append(file.split('.csv')[0])
        return filepaths

    # retrieve raw and processed file path names
    @property
    def raw_file_names(self):
        return self.retrieve_file_path(f'{self.root}/raw', 'csv')

    @property
    def processed_file_names(self):
        return self.retrieve_file_path(f'{self.root}/processed', 'pt')

    def download(self):
        pass

    def process(self):
        for idx, matrix_file in enumerate(self.retrieve_file_path(f'{self.root}/raw', 'csv')):
            graph = MakeGraph.MakeGraph(f'{self.root}/raw',
                                        matrix_file,
                                        ['Stain 1 Nucleus OD', 'Stain 2 Cytoplasm OD'],
                                        75).graph(torch_graph=True)

            if self.pre_filter is not None and not self.pre_filter(graph):
                continue

            if self.pre_transform is not None:
                graph = self.pre_transform(graph)

            torch.save(graph,
                       f'{self.root}/processed/graph_{idx}.pt')

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        graph = torch.load(f'{self.root}/processed/graph_{idx}.pt')
        return graph


NuclearSubtypeDataset(root='/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data')

x = DataLoader(NuclearSubtypeDataset(root='/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data'),
               batch_size=2,
               shuffle=True)
