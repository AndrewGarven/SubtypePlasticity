from Funcs import MakeGraph
import os
from collections import defaultdict
import json
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
import torch
from torch_geometric.data import Dataset
from torch_geometric.loader import DataLoader


class NuclearSubtypeDataset(Dataset):

    """
    NuclearSubtypeDataset is a torch_geometric.data.Dataset generated dynamically by torch_geometric.loader.DataLoader()
    used in torch_geometric GNN classification.

        Typical uses:

            Dataset : NuclearSubtypeDataset(root='path/to/root')

                        - This will dynamically generate graph_{idx}.pt file for any new 'root/raw/*.csv' files.
                            - These graphs are generated from .csv files in 'root/raw' according to MakeGraph().graph()
                              function
                            - This function will simply return requested torch_geometric.data.Dataset when all
                              'root/raw/*.csv' files have generated 'root/processed/*.pt' files.

            DataLoader : torch_geometric.loader.DataLoader(NuclearSubtypeDataset(root='path/to/root'),
                                                           batch_size=batch_size)



        :root: [str] path to root directory which contains 'raw' and 'processed' subdirectories
            ex. '~/SubtypePlasticity/Data/'
    """

    __metaclass__ = MakeGraph

    def __init__(self, root, transform=None, pre_transform=None, pre_filter=None):
        super(NuclearSubtypeDataset, self).__init__(root, transform, pre_transform, pre_filter)

    # create a generator object to retrieve filenames within a directory that have a certain filetype
    def retrieve_file_path(self, directory, filetype):

        """
            retrieve_file_path() generates list of files within 'directory' that have file type 'filetype'

            :self: NuclearSubtypeDataset Object described above
            :directory: [str] path to subdirectory of 'root'
            :filetype: [str] type of files you wish generate list from ('csv', 'pt')
            :return: list containing networkx node data [(node id, {weight: [], size: [], pos: ()}), ...]
        """

        filepaths = []
        for file in os.listdir(directory):
            if os.path.isfile(os.path.join(directory, file)):
                if file.endswith(f'.{filetype}'):
                    filepaths.append(file.split(f'.{filetype}')[0])
        return filepaths

    # retrieve raw and processed file path names
    @property
    def raw_file_names(self):
        return self.retrieve_file_path(f'{self.root}/raw', 'csv')

    @property
    def processed_file_names(self):
        return self.retrieve_file_path(f'{self.root}/processed', 'pt')

    # option to download object url data in torch_geometric.data.Dataset() [passed]
    def download(self):
        pass

    def process(self):

        """
            process() generates a 'processed/*.pt' file from 'raw/*.csv' file use MakeGraph().graph() function

            :self: NuclearSubtypeDataset Object described above
            :return: '.pt' file saved in 'root/processed/*.pt'
        """
        idx_id = defaultdict(dict)
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

            idx_id[f'{str(idx)}']: matrix_file
        idx_id = dict(idx_id)
        with open('PositionInputData/idx_id.json', 'wb') as fp:
            json.dump(idx_id, fp)

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        """
            get() generates a torch_geoetric.data.data.Data() object of 'root/processed/graph_{idx}.pt'

            :self: NuclearSubtypeDataset Object described above
            :idx: index generated by
            :return: '.pt' file saved in 'root/processed/*.pt'
        """
        graph = torch.load(f'{self.root}/processed/graph_{idx}.pt')
        return graph


#NuclearSubtypeDataset(root='/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data')

x = DataLoader(NuclearSubtypeDataset(root='/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data'),
               batch_size=2)

print(x)
