import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
from scipy.stats import pearsonr


class RNASeqExplore:

    def __init__(self, datadir, rna, clinical):
        self.rna = pd.read_csv(f'{datadir}/{rna}.csv', index_col=0)
        self.clinical = pd.read_csv(f'{datadir}/{clinical}.csv', index_col=0)

    def get_subset(self, subset_col, subset_names):

        subset_count = len(subset_names)

        subset1 = []
        subset2 = []
        subset3 = []
        subset4 = []
        subset5 = []

        # Get all shared IDs in both rna and clinical data
        shared_ids = list(set(self.rna.columns.values.tolist()).intersection(set(self.clinical.index.values.tolist())))
        for shared_id in shared_ids:
            info = self.clinical.loc[shared_id, subset_col]
            if subset_count >= 1 and info == subset_names[0]:
                subset1.append(shared_id)
            elif subset_count >= 2 and info == subset_names[1]:
                subset2.append(shared_id)
            elif subset_count >= 3 and info == subset_names[2]:
                subset3.append(shared_id)
            elif subset_count >= 4 and info == subset_names[3]:
                subset4.append(shared_id)
            elif subset_count >= 5 and info == subset_names[4]:
                subset5.append(shared_id)
            else:
                continue

        subset_list = [subset1, subset2, subset3, subset4, subset5][:subset_count]
        for i, subset in enumerate(subset_list):
            for patient_id in subset:
                self.rna.loc[f'{subset_col}', patient_id] = i

    def make_legend(self, label):

        label_conv = {
            'Ba/Sq': 'x',
            'Stroma-rich': 'o',
            'LumP': '^',
            'LumU': 's',
            'LumNS': '.'
        }

        lab_legend = mlines.Line2D([], [], color='black',
                                   marker=f'{label_conv[label]}',
                                   linestyle='None',
                                   markersize=5,
                                   label=f'{label}')

        return lab_legend

    def plot_pca(self, gene, subtype_col=None, subtype_names=None):
        self.rna.columns = self.rna.columns.astype(str)
        self.rna.index = self.rna.index.astype(str)
        pca = PCA(n_components=2)
        pca_transformed = pca.fit_transform(self.rna.transpose())
        gray = plt.get_cmap('gray')

        if subtype_col and subtype_names:

            marker_dict = {
                0.0: 'x',
                1.0: 'o',
                2.0: '^',
                3.0: 's',
                4.0: '.'
            }

            label_dict = {
                0.0: 'Ba/Sq',
                1.0: 'Stroma-rich',
                2.0: 'LumP',
                3.0: 'LumU',
                4.0: 'LumNS'
            }
            self.get_subset(subtype_col, subtype_names)
            rna = self.rna.dropna(axis='columns')

            subtypes = rna.loc[f'{subtype_col}', :]
            rna.drop(f'{subtype_col}', axis='rows')

            gene_exp = np.asarray(rna.loc[f'{gene}', :]) / np.asarray(rna.loc[f'{gene}', :]).max()
            marker = [marker_dict[x] for x in subtypes]

            x = np.asarray(pca_transformed[:, 0])
            y = np.asarray(pca_transformed[:, 1])

            ax = plt.subplot(111)

            for i in range(len(x)):
                ax.plot(x[i], y[i], marker=marker[i], color=gray(gene_exp[i]))
            plt.xlabel('component 1')
            plt.ylabel('component 2')

        handle = [self.make_legend('Ba/Sq'),
                  self.make_legend('Stroma-rich'),
                  self.make_legend('LumP'),
                  self.make_legend('LumU'),
                  self.make_legend('LumNS')]

        plt.legend(loc='best', handles=handle)
        plt.show()

    def compare_genes(self, gene1, gene2):
        g1 = self.rna.loc[gene1, :]
        g2 = self.rna.loc[gene2, :]
        r, pval = pearsonr(g1, g2)
        print(f'pearson correlation coeff: {r}, p_value: {pval}')
        g1.hist()
        plt.show()


RNASeqExplore('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data/RNAseq',
              'UROMOL-GENE-VST-subtypeswitch',
              'UROMOL-HistoPath',).plot_pca(gene='ENSG00000205420',
                                            subtype_col='MIBC consensus class',
                                            subtype_names=['Ba/Sq',
                                                           'Stroma-rich',
                                                           'LumP',
                                                           'LumU',
                                                           'LumNS'])

RNASeqExplore('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data/RNAseq',
              'UROMOL-GENE-VST-subtypeswitch',
              'UROMOL-HistoPath',).compare_genes('ENSG00000100644', 'ENSG00000186081')
