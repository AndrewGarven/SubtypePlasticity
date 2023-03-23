import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt


class HALOSummary:

    """HALOSummary is used to generate a summary data from a HALO image analysis software output .csv file

        Typical use:

            report : HALOSummary(datadir, filename, stain1, stain2, loc1, loc2).report()

        returns a summarized report including the following:
            1. histogram of the relative average optical density values on the Tissue Micro-Array (TMA)

        Alt use:
             plot_hist : HALOSummary(datadir, filename, stain1, stain2, loc1, loc2).plot_hist()

        :datadir: [str] of directory containing nuclear position data
        :filename: [str] of (.csv) file containing HALO image analysis output
        :stain1: [str] name of first stain (ex. 'Gata3')
        :stain2: [str] name of second stain (ex. 'Ck5')
        :loc1: [str] cellular staining location nuclear: 'nuclear' or 'n',
                                                cytoplasmic: 'cytoplasmic' or 'c'
        :loc2: [str] cellular staining location nuclear: 'nuclear' or 'n',
                                                cytoplasmic: 'cytoplasmic' or 'c'

        """

    def __init__(self, datadir, file_name, stain1='stain1', stain2='stain2', loc1='n', loc2='n', anno=True):

        # Initializing preprocessing class object
        # automatically generates a pandas dataframe
        self.datadir = datadir
        self.file_name = file_name
        self.df = pd.read_csv(f'{datadir}/{file_name}.csv', index_col=0)
        self.stain1 = stain1
        self.stain2 = stain2
        self.loc1 = loc1
        self.loc2 = loc2
        self.stain_data = defaultdict(dict)

        if anno:
            cell_data = self.df[self.df['Analysis Region'] == 'Layer 1'].values.tolist()
        else:
            cell_data = self.df[self.df['Analysis Region'] == 'entire image'].values.tolist()

        for tma in cell_data:
            spot_id = tma[10]
            core_id = tma[6]
            total_cells = tma[16]

            stain1_0 = tma[18]
            stain1_1 = tma[19]
            stain1_2 = tma[20]
            stain1_3 = tma[21]
            stain1_h_score = tma[34]

            stain2_0 = tma[23]
            stain2_1 = tma[24]
            stain2_2 = tma[25]
            stain2_3 = tma[26]
            stain2_h_score = tma[40]

            dual_positive = tma[27]
            dual_negative = tma[28]

            if loc1.lower() in ['cytoplasmic', 'c']:
                stain1_optical_density = tma[45]
            else:
                stain1_optical_density = tma[43]

            if loc2.lower() in ['cytoplasmic', 'c']:
                stain2_optical_density = tma[46]
            else:
                stain2_optical_density = tma[44]

            self.stain_data[f'{spot_id}'][f'{core_id}'] = {'total_cells': total_cells,
                                                           self.stain1: {'score_count': [stain1_0,
                                                                                         stain1_1,
                                                                                         stain1_2,
                                                                                         stain1_3],
                                                                         'h_score': stain1_h_score,
                                                                         'optical_density': stain1_optical_density},
                                                           self.stain2: {'score_count': [stain2_0,
                                                                                         stain2_1,
                                                                                         stain2_2,
                                                                                         stain2_3],
                                                                         'h_score': stain2_h_score,
                                                                         'optical_density': stain2_optical_density},
                                                           'dual_positive': dual_positive,
                                                           'dual_negative': dual_negative}

        self.stain_data = dict(self.stain_data)
        print(self.stain_data['211'])

    def plot_hist(self, save_image=False, image_dir=None):

        """
            plot_hist() generates a matplotlib.pyplot.hist() plot from HALO optical density values

            :self: HALOSummary Object described above
            :save_image: set True if you would like a png image of histogram
            :image_path: [str] contain path to directory to save png image file

            :return: 1. (default) matplotlib.pyplot.hist() plot
                     2. (if image = True): PNG image file of histogram
                        and/or size of node image saved to (default): 'datadir'/'filename'.png
                                                           (if image_path specified): :image_path:/'filename'.png
        """

        stain_1 = []
        stain_2 = []
        plot_title = f'Optical Density Histogram [{self.stain1} and {self.stain2}]'
        for patient in list(self.stain_data.keys()):
            patient_stain1 = []
            patient_stain2 = []
            for tma in list(self.stain_data[patient].keys()):
                od1 = self.stain_data[patient][tma][self.stain1]['optical_density']
                patient_stain1.append(od1)
                od2 = self.stain_data[patient][tma][self.stain2]['optical_density']
                patient_stain2.append(od2)
            stain_1.extend(patient_stain1)
            stain_2.extend(patient_stain2)

        plt.hist(stain_1, alpha=0.5, label=self.stain1)
        plt.hist(stain_2, alpha=0.5, label=self.stain2)
        plt.title(plot_title)
        plt.xlabel('Optical Density Value')
        plt.ylabel('Number of Patients')
        plt.legend(loc='upper right')

        if save_image:
            if image_dir:
                plt.savefig(f'{image_dir}/{self.file_name}.png')
            else:
                plt.savefig(f'{self.datadir}/{self.file_name}.png')
        else:
            plt.show()


HALOSummary('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data',
            'TMA1_GATA3CK5_summary',
            stain1='Gata3',
            stain2='Ck5',
            loc1='n',
            loc2='c')
