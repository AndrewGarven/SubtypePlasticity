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
        :clinical: [str] of (.csv) file name [containing clinical data (including Subtype, Grade, etc)]
                  - rows: Patients
                  - columns: Clinical Features (ex. Subtype, Grade, etc)
        :stain1: [str] name of first stain (ex. 'Gata3')
        :stain2: [str] name of second stain (ex. 'Ck5')
        :loc1: [str] cellular staining location nuclear: 'nuclear' or 'n',
                                                cytoplasmic: 'cytoplasmic' or 'c'
        :loc2: [str] cellular staining location nuclear: 'nuclear' or 'n',
                                                cytoplasmic: 'cytoplasmic' or 'c'

        """

    def __init__(self, datadir, file_name, clinical, stain1='stain1', stain2='stain2', loc1='n', loc2='n', anno=True):

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

            stain1_1 = tma[19]
            stain1_2 = tma[20]
            stain1_3 = tma[21]
            stain1_h_score = tma[34]

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
                                                           self.stain1: {'score_count': [stain1_1,
                                                                                         stain1_2,
                                                                                         stain1_3],
                                                                         'h_score': stain1_h_score,
                                                                         'optical_density': stain1_optical_density},
                                                           self.stain2: {'score_count': [stain2_1,
                                                                                         stain2_2,
                                                                                         stain2_3],
                                                                         'h_score': stain2_h_score,
                                                                         'optical_density': stain2_optical_density},
                                                           'dual_positive': dual_positive,
                                                           'dual_negative': dual_negative}

        self.stain_data = dict(self.stain_data)
        self.clinical = pd.read_csv(f'{datadir}/{clinical}.csv', index_col='Sample_ID')

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
                od1 = self.stain_data[int(patient)][tma][self.stain1]['optical_density']
                patient_stain1.append(od1)
                od2 = self.stain_data[int(patient)][tma][self.stain2]['optical_density']
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

    def get_data(self, attribute: str, group_names: list):

        """
            get_data() retrieves cell type (stratified by group) data from HALO tma summary dataframe

            :self: HALOSummary Object described above
            :attribute: [str] column name within clinical dataframe that you want to subset (ex. 'tumor_recurrence')
            :group_names: [list] of subgroups within :attribute: (ex. ['yes', 'no', 'maybe'])

            :return: [dict] containing cell type data stratified by group
        """

        int_id = list(dict.fromkeys(self.df['Spot Id 1'].values.tolist()))
        strid = [int(x) for x in int_id if str(x) != 'nan' and int(x) not in [868, 17, 220, 314, 8, 666]]
        property_df = self.clinical.loc[strid, attribute]

        group_ids = {}
        cell_data = {}

        for i, group in enumerate(group_names):
            stratified_patient_id = property_df[property_df == group]
            group_ids[group] = stratified_patient_id

        for group in group_ids.keys():
            for patient_id in group_ids[group].index.values.tolist():
                if str(float(patient_id)) in self.stain_data.keys():
                    for i, stain in enumerate(list(self.stain_data[str(float(patient_id))].keys())):

                        stain_dict = self.stain_data[str(float(patient_id))][stain]
                        if stain_dict['total_cells'] == 0:
                            continue
                        else:
                            norm_gata3 = [int(x)/int(stain_dict['total_cells']) for x in stain_dict['Gata3']['score_count']]
                            norm_ck5 = [int(x)/int(stain_dict['total_cells']) for x in stain_dict['Ck5']['score_count']]
                            norm_dual_neg = [int(stain_dict['dual_negative'])/int(stain_dict['total_cells'])]
                            norm_dual_pos = [int(stain_dict['dual_positive'])/int(stain_dict['total_cells'])]

                            cell_data[str(patient_id) + '-' + stain] = norm_gata3+norm_ck5+norm_dual_neg+norm_dual_pos
        return cell_data

    def plot_compare(self, attribute: str, group_names: list):
        """
            plot_compare() generates a matplotlib.pyplot.bar() graph object from HALO TMA summary data for each group
            within 'group_names'

            :self: HALOSummary Object described above
            :attribute: [str] column name within clinical dataframe that you want to subset (ex. 'tumor_recurrence')
            :group_names: [list] of subgroups within :attribute: (ex. ['yes', 'no', 'maybe'])

            :return: (display) matplotlib.pyplot.bar() graph for each group
        """

        cell_data = self.get_data(attribute=attribute, group_names=group_names)
        cell_data_dict = {}

        index_labels = ['Gata3-1', 'Gata3-2', 'Gata3-3', 'Ck5-1', 'Ck5-2', 'Ck5-3', 'dual negative', 'dual positive']
        color_spectrum = ['limegreen', 'green', 'darkgreen', 'blue', 'mediumblue', 'darkblue', 'white', 'black']
        x = [1, 2, 3, 4, 5, 6, 7, 8]
        for group in cell_data.keys():
            cell_data_dict[group] = pd.DataFrame(cell_data[group], index=index_labels)
            plt.bar(x=x,
                    height=(cell_data_dict[group].mean(axis=1) * 100),
                    color=color_spectrum)
            plt.title(f'Cell Type Breakdown Stratified by {attribute}')
            plt.xticks(x, labels=index_labels)
            plt.xlabel('Cell Type')
            plt.ylabel('Percentage Cell Count')
            plt.show()

        return cell_data_dict


low1 = HALOSummary('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data/PositionInputData',
                   'TMA1_GATA3CK5_summary',
                   'NMIBC_clinical',
                   stain1='Gata3',
                   stain2='Ck5',
                   loc1='n',
                   loc2='c').plot_compare('Grade_1.x', [1, 2])


