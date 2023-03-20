import pandas as pd
from collections import defaultdict


class Summary:

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

    def __init__(self, datadir, filename):
        # Initializing preprocessing class object
        # automatically generates a pandas dataframe
        self.datadir = datadir
        self.filename = filename
        self.df = pd.read_csv(f'{datadir}/{filename}.csv', index_col=0)

        staindata = defaultdict(dict)
        annotatedCellData = self.df[self.df['Analysis Region'] == 'Layer 1'].values.tolist()
        for celldata in annotatedCellData:
            spotid = celldata[10]
            coreid = celldata[6]
            totalcells = celldata[16]
            Gata0 = celldata[18]
            Gata1 = celldata[19]
            Gata2 = celldata[20]
            KRT0 = celldata[22]
            KRT1 = celldata[23]
            KRT2 = celldata[24]
            Both = celldata[25]
            Neither = celldata[26]

            staindata[f'{spotid}'][f'{coreid}'] = {'totalcells': totalcells,
                                                   'gata3': [Gata0, Gata1, Gata2],
                                                   'ck5': [KRT0, KRT1, KRT2],
                                                   'Dual Positive': Both,
                                                   'Dual Negative': Neither}
        self.celldata = dict(staindata)


Summary('/Users/andrewgarven/PycharmProjects/SubtypePlasticity/Data', 'TMA1_GATA3CK5_summary')

