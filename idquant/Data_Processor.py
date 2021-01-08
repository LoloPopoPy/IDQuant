import pandas as pd
import logging


class Processor:
    """
    Class to prepare IDQ data before passing it to Polynomial_Calculator

    :param data: dataframe containing MS data (resulting from a pandas.read_excel or .read_csv)
    :type data: class: 'Pandas.Dataframe'
    :param conc_calib_df: dataframe containing calibration table
    :type conc_calib_df: class: 'Pandas.Dataframe'
    """

    def __init__(self, data, conc_calib_df):

        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        if not self.logger.handlers:
            handle = logging.StreamHandler()
            formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handle.setFormatter(formatter)
            self.logger.addHandler(handle)

        self.logger.debug('Initializing Data Object')

        self.data = data
        self.conc_calib_df = conc_calib_df
        self.data = self.data[["compound", "mz", "mz0", "mi", "area", "source"]]
        self.data.loc[:, 'compound'] = self.data.loc[:, 'compound'].str.replace('/', '_')
        self.conc_calib_df.loc[:, 'Compound'] = self.conc_calib_df.loc[:, 'Compound'].str.replace('/', '_')
        self.conc_calib_df.set_index("Compound", inplace=True)

        self.missing_values_df = None
        self.ready_cal_df = None
        self.ready_sample_df = None

    # Fonction pour isoler les M0
    @staticmethod
    def isolate_sample_M0(sample_df):
        """
        Get all the M0 from the sample data

        :param sample_df: Dataframe with MS data
        :type sample_df: class: 'Pandas.Dataframe'

        :return: dataframe containing M0 from each metabolite and sample
        :rtype: class: 'Pandas.Dataframe'
        """

        M0_df = sample_df[sample_df["mi"] == 0]
        M0_df.reset_index(drop=True, inplace=True)

        return M0_df

    # Fonction pour isoler les Mn
    @staticmethod
    def isolate_sample_Mn(sample_df):
        """
        Isolate highest isotopologue for each metabolite in sample data

        :param sample_df: dataframe with MS data
        :type sample_df: class: 'Pandas.Dataframe'

        :return: dataframe containing Mn for each sample and metabolite
        :rtype: class: 'Pandas.Dataframe'
        """

        list_of_Mn = []  # Container for selected rows
        samples = list(sample_df["source"].unique())
        compounds = list(sample_df["compound"].unique())

        for sample in samples:
            for compound in compounds:
                mi_no = len(sample_df[(sample_df["source"] == sample) & (sample_df["compound"] == compound)]["mi"]) - 1
                tmpdf = sample_df[
                    (sample_df["source"] == sample) & (sample_df["compound"] == compound) & (sample_df["mi"] == mi_no)]
                list_of_Mn.append(tmpdf)

        sample_Mn_df = pd.concat(list_of_Mn, ignore_index=True)

        return sample_Mn_df

    # Fonction pour séparer les points de calibration des points des échantillons
    @staticmethod
    def split_sample_cal(df_to_split):
        """
        Separate calibration points and sample points in MS measured data

        :param df_to_split: dataframe to split
        :type df_to_split: class: 'Pandas.Dataframe'

        :return splitted_calib_df: dataframe containing MS recorded data (calibration points)
        :rtype splitted_calib_df: class: 'Pandas.Dataframe'
        :return splitted_sample_df: dataframe containing MS recorded data (sample points)
        :rtype splitted_sample_df: class: 'Pandas.Dataframe'
        """

        splitted_calib_df = df_to_split[df_to_split["source"].str.contains("Cal")].copy()
        splitted_sample_df = df_to_split[df_to_split["source"].str.contains("Cal") is False].copy()

        return splitted_calib_df, splitted_sample_df

    # Fonction pour récupérer le numéro du point de gamme
    @staticmethod
    def get_cal_num(splitted_calib_df):
        """
        Get number for each calibration point in sample data

        :param splitted_calib_df: dataframe containing MS measured data (calibration points)
        :type splitted_calib_df: class: 'Pandas.Dataframe'

        :return numbered_cal_df: dataframe containing calibration points and associated number from calibration table
        :rtype: class: 'Pandas.Dataframe'
        """

        numbered_cal_df = splitted_calib_df.copy()

        for ind, row in splitted_calib_df.iterrows():
            val = splitted_calib_df.loc[ind, "source"][3]
            numbered_cal_df.at[ind, "cal_num"] = val

        numbered_cal_df = numbered_cal_df.sort_values(by=["compound", "cal_num"]).reset_index(drop=True)

        return numbered_cal_df

    # Fonction pour récupérer les concentrations
    def get_cal_concentration(self, numbered_cal_df):
        """
        Get concentration of calibration point in calibration datafile

        :param numbered_cal_df: dataframe containing calibration points and associated number from calibration table
        :type numbered_cal_df: class: 'Pandas.Dataframe'

        :return final_cal_df: final clean calibration dataframe
        :rtype: class: 'Pandas.Dataframe'
        """

        final_cal_df = numbered_cal_df.copy(deep=True)

        self.missing_values_df = pd.DataFrame(columns=final_cal_df.columns)
        missing_indices = []
        removed_metabolites = set()

        for ind, row in numbered_cal_df.iterrows():
            metabolite = row["compound"]
            cal_num = numbered_cal_df.at[ind, "cal_num"]

            try:
                final_cal_df.at[ind, "calibration concentration"] = self.conc_calib_df.loc[metabolite, str(cal_num)]

            except KeyError:
                missing_indices.append(ind)
                self.missing_values_df = self.missing_values_df.append(row)
                removed_metabolites.add(metabolite)
                final_cal_df.drop(ind, inplace=True)

                continue

            except Exception as e:
                self.logger.error('There was an unexpected problem: {}'.format(e))

        self.logger.warning('Concentrations for: {} are missing from calibration file\n'
                            .format(removed_metabolites))

        return final_cal_df

    # Fonction pour calculer les ratios M0/Mn
    @staticmethod
    def calculate_ratio(final_cal):
        """
        Calculate ratios between M0 and Mn and store in separate column

        :param final_cal: final clean calibration dataframe
        :type final_cal: class: 'Pandas.Dataframe'

        :return calculated_df: final clean calibration dataframe with M0/Mn ratio
        :rtype: class: 'Pandas.Dataframe'
        """

        calculated_df = final_cal.copy(deep=True)
        calculated_df["M0/Mn"] = final_cal.apply(
            lambda row: row.M0_area / row.Mn_area, axis=1)

        return calculated_df

    def final_cleanup(self, ready_cal_df, ready_sample_df):
        """
        Cleanup of NaN values before passing data to calculator

        :param ready_cal_df: final clean calibration dataframe
        :type ready_cal_df: class: 'Pandas.Dataframe'
        :param ready_sample_df: final clean sample dataframe
        :type ready_sample_df: class: 'Pandas.Dataframe'

        :return finished_cal_df: finished calibration dataframe ready for calculator
        :rtype finished_cal_df: class: 'Pandas.Dataframe'
        :return finished_sample_df: finished sample dataframe ready for calculator
        :rtype finished_sample_df: class: 'Pandas.Dataframe'
        """

        self.logger.info('Final cleanup...')
        self.logger.debug('Working on cal data')

        null_df = ready_cal_df[ready_cal_df.isnull().values]
        null_df = null_df.append(ready_sample_df[ready_sample_df.isnull().values])

        self.missing_values_df = self.missing_values_df.append(null_df)

        finished_cal_df = ready_cal_df.dropna()
        finished_sample_df = ready_sample_df.dropna()

        removed_metabolites = {met for met in self.missing_values_df['compound']}

        self.logger.warning('List of metabolites with missing data: {}\n For more details, check Missing values.xlsx'
                            .format(removed_metabolites))

        return finished_cal_df, finished_sample_df

    # Fonction appelé par le bouton pour contrôler les autres fonction et faire le traitement de données
    def prep_data(self):
        """Data Processing main function"""

        self.logger.info("Preparing Data... ")

        self.logger.debug("Getting M0 and MN")
        M0_df = Processor.isolate_sample_M0(self.data)
        Mn_df = Processor.isolate_sample_Mn(self.data)

        self.logger.debug("Splitting dataframes into cal and sample")
        cal_Mn, sample_Mn = Processor.split_sample_cal(Mn_df)
        cal_M0, sample_M0 = Processor.split_sample_cal(M0_df)

        self.logger.debug('Getting calibration number')
        numbered_cal_Mn = Processor.get_cal_num(cal_Mn)
        numbered_cal_M0 = Processor.get_cal_num(cal_M0)

        self.logger.debug('Getting concentrations')
        final_cal_Mn = self.get_cal_concentration(numbered_cal_Mn)
        final_cal_M0 = self.get_cal_concentration(numbered_cal_M0)

        self.logger.debug("Cleaning up")
        sample_M0.rename(columns={'area': 'M0_area'}, inplace=True)
        sample_M0.drop(["mz", "mi", "mz0"], axis=1, inplace=True)
        sample_Mn.rename(columns={'area': 'Mn_area'}, inplace=True)
        sample_Mn.drop(["mz", "mi", "mz0"], axis=1, inplace=True)

        final_cal_M0.rename(columns={'area': 'M0_area'}, inplace=True)
        final_cal_M0.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)
        final_cal_Mn.rename(columns={'area': 'Mn_area'}, inplace=True)
        final_cal_Mn.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)

        self.logger.debug("Merging dataframes")
        merge_cal_df = pd.merge(final_cal_M0, final_cal_Mn, on=['compound', 'source', "calibration concentration"])
        merge_sample_df = pd.merge(sample_M0, sample_Mn, on=['compound', 'source'])

        self.logger.debug("Calculating ratios")
        ready_sample_df = Processor.calculate_ratio(merge_sample_df)
        ready_cal_df = Processor.calculate_ratio(merge_cal_df)

        self.ready_cal_df, self.ready_sample_df = self.final_cleanup(ready_cal_df, ready_sample_df)

        self.ready_sample_df.to_excel('Sample Data.xlsx')
        self.ready_cal_df.to_excel('Calibration Sample data.xlsx')
        self.missing_values_df.to_excel('Missing values.xlsx')

        self.logger.info("Passing data to calculator...")
