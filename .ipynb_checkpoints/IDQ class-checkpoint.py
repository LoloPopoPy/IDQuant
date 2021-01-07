import pandas as pd
import seaborn as sns
import logging


logger = logging.getLogger(__name__)

formatter = logging.Formatter(fmt='%(asctime)s - %(levelname)s - %(message)s')
handle = logging.StreamHandler()
handle.setFormatter(formatter)

logger.addHandler(handle)
handle.setLevel(logging.DEBUG)

logger.warning('Testing')

class DataSet():
    
    
    def __init__(self, data):
        logger.debug('Initializing Data Object')
        self.data = data
        self.conc_calib_df = pd.read_csv(r'C:\Users\legregam\Documents\idquant-master\test_data\calibration.csv') 
        print(self.data)
        
    #Fonction pour isoler les M0
    def isolate_sample_M0(self, sample_df):
        M0_df = sample_df[sample_df["mi"] == 0]
        M0_df.reset_index(drop=True, inplace=True)
        return M0_df
    
    
    #Fonction pour isoler les Mn
    def isolate_sample_Mn(self, sample_df):
        list_of_Mn = []
        samples = list(sample_df["source"].unique())
        compounds = list(sample_df["compound"].unique())
        for sample in samples:
            for compound in compounds:
                mi_no = len(sample_df[(sample_df["source"] == sample) & (sample_df["compound"] == compound)]["mi"])-1
                tmpdf = sample_df[(sample_df["source"] == sample) & (sample_df["compound"] == compound)& (sample_df["mi"] == mi_no)]
                list_of_Mn.append(tmpdf)
        sample_Mn_df = pd.concat(list_of_Mn, ignore_index=True)
        return sample_Mn_df
    
    
    #Fonction pour séparer les points de calibration des points des échantillons 
    def split_sample_cal(self, df_to_split):
        splitted_calib_df = df_to_split[df_to_split["source"].str.contains("Cal")].copy()
        splitted_sample_df = df_to_split[df_to_split["source"].str.contains("Cal") == False].copy()
        return splitted_calib_df, splitted_sample_df
    
    
    #Fonction pour récupérer le numéro du point de gamme
    def get_cal_num(self, splitted_calib_df):
        numbered_cal_df = splitted_calib_df.copy()
        for ind, row in splitted_calib_df.iterrows():
            val = splitted_calib_df.loc[ind, "source"][3]
            numbered_cal_df.at[ind, "cal_num"] = val
        numbered_cal_df = numbered_cal_df.sort_values(by=["compound", "cal_num"]).reset_index(drop=True)
        return numbered_cal_df
    
    #Fonction pour récupérer les concentrations 
    def get_cal_concentration(self, numbered_cal_df, conc_calib_df):
        final_cal_df = numbered_cal_df.copy(deep=True)
        for ind, row in numbered_cal_df.iterrows():
            metabolite = row["compound"]
            cal_num = numbered_cal_df.at[ind, "cal_num"]
            final_cal_df.at[ind, "calibration concentration"] = conc_calib_df.loc[metabolite, str(cal_num)]
    
        return final_cal_df
    
    
    #Fonction pour calculer les ratios M0/Mn
    def calculate_ratio(self, final_cal):
        calculated_df = final_cal.copy(deep=True)
        calculated_df["M0/Mn"] = final_cal.apply(lambda row: row.M0_area / row.Mn_area, axis = 1)
    
        return calculated_df
    
    
    
    #Fonction appelé par le bouton pour contrôler les autres fonction et faire le traitement de données
    def prep_data(self):
    
        logger.info("Processing Data... ")
        
        logger.debug("Getting M0 and MN")
        M0_df = self.isolate_sample_M0(self.data)
        Mn_df = self.isolate_sample_Mn(self.data)
        
        logger.debug("Splitting dataframes into cal and sample")
        cal_Mn, sample_Mn = self.split_sample_cal(Mn_df)
        cal_M0, sample_M0 = self.split_sample_cal(M0_df)
        
        logger.debug('Getting calibration number')
        numbered_cal_Mn = self.get_cal_num(cal_Mn)
        numbered_cal_M0 = self.get_cal_num(cal_M0)
        
        logger.debug('Getting concentrations')
        final_cal_Mn, final_cal_M0 = self.get_cal_concentration(numbered_cal_Mn, self.conc_calib_df), self.get_cal_concentration(numbered_cal_M0, self.conc_calib_df)
        
        logger.debug("Cleaning up")
        sample_M0.rename(columns={'area' : 'M0_area'}, inplace=True)
        sample_M0.drop(["mz", "mi", "mz0"], axis=1, inplace=True)
        sample_Mn.rename(columns={'area' : 'Mn_area'}, inplace=True)
        sample_Mn.drop(["mz", "mi", "mz0"], axis=1, inplace=True)
    
        cal_M0.rename(columns={'area' : 'M0_area'}, inplace=True)
        cal_M0.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)
        cal_Mn.rename(columns={'area' : 'Mn_area'}, inplace=True)
        cal_Mn.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)
        
        logger.debug("Merging dataframes")
        merge_cal_df = pd.merge(final_cal_M0, final_cal_Mn, on=['compound', 'source', "calibration concentration" ])
        merge_sample_df = pd.merge(sample_M0, sample_Mn, on=['compound', 'source'])
        
        logger.debug("Calculating ratios")
        self.ready_sample_df = self.calculate_ratio(merge_sample_df)
        self.ready_cal_df = self.calculate_ratio(merge_cal_df)
        
    
        logger.info("Done!")          
'''         
    def plot_cal_metabolite(df, metabolite):
        tmpdf = df[df["compound"] == metabolite]
    
        myplot = sns.regplot(data=tmpdf,
                             x="calibration concentration", 
                             y="M0/Mn", 
                             line_kws={"color" : 'red'},
                             order=2
                            )
        plt.title(metabolite)                    
        plt.show(myplot)
 '''  
            


mydata = pd.read_csv(r'C:\Users\legregam\Documents\idquant-master\test_data\Test_Data.csv', sep=';')
test = DataSet(mydata)
#test.prep_data()
#test.ready_sample_df






