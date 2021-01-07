try:
    import pandas as pd
    import numpy as np
    import ipywidgets as widgets
    from ipywidgets import Layout
    import IPython
    import io
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import linregress
    from Metatoul_IDQuant import *
except ModuleNotFoundError as err:
    raise ModuleNotFoundError('Some dependencies might be missing. Check installation and try again')
except Exception:
    print('Unexpected Error')


#Création d'un output pour attraper les erreurs sur click de bouton
debug_view = widgets.Output(layout={'border': '1px solid black'})

#Class pour récupérer la valeur en sortie de fonction appelée par un bouton
class ValueHolder():
    x: int = None
vh = ValueHolder()


#Evenement pour soumettre les données MS et insérer dans un DF
@debug_view.capture(clear_output=True)
def upload_data(event):
    
    #Récupération des datas mis en ligne par le bouton upload
    data_upload_filename = next(iter(upload_data_btn.value))
    data_content = upload_data_btn.value[data_upload_filename]['content']
    with open('myfile', 'wb') as f: f.write(data_content)
        
    #Entrons les datas dans un dataframe
    sample_df = pd.read_csv(io.BytesIO(data_content), sep=";")

    #Sélectionnons les colonnes utiles pour nous
    sample_df = sample_df[["compound", "mz", "mz0", "mi", "area", "source"]]
    
    #Mettons le df dans un ValueHolder pour l'avoir après l'utilisation du bouton
    vh.sample_df = sample_df
    
    with out2:
        print("MS Data has been loaded")
    
    return vh.sample_df


#Evenement pour soumettre la gamme de concentration et insérer dans un DF
def upload_conc_calib_data(event):
    
    #Récupération des datas mis en ligne par le bouton upload
    conc_calib_upload_filename = next(iter(upload_conc_calib_btn.value))
    conc_calib_content = upload_conc_calib_btn.value[conc_calib_upload_filename]['content']
    with open('myfile', 'wb') as f: f.write(conc_calib_content)
    
    #Entrons les datas dans un dataframe
    conc_calib_df = pd.read_csv(io.BytesIO(conc_calib_content), sep="\t")
    
    #Un peu de mise en forme
    conc_calib_df.set_index("Compound", drop=True, inplace=True)
    conc_calib_df.drop("RT", axis=1, inplace=True)
    conc_calib_df.drop("Unnamed: 10", axis=1, inplace=True)
    
    #Mettons le df dans un ValueHolder pour l'avoir après l'utilisation du bouton
    vh.conc_calib_df = conc_calib_df
    
    with out2:
        print("Calibration data has been loaded")
    
    return vh.conc_calib_df
    
    
#Fonction pour isoler les M0
def isolate_sample_M0(sample_df):
    M0_df = sample_df[sample_df["mi"] == 0]
    M0_df.reset_index(drop=True, inplace=True)
    return M0_df


#Fonction pour isoler les Mn
def isolate_sample_Mn(sample_df):
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
def split_sample_cal(df_to_split):
    splitted_calib_df = df_to_split[df_to_split["source"].str.contains("Cal")].copy()
    splitted_sample_df = df_to_split[df_to_split["source"].str.contains("Cal") == False].copy()
    return splitted_calib_df, splitted_sample_df


#Fonction pour récupérer le numéro du point de gamme
def get_cal_num(splitted_calib_df):
    numbered_cal_df = splitted_calib_df.copy()
    for ind, row in splitted_calib_df.iterrows():
        val = splitted_calib_df.loc[ind, "source"][3]
        numbered_cal_df.at[ind, "cal_num"] = val
    numbered_cal_df = numbered_cal_df.sort_values(by=["compound", "cal_num"]).reset_index(drop=True)
    return numbered_cal_df

#Fonction pour récupérer les concentrations 
def get_cal_concentration(numbered_cal_df, conc_calib_df):
    final_cal_df = numbered_cal_df.copy(deep=True)
    for ind, row in numbered_cal_df.iterrows():
        metabolite = row["compound"]
        cal_num = numbered_cal_df.at[ind, "cal_num"]
        final_cal_df.at[ind, "calibration concentration"] = vh.conc_calib_df.loc[metabolite, str(cal_num)]
        
    return final_cal_df


#Fonction pour calculer les ratios M0/Mn
def calculate_ratio(final_cal):
    calculated_df = final_cal.copy(deep=True)
    calculated_df["M0/Mn"] = final_cal.apply(lambda row: row.M0_area / row.Mn_area, axis = 1)
    
    return calculated_df

def plot_cal_metabolite(df, metabolite):
    tmpdf = df[df["compound"] == metabolite]

    myplot = sns.regplot(data=tmpdf,
                         x="calibration concentration", 
                         y="M0/Mn", 
                         line_kws={"color" : 'red'},
                         order=1
                        )
    plt.title(metabolite)                    
    plt.show(myplot)
    
def calculate_regression_values(df, metabolite):
    
    x = df['calibration concentration'].to_numpy()
    y = df['M0/Mn'].to_numpy()
        
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    
    return slope, intercept, r_value, std_err

def predict_concentration(intercept, slope, M0_Mn):
    
    concentration = (M0_Mn-intercept)/slope
    return concentration

#Fonction appelé par le bouton pour contrôler les autres fonction et faire le traitement de données
@debug_view.capture(clear_output=True)
def main(event):
    
    with out2:
        print("Processing Data... ")
    
    M0_df = isolate_sample_M0(vh.sample_df)
    Mn_df = isolate_sample_Mn(vh.sample_df)
    
    splitted_cal_Mn, splitted_sample_Mn = split_sample_cal(Mn_df)
    splitted_cal_M0, splitted_sample_M0 = split_sample_cal(M0_df)
    
    numbered_cal_Mn = get_cal_num(splitted_cal_Mn)
    numbered_cal_M0 = get_cal_num(splitted_cal_M0)
    
    final_cal_Mn, final_cal_M0 = get_cal_concentration(numbered_cal_Mn, vh.conc_calib_df), get_cal_concentration(numbered_cal_M0, vh.conc_calib_df)
    
    splitted_sample_M0.rename(columns={'area' : 'M0_area'}, inplace=True)
    splitted_sample_M0.drop(["mz", "mi", "mz0"], axis=1, inplace=True)
    splitted_sample_Mn.rename(columns={'area' : 'Mn_area'}, inplace=True)
    splitted_sample_Mn.drop(["mz", "mi", "mz0"], axis=1, inplace=True)
    
    final_cal_M0.rename(columns={'area' : 'M0_area'}, inplace=True)
    final_cal_M0.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)
    final_cal_Mn.rename(columns={'area' : 'Mn_area'}, inplace=True)
    final_cal_Mn.drop(["mz", "mi", "mz0", "cal_num"], axis=1, inplace=True)
    
    merge_cal_df = pd.merge(final_cal_M0, final_cal_Mn, on=['compound', 'source', "calibration concentration" ])
    merge_sample_df = pd.merge(splitted_sample_M0, splitted_sample_Mn, on=['compound', 'source'])
    
    calculated_sample_df = calculate_ratio(merge_sample_df)
    calculated_cal_df = calculate_ratio(merge_cal_df)
    
    finished_sample_df_list = []
    
    for metabolite in calculated_cal_df['compound'].unique():
        
        tmpdf = calculated_sample_df[calculated_sample_df["compound"] == metabolite].copy()
        cal_tmpdf = calculated_cal_df[calculated_cal_df["compound"] == metabolite].copy()
        
        slope, intercept, r_value, SD = calculate_regression_values(cal_tmpdf, metabolite)
        with out2:
            print("slope: {}, intercept: {}, R²: {},  SD: {}".format(slope, intercept, r_value, SD) )
            
        with out2:
            plot_cal_metabolite(cal_tmpdf, metabolite)
            
        for ind, row in tmpdf.iterrows():
            tmpdf.at[ind, "Concentration calculée"] = predict_concentration(intercept, slope, tmpdf.at[ind, "M0/Mn"])
        
        finished_sample_df_list.append(tmpdf)
        
    vh.finished_sample_df = pd.concat(finished_sample_df_list)
    
    with out2:
        print("Done!")

#Création du dashboard
out = widgets.Output()
out2 = widgets.Output()


submit_btn.on_click(upload_data)
submit_btn.on_click(upload_conc_calib_data)
processing_btn.on_click(main)    
