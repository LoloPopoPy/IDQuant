import io
import logging
import os
import datetime

import ipywidgets as widgets
import pandas as pd

from data_processing import DataSet
from class_calculate_polynomial_v2 import Calculator



#Class pour récupérer la valeur en sortie de fonction appelée par un bouton
class ValueHolder():
    x: int = None
vh = ValueHolder()


vh.out = widgets.Output()
vh.debug_view = widgets.Output()

def make_upload_btn():
    
    global upload_data_btn
    
    upload_data_btn = widgets.FileUpload(
        accept='',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
        multiple=False, # True to accept multiple files upload else False
        description="Upload data"
    )
    
    return upload_data_btn
    
calib_layout = widgets.Layout(width='auto', height='auto')
    
def make_upload_conc_calib_btn():
    
    global upload_conc_calib_btn
    
    upload_conc_calib_btn = widgets.FileUpload(
        accept='',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
        multiple=False,  # True to accept multiple files upload else False
        description = 'Calib Data Upload',
        layout = calib_layout
    )

    return upload_conc_calib_btn

def make_text_box():
    
    global text_box
    
    text_box = widgets.Text(value='',
                           description='Run name:')
    return text_box

def make_log_lvl():
    
    global log_lvl
    
    log_lvl = widgets.Checkbox(value = False,
                              description = 'Debug mode'
                              )
    
    return log_lvl

submit_btn = widgets.Button(description='Submit data',
                            disabled=False,
                            button_style='',
                            tooltip='Click me',
                            icon='')

processing_btn = widgets.Button(description='Calculate',
                                disabled=False,
                                button_style='',
                                tooltip='Test',
                                icon='' )   


#Evenement pour soumettre les données MS et insérer dans un DF
@vh.debug_view.capture(clear_output=True)
def upload_data(event):
    
    #Récupération des datas mis en ligne par le bouton upload
    data_upload_filename = next(iter(upload_data_btn.value))
    data_content = upload_data_btn.value[data_upload_filename]['content']
    with open('myfile', 'wb') as f: f.write(data_content)
        
    #Entrons les datas dans un dataframe
    try:
        sample_df = pd.read_csv(io.BytesIO(data_content), sep=";")
    except:
        try:
            sample_df = pd.read_excel(io.BytesIO(data_content))
        except Exception as e:
            print('There was a problem reading file')
            print(e)
            
    #Mettons le df dans un ValueHolder pour l'avoir après l'utilisation du bouton
    vh.sample_df = sample_df
    
    with vh.out:
        print("MS Data has been loaded")
    
    return vh.sample_df


#Evenement pour soumettre la gamme de concentration et insérer dans un DF
def upload_conc_calib_data(event):
    
    #Récupération des datas mis en ligne par le bouton upload
    conc_calib_upload_filename = next(iter(upload_conc_calib_btn.value))
    conc_calib_content = upload_conc_calib_btn.value[conc_calib_upload_filename]['content']
    with open('myfile', 'wb') as f: f.write(conc_calib_content)
    
    #Entrons les datas dans un dataframe
    conc_calib_df = pd.read_csv(io.BytesIO(conc_calib_content), sep=";")
    
    #Mettons le df dans un ValueHolder pour l'avoir après l'utilisation du bouton
    vh.conc_calib_df = conc_calib_df
    
    with vh.out:
        print("Calibration data has been loaded")
    
    return vh.conc_calib_df

def main_process(event):
    
    now = datetime.datetime.now()
    date_time = now.strftime("%d%m%Y_%H%M%S") #Récupération date et heure
    mydir = os.getcwd()
    os.mkdir(date_time) #Créons le dir
    os.chdir(date_time) #Rentrons dans le dir
    
    if log_lvl.value:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
        
    data = DataSet(vh.sample_df, vh.conc_calib_df)
    data.prep_data()
    calculator = Calculator(text_box.value, log_level, data.ready_cal_df, data.ready_sample_df)
    calculator.main()
    
    os.chdir(mydir) #Revenir au dir initial   

submit_btn.on_click(upload_data)
submit_btn.on_click(upload_conc_calib_data)

with vh.out:
    processing_btn.on_click(main_process) 
