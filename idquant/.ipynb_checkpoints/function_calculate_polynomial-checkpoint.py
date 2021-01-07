from data_processing import *
import pandas as pd
import numpy as np
from natsort import natsorted
import logging
import matplotlib.pyplot as plt
import seaborn as sns
import os
import datetime

        
now = datetime.datetime.now()
date_time = now.strftime("%d%m%Y_%H%M%S") #Récupération date et heure
mydir = os.getcwd()
os.mkdir(date_time) #Créons le dir
os.chdir(date_time) #Rentrons dans le dir
    

#Initializing logger
logger2 = logging.getLogger(__name__)
logger2.setLevel(logging.DEBUG)

if not logger2.handlers:
    info_handle = logging.FileHandler('IDQuant.log') #user log goes to file
    debug_handle = logging.StreamHandler() #Debug log is in console
    
    formatter = logging.Formatter(fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    info_handle.setFormatter(formatter)
    debug_handle.setFormatter(formatter)
    
    info_handle.setLevel(logging.INFO)
    debug_handle.setLevel(logging.DEBUG)
    
    logger2.addHandler(info_handle)
    logger2.addHandler(debug_handle)
        
def get_data(cal_data, sample_data, metabolite):
    '''Getting data from dataframes
    
    :param cal_dict: Dictionnary containing cal data for one met
    :type cal_dict: dict
    :param data_dict: Dict containing sample data with y to predict
    :type data_dict: dict
    :param y_limits: upper and lower limits from cal data 
    :type y_limits: list of floats
    '''
    
    logger2.info("Getting_data for {}\n".format(metabolite))
    
    #Filter by metabolite
    tmp_cal_df = cal_data[
        cal_data["compound"] == metabolite]
    tmp_sample_df = sample_data[
        sample_data["compound"] == metabolite]
    
    #One dictionnary for calibration datas
    cal_dict, data_dict = {}, {}
    cal_dict["names"] = tmp_cal_df["source"].to_list()
    cal_dict["x"] = tmp_cal_df["calibration concentration"].to_list()
    cal_dict['y'] = tmp_cal_df["M0/Mn"].to_list()
    
    #One dictionnary for sample data and predicted values
    data_dict["sources"] = tmp_sample_df["source"].to_list()
    data_dict["y_to_pred"] = tmp_sample_df["M0/Mn"].to_list()
    data_dict['calculated concentrations'] = []
    
    #We get the lowest and highest M0/Mn for determining if values to predict are in range
    y_limits = [min(cal_dict["y"]),
                max(cal_dict["y"])]
    
    logger2.debug("checking limits for {}\n".format(metabolite))
    
    for ind, val in enumerate(data_dict["y_to_pred"]):
        if val < y_limits[0]:
            
            data_dict["y_to_pred"][ind] = "Under range"
            
            logger2.info("For {} the value to predict ({}) is under range\n".format(metabolite, val))
        
        elif val > y_limits[1]:
            
            data_dict["y_to_pred"][ind] = "Over range"
            
            logger2.info("For {} the value to predict ({}) is over range\n".format(metabolite, val))
    
    return cal_dict, data_dict

def get_residuals(cal_dict, equation):
    '''Function to calculate residuals in calibration data from equation'''
    
    cal_dict["relative_residuals"] = []
    for x, y in zip(cal_dict["x"], cal_dict['y']):
        if x == 0 or y == 0:
            x= np.nan
            y = np.nan
            
        cal_dict["relative_residuals"].append((equation(x)-y)/y)
    
    return cal_dict

def equation_and_plot(cal_dict, metabolite):
    '''Calculating value of coefficients of the quadratic polynomial
    
    :param nans: NaN values that must be removed from cal_dict
    :type nans: list of integers
    :param polyfit: series instance that is the least squares fit to the data y sampled at x
    :type polyfit: classmethod: Polynomial.fit
    :param convert: coefficients for the unscaled and unshifted basis polynomials
    :type convert: list of coefficients
    :param equation: equation of the polynomial fit to calibration data for metabolite of interest
    :type equation: class: numpy.poly1D
    :param r2: coefficient of determination
    :type r2: float
    '''
    
    #Checking values and removing NaNs
    nans = [ind for ind, val in enumerate(cal_dict["y"]) if np.isnan(val)]
    
    logger2.debug('nans list: {}'.format(nans))
    logger2.debug('Cal dict x: {}\n cal dict y: {}\n cal dict names: {}\n'
                  .format(cal_dict['x'], cal_dict['y'], cal_dict['names']))
    
    nancount = 0 #Counter to substract from indices 
    
    for ind in nans:
        
        logger2.info("Removing {} from {} calibration datas because nan value detected"
                     .format(cal_dict["names"][ind - nancount], metabolite))
        
        del(cal_dict['x'][ind - nancount])
        del(cal_dict['y'][ind - nancount])
        del(cal_dict['names'][ind - nancount])
        nancount += 1
        
    if len(cal_dict['names']) < 4:
        logger2.error("Couldn't do calibration for {} because not enough calibration points (less than 4)"
                     .format(metabolite))
    
    else:
        logger2.debug('Cal dict x: {}\n cal dict y: {}\n cal dict names: {}\n'
                  .format(cal_dict['x'], cal_dict['y'], cal_dict['names']))
    
        #Transform lists to arrays for calculations
        cal_dict["x"] = np.array(cal_dict["x"])
        cal_dict['y'] = np.array(cal_dict['y'])

        #Equation part
        x = cal_dict["x"]
        y = cal_dict["y"]            

        logger2.debug("Performing calculations for equation\n x={} \n y={}".format(x, y))

        #Better than polyfit for numerical stability reasons (check docu for details)
        try:
            polyfit = np.polynomial.polynomial.Polynomial.fit(x, y, 2) 
            convert = polyfit.convert().coef
            equation = np.poly1d(convert)
            r2 = round(np.corrcoef(cal_dict['x'], cal_dict['y'])[0,1]**2, 3)

        except Exception as e:
            logger2.error('Error during calculation of equation for {}. Error: {}\n'
                         .format(metabolite, e))

        logger2.info("Equation for {} is equal to: \n{}\n ".format(metabolite, equation))

        logger2.debug("Calculating residuals \n")

        #We get relative residuals, will be used later for residual plot
        cal_dict = get_residuals(cal_dict, equation)

        #Transform to list for deletions
        cal_dict["x"] = list(cal_dict["x"])
        cal_dict["y"] = list(cal_dict["y"])
        cal_dict["relative_residuals"] = list(cal_dict["relative_residuals"])

        #Plotting part
        ipass = 0 #Counter for removing residuals
        npass = 1 #Plot number

        #Create plots a first time
        plot_reg(cal_dict, metabolite, r2, npass)

        plot_res(cal_dict, metabolite, r2, npass)

        #Remove data points if residuals > or < by 20%    
        for ind, val in enumerate(cal_dict["relative_residuals"]):
            if val < -20 or val > 20:
                ipass+=1

        #Keep plotting if some residuals are still too high or low
        while ipass > 0:
            logger2.info("Removed data point {} from {} because residual too large\n".format(
                cal_dict["names"][0], metabolite))  

            del(cal_dict["relative_residuals"][0])
            del(cal_dict["x"][0])
            del(cal_dict["y"][0])
            del(cal_dict["names"][0])

            npass +=1
            plot_reg(cal_dict, metabolite, r2, npass)
            plot_res(cal_dict, metabolite, r2, npass)
            ipass-=1

        return cal_dict, equation, r2

def predict_x_value(equation, y_value):
    '''Function to calculate concentration from value using the roots of the polynomial'''
    
    #To get roots from polynomial we must have equation in form ax² + bx + c - y = 0
    logger2.debug('equation is equal to: {}\n y_value is equal to: {}\n'.format(equation, y_value))
    
    if not np.isfinite(y_value):
        logger2.error('Error: exp. value ({}) must be a number.'.format(y))
    else:
        nul_eq = equation - y_value 
        logger2.debug('roots are equal to: {}\n'.format(nul_eq.roots))
        x_pred = nul_eq.roots[1] #Remember we are solving for x values, not y
        x_pred = x_pred[np.isreal(x_pred)]
    
    return x_pred

def calculate_predictions(data_dict, equation, metabolite):
    '''Function to coordinate predictions and check that datas are conform'''
    
    logger2.debug("Predicting concentrations\n")
    
    for val in data_dict["y_to_pred"]:
        logger2.debug("Trying to predict {} of type {} \n".format(val, type(val)))
        
        if isinstance(val, str):
            
            if val == "Under range":
                data_dict["calculated concentrations"].append("Under range")
            elif val == "Over range":
                data_dict["calculated concentrations"].append("Over range") 
            else:
                logger2.warning('Recieved non admissable value to predict: {} \n'.format(val))
                continue
            
        elif np.isnan(val):
            logger2.warning('NaN detected in values to predict')
            data_dict["calculated concentrations"].append("NaN")
            
        else: 
            try:
                logger2.debug('The calculated value is equal to: {}\n'.format(predict_x_value(equation, val)))
                data_dict["calculated concentrations"].append(predict_x_value(equation, val))
            except Exception as e:
                logger2.error('There was a problem while predicting the value from "{}".\nError: {}'
                            .format(val, e))
    
    return data_dict
    
def rebuild_dataframes(metabolite, data_dict, cal_dict):
    '''Function to assemble datas after calculations into final dataframe'''
    
    logger2.debug("Rebuilding dataframes\n")
    
    names, concentrations = data_dict["sources"], data_dict["calculated concentrations"]
    
    data_df, cal_df = pd.DataFrame.from_dict({"names":names,"concentrations (µM)":concentrations}), pd.DataFrame.from_dict(cal_dict)
    
    data_df = data_df.assign(metabolite = metabolite)
    cal_df = cal_df.assign(metabolite = metabolite)
    
    data_df.set_index(['metabolite', 'names'], inplace=True)
    cal_df.set_index(['metabolite', 'names'], inplace=True)
    
    data_df.index = pd.MultiIndex.from_tuples(natsorted(data_df.index), names=['metabolite', 'names'])
    cal_df.index = pd.MultiIndex.from_tuples(natsorted(cal_df.index), names=['metabolite', 'names'])
    
    return data_df, cal_df

def plot_reg(cal_dict, metabolite, r2, npass):
    '''Function to plot the polynomial regression curve'''

    regression_plot = sns.regplot(data=cal_dict,
                         x=cal_dict["x"], 
                         y=cal_dict["y"], 
                         ci = None,
                         line_kws={"color" : 'red'},
                         order=2
                        )

    plt.xlabel("Concentration (µM)")
    plt.ylabel("M0/Mn")
    plt.title('{} regression plot\nR²={}\nPass no={}'.format(metabolite, r2, npass))

    plt.savefig('{} regression plot {}.pdf'.format(metabolite, npass), format = 'pdf')     

    plt.show(regression_plot)

def plot_res(cal_dict, metabolite, r2, npass):
    
    residual_plot = plt.scatter(x=cal_dict['x'],
                                y=cal_dict['relative_residuals'],
                                c='b')
    
    #Get limits for plot y size depending on value of minimum and maximum relatiive residuals
    if min(cal_dict['relative_residuals']) < -30:
        min_ylim = min(cal_dict['relative_residuals']) * 1.25
    else:
        min_ylim = -30
    
    if max(cal_dict['relative_residuals']) > 30:
        max_ylim = max(cal_dict['relative_residuals']) * 1.25
    else:
        max_ylim = 30
        
    plt.ylim(min_ylim, max_ylim)        
    plt.hlines([-20, 20], -10, max(cal_dict['x'])+10, 'r')
    plt.hlines(0, -10, max(cal_dict['x'])+10, linestyles='dashed', colors='k')
    
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.title('{} residual plot\nR²={}\nPass no= {}\n'.format(metabolite, r2, npass))          

    plt.savefig('{} residual plot {}.pdf'.format(metabolite, npass), format = 'pdf',
               bbox_inches='tight')    

    plt.show(residual_plot)    

        
def main(cal_data, sample_data):
    '''
    Main function to create a dataframe with quadratic polynomial 
    predictions of concentrations from MS experiments 
    '''
    
    logger2.info("Starting to process calculations...\n")
    
    list_of_cal_dfs, list_of_data_dfs = [], []

    for metabolite in cal_data["compound"].unique():
        
        cal_dict, data_dict = get_data(cal_data, sample_data, metabolite)
        
        try:
            
            cal_dict, equation, r2 = equation_and_plot(cal_dict, metabolite)
            
            data_dict = calculate_predictions(data_dict, equation, metabolite)
            
            data_df, cal_df = rebuild_dataframes(metabolite, data_dict, cal_dict)
            
            list_of_cal_dfs.append(cal_df)
            list_of_data_dfs.append(data_df)
            

        
        except Exception as err:
            
            logger2.error("There was a problem calculating for {}\n Error: {} \n ".format(metabolite, err))
            
            continue
            
    final_data_df, final_cal_df = pd.concat(list_of_data_dfs), pd.concat(list_of_cal_dfs)
    
    final_data_df.to_excel(r'Calculated datas.xlsx', index=True)
    final_cal_df.to_excel(r'Calibration datas.xlsx', index=True)

    logger2.info("Done!")
    
    os.chdir(mydir) #Revenir au dir initial    
    
    return final_cal_df

        



