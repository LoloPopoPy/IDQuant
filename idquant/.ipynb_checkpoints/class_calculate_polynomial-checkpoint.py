import data_processing
import pandas as pd
import numpy as np
import logging

#Initializing logger
logger2 = logging.getLogger(__name__)
logger2.setLevel(logging.INFO)

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


mydata = pd.read_csv(r'C:\Users\legregam\Documents\idquant-master\test_data\Test_Data.csv', sep=';')

test = data_processing.DataSet(mydata)
test.prep_data()


class CalculateQuadraticPolynomial():
    '''Class to create a dataframe with quadratic polynomial predictions of concentrations from MS experiments '''
    
    def __init__(self, cal_data, sample_data):
        
        self.cal_data = cal_data
        self.sample_data = sample_data
        
    def get_data(self, metabolite):
        '''Getting data from dataframes'''
        
        logger2.debug("Getting_data for {}\n".format(metabolite))
        
        #Filter by metabolite
        tmp_cal_df = self.cal_data[
            self.cal_data["compound"] == metabolite]
        tmp_sample_df = self.sample_data[
            self.sample_data["compound"] == metabolite]
        
        cal_dict, data_dict = {}, {}
        cal_dict["names"] = tmp_cal_df["source"]
        cal_dict["x"] = tmp_cal_df["calibration concentration"].to_numpy()
        cal_dict['y'] = tmp_cal_df["M0/Mn"].to_numpy()
        
        data_dict["sources"] = tmp_sample_df["source"].to_list()
        data_dict["y_to_pred"] = tmp_sample_df["M0/Mn"].to_list()
        data_dict['calculated concentrations'] = []
        
        #We get the lowest and highest M0/Mn for determining if values to predict are in range
        self.x_limits = [min(cal_dict["y"]),
                    max(cal_dict["y"])]
        
        logger2.debug("checking limits for {}\n".format(metabolite))
        
        for ind, val in enumerate(data_dict["y_to_pred"]):
            if val < self.x_limits[0]:
                
                data_dict["y_to_pred"][ind] = "Under range"
                
                logger2.info("For {} the value {} is under range\n".format(metabolite, val))
            
            elif val > self.x_limits[1]:
                
                data_dict["y_to_pred"][ind] = "Over range"
                
                logger2.info("For {} the value {} is over range\n".format(metabolite, val))
        
        return cal_dict, data_dict
    
    
    def get_equation(self, cal_dict, metabolite):
        '''Calculating value of coefficients of the quadratic polynomial'''
        
        logger2.debug("Performing calculations for equation\n")
        
        x = cal_dict["x"]
        y = cal_dict["y"]
        
        polyfit = np.polynomial.polynomial.Polynomial.fit(x, y, 2) #Better than polyfit for numerical stability reasons (check docu for details)
        convert = polyfit.convert().coef
        equation = np.poly1d(convert)
        
        logger2.info("Equation for {} is equal to: \n{}\n ".format(metabolite, equation))
        
        #We get relative residuals, will be used later for residual plot
        cal_dict["relative_residuals"] = (equation(cal_dict["x"])-cal_dict["y"])/cal_dict["y"] 
        
        #Transform to list for deletions
        cal_dict["names"] = list(cal_dict["names"])
        cal_dict["x"] = list(cal_dict["x"])
        cal_dict["y"] = list(cal_dict["y"])
        cal_dict["relative_residuals"] = list(cal_dict["relative_residuals"])
        
        #Remove data points if residuals > or < by 20%
        for ind, val in enumerate(zip(cal_dict["names"], cal_dict["x"], cal_dict["y"], cal_dict["relative_residuals"])):
            if cal_dict["relative_residuals"][ind] < -20 or cal_dict["relative_residuals"][ind] > 20:
                del cal_dict["names"][ind]
                del cal_dict["x"][ind]
                del cal_dict["y"][ind]
                del cal_dict["relative_residuals"][ind]
                
                logger2.info("Removed data point x = {} from {} because residual too large\n".format(
                    cal_dict["x"][ind], metabolite))
        
        return cal_dict, equation
    
    def predict_x_value(self, equation, y_value):
        
        #To get roots from polynomial we must have equation in form ax² + bx + c - y = 0
        nul_eq = equation - y_value 
        x_pred = nul_eq.roots[1] #Remember we are solving for x values, not y
        
        return x_pred
    
    def calculate_predictions(self, data_dict, equation, metabolite):
        
        logger2.debug("Predicting concentrations\n")
        
        for val in data_dict["y_to_pred"]:
            if isinstance(val, str) == False: #Because if True that means that value is out of range
                data_dict["calculated concentrations"].append(self.predict_x_value(equation, val))
            elif val == "Under range":
                data_dict["calculated concentrations"].append("Under range") #See what I mean?
            elif val == "Over range":
                data_dict["calculated concentrations"].append("Over range") #See what I mean?
        
        return data_dict
        
    def rebuild_dataframes(self, metabolite, data_dict, cal_dict):
        
        logger2.debug("Rebuilding dataframes\n")
        
        names, concentrations = data_dict["sources"], data_dict["calculated concentrations"]
        
        data_df, cal_df = pd.DataFrame.from_dict({"names":names,"concentrations (µM)":concentrations}), pd.DataFrame.from_dict(cal_dict)
        
        data_df["metabolite"], cal_df["metabolite"] = metabolite, metabolite

        data_df.set_index(["metabolite", "names"], inplace=True)
        cal_df.set_index(["metabolite", "names"], inplace=True)
        
        return data_df, cal_df
        
    def processing(self):
        
        logger2.info("Starting to process calculations...\n")
        
        list_of_cal_dfs, list_of_data_dfs = [], []
        
        for metabolite in self.cal_data["compound"].unique():
            
            cal_dict, data_dict = self.get_data(metabolite)
            
            try:
                
                cal_dict, equation = self.get_equation(cal_dict, metabolite)
                
                data_dict = self.calculate_predictions(data_dict, equation, metabolite)
                
                data_df, cal_df = self.rebuild_dataframes(metabolite, data_dict, cal_dict)
                
                list_of_cal_dfs.append(cal_df)
                list_of_data_dfs.append(data_df)
            
            
            except Exception as err:
                
                logger2.error("There was a problem calculating for {}\n ".format(metabolite))
                logger2.error(err)
                
                continue
                
        self.final_data_df, self.final_cal_df = pd.concat(list_of_data_dfs), pd.concat(list_of_cal_dfs)
        
        self.final_data_df.to_excel(r'Calculated datas.xlsx', index=True)
        self.final_cal_df.to_excel(r'Calibration datas.xlsx', index=True)
        
        logging.info("Done!")
        
        
        
test2 = CalculateQuadraticPolynomial(test.ready_cal_df, test.ready_sample_df)
test2.processing()

