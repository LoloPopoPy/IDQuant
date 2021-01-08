import logging
import copy

import pandas as pd
import numpy as np
from natsort import index_natsorted
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns


class Calculator:
    """
    IDQ calculator for building regressions and plots

    :param run_name: Name of the run
    :type run_name: str
    :param log_level: Level to set the logger to
    :type log_level: int
    :param cal_data: Data from calibration table
    :type cal_data: class: 'pandas.Dataframe'
    :param sample_data: Data from MS Experiments
    :type sample_data: class: 'pandas.Dataframe'

    """

    def __init__(self, run_name, log_level, cal_data, sample_data):

        # Initialize object parameters
        self.run_name = run_name
        self.log_level = log_level
        self.cal_data = cal_data
        self.sample_data = sample_data
        self.equations = {}
        self.rsquares = {}
        self.linear_convert = None
        self.linear_equation = None
        self.linear_r2 = None
        self.quad_convert = None
        self.quad_equation = None
        self.quad_r2 = None

        self.tmp_cal_dict = {}
        self.tmp_sample_dict = {}
        self.final_data_df = None
        self.final_cal_df = None

        # Initializing logger
        self._logger = logging.getLogger(self.run_name)
        self._logger.setLevel(self.log_level)

        self._file_handle = logging.FileHandler(run_name + '.log')
        self._stream_handle = logging.StreamHandler()

        self._log_formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        self._file_handle.setFormatter(self._log_formatter)
        self._stream_handle.setFormatter(self._log_formatter)

        self._file_handle.setLevel(logging.DEBUG)
        self._stream_handle.setLevel(logging.DEBUG)

        self._logger.addHandler(self._file_handle)
        self._logger.addHandler(self._stream_handle)

        self._logger.info('....................\n')
        self._logger.info(f'CHOSEN PARAMETERS:\nRun name: {self.run_name} \nLogger level: {self.log_level} \n')
        self._logger.info('....................\n')

    def equation_list(self):
        print(self.equations)

    def get_data(self, metabolite):
        """
        Getting data from dataframes

        :param metabolite: name of metabolite to process
        :type metabolite: str

        """

        self._logger.info(f"Getting data for {metabolite}\n")

        # Filter by metabolite
        tmp_cal_df = self.cal_data[
            self.cal_data["compound"] == metabolite]
        tmp_sample_df = self.sample_data[
            self.sample_data["compound"] == metabolite]

        # One dictionary for calibration datas
        self.tmp_cal_dict["sources"] = tmp_cal_df["source"].to_list()
        self.tmp_cal_dict["x"] = tmp_cal_df["calibration concentration"].to_list()
        self.tmp_cal_dict['y'] = tmp_cal_df["M0/Mn"].to_list()

        # One dictionnary for sample data and predicted values
        self.tmp_sample_dict["sources"] = tmp_sample_df["source"].to_list()
        self.tmp_sample_dict["y_to_pred"] = tmp_sample_df["M0/Mn"].to_list()
        self.tmp_sample_dict['Linear result'] = []
        self.tmp_sample_dict['Quadratic result'] = []

        # We get the lowest and highest M0/Mn for determining if values to predict are in range
        y_limits = [min(self.tmp_cal_dict["y"]),
                    max(self.tmp_cal_dict["y"])]

        self._logger.debug(f"Checking limits for {metabolite}\n")

        for ind, val in enumerate(self.tmp_sample_dict["y_to_pred"]):
            if val < y_limits[0]:

                self.tmp_sample_dict["y_to_pred"][ind] = "Under range"

                self._logger.info(f"For {metabolite} the value to predict ({val}) is under range\n")

            elif val > y_limits[1]:

                self.tmp_sample_dict["y_to_pred"][ind] = "Over range"

                self._logger.info(f"For {metabolite} the value to predict ({val}) is over range\n")

    def get_residuals(self, cal_x, cal_y, equation):
        """
        Function to calculate residuals in calibration data from equation

        :param cal_x: x values from measured MS calibration data
        :type cal_x: class: numpy.array
        :param cal_y: y values from measured MS calibration data
        :type cal_y: class: numpy.array
        :param equation: polynomial equation in numpy one-dimensional polynomial class
        :type equation: class: numpy.poly1d

        :return: residual values for each x,y pair
        :rtype:list

        """

        residuals = []

        for x, y in zip(cal_x, cal_y):
            if x == 0 or y == 0:
                x = np.nan
                y = np.nan

            residuals.append((equation(x) - y) / y)

        return residuals

    def nan_clean(self, metabolite):
        """
        Check data for nans and process them

        :param metabolite: metabolite to process
        :type metabolite: str
        """
        # Checking values and removing NaNs
        nans = [ind for ind, val in enumerate(self.tmp_cal_dict["y"]) if np.isnan(val)]

        self._logger.debug(f'nans list: {nans}')
        self._logger.debug(
            f"Cal dict x: {self.tmp_cal_dict['x']}\n cal" 
            f"dict y: {self.tmp_cal_dict['y']}\n cal dict names: {self.tmp_cal_dict['sources']}\n")

        nancount = 0  # Counter to substract from indices

        for ind in nans:
            self._logger.info(
                f"Removing {self.tmp_cal_dict['sources'][ind - nancount]} from {metabolite} calibration datas because "
                "nan value detected")

            del (self.tmp_cal_dict['x'][ind - nancount])
            del (self.tmp_cal_dict['y'][ind - nancount])
            del (self.tmp_cal_dict['sources'][ind - nancount])
            nancount += 1

        if len(self.tmp_cal_dict['sources']) < 4:
            self._logger.error(
                f"Couldn't do calibration for {metabolite} because not enough calibration points (less than 4)")

        else:
            self._logger.debug('Cal dict x: {}\n cal dict y: {}\n cal dict names: {}\n'
                               .format(self.tmp_cal_dict['x'],
                                       self.tmp_cal_dict['y'],
                                       self.tmp_cal_dict['sources']))

    def calculate_equations(self, x, y, pol_order):
        """
        Calculate polynomial equations and r²

        :param x: x values from measured MS calibration data
        :type x: class: 'numpy.array'
        :param y: y values from measured MS calibration data
        :type y: class: 'numpy.array'
        :param pol_order: order for calculated polynomial equation
        :type pol_order: int
        :return convert: polynomial equation coefficients
        :rtype convert: class: 'numpy.array'
        :return equation: one-dimensional polynomial class to encapsulate 'natural' operations on polynomials
        :rtype equation: class: 'numpy.poly1D'
        :return r2: calculated r² for given polynomial equation
        :rtype r2: float

        """

        polyfit, stats = np.polynomial.polynomial.Polynomial.fit(x, y, pol_order, full=True)
        convert = list(polyfit.convert().coef)
        convert.reverse()
        equation = np.poly1d(convert)

        # Calculate r²
        yhat = equation(x)
        ybar = np.sum(y) / len(y)
        ssreg = np.sum((yhat - ybar) ** 2)
        sstot = np.sum((y - ybar) ** 2)
        r2 = round((ssreg / sstot), 4)

        return convert, equation, r2

    def build_polynome(self, metabolite):
        """Calculating value of coefficients of the quadratic polynomial

        :param metabolite: metabolite to process
        :type metabolite: str

        """

        # Equation part
        x = np.array(self.tmp_cal_dict["x"])
        y = np.array(self.tmp_cal_dict['y'])

        self._logger.debug(f"Performing calculations for equation\n x={x} \n y={y}")

        # Better than polyfit for numerical stability reasons (check docu for details)
        try:
            self.linear_convert, self.linear_equation, self.linear_r2 = self.calculate_equations(x, y, 1)
            self._logger.info(f'Linear r² = {self.linear_r2}\n')

        except Exception as e:
            self._logger.error(f'Error during calculation of linear equation for {metabolite}. Error: {e}\n')

        try:
            self.quad_convert, self.quad_equation, self.quad_r2 = self.calculate_equations(x, y, 2)
            self._logger.info(f'Quadratic r² = {self.quad_r2}\n')

        except Exception as e:
            self._logger.error(f'Error during calculation of quadratic equation for {metabolite}. Error: {e}\n')

        self._logger.info(f"Equations for {metabolite} are equal to:\nLinear: {self.linear_equation}\n ")
        self._logger.info(f'Quadratic: {self.quad_equation}\n')

        self._logger.debug("Calculating residuals \n")

        # We get relative residuals, will be used later for residual plot
        linear_residuals = self.get_residuals(x, y, self.linear_equation)
        quad_residuals = self.get_residuals(x, y, self.quad_equation)

        self.tmp_cal_dict["linear_relative_residuals"] = linear_residuals
        self.tmp_cal_dict["quad_relative_residuals"] = quad_residuals

        # Transform to list for deletions
        self.tmp_cal_dict["x"] = list(x)
        self.tmp_cal_dict["y"] = list(y)

    def build_plots(self, metabolite, pdf):
        """
        Coordinates the building of different plots and outputs them to pdf file

        :param metabolite: metabolite to process
        :type metabolite: str
        :param pdf: name of pdf file
        :type pdf: str

        """
        lin_cal_dict = copy.deepcopy(self.tmp_cal_dict)
        quad_cal_dict = copy.deepcopy(self.tmp_cal_dict)

        # Plotting part
        lin_ipass = 0  # Counter for removing residuals
        quad_ipass = 0

        lin_npass = 1  # Plot number
        quad_npass = 1

        # Create plots a first time
        fig, axs = plt.subplots(2, 2, figsize=(12, 10))

        self.plot_reg(lin_cal_dict,
                      metabolite,
                      self.linear_r2,
                      lin_npass,
                      1, 'linear',
                      ax=axs[0, 0])

        self.plot_reg(quad_cal_dict,
                      metabolite,
                      self.quad_r2,
                      quad_npass,
                      2, 'quadratic',
                      ax=axs[0, 1])

        self.plot_res(lin_cal_dict["x"],
                      lin_cal_dict["linear_relative_residuals"],
                      metabolite, self.linear_r2,
                      lin_npass, "linear",
                      ax=axs[1, 0])

        self.plot_res(quad_cal_dict["x"],
                      quad_cal_dict["quad_relative_residuals"],
                      metabolite, self.quad_r2,
                      quad_npass, "quadratic",
                      ax=axs[1, 1])

        fig.tight_layout()

        pdf.savefig(fig, bbox_inches='tight')

        # Remove data points if residuals > or < by 20%
        for lin_val in lin_cal_dict["linear_relative_residuals"]:
            if lin_val < -20 or lin_val > 20:
                lin_ipass += 1

        for quad_val in quad_cal_dict['quad_relative_residuals']:
            if quad_val < -20 or quad_val > 20:
                quad_ipass += 1

        # Keep plotting if some residuals are still too high or low
        while lin_ipass > 0 or quad_ipass > 0:

            fig, axs = plt.subplots(2, 2, figsize=(12, 10))

            if lin_ipass > 0:
                del (lin_cal_dict["linear_relative_residuals"][0])
                del (lin_cal_dict["x"][0])
                del (lin_cal_dict["y"][0])
                del (lin_cal_dict["sources"][0])

                lin_npass += 1

                self._logger.info(
                    "Removed {} data point {} because linear regression residual too large\n".format(
                        metabolite, lin_cal_dict["sources"][0])
                )

                self.plot_reg(lin_cal_dict,
                              metabolite,
                              self.linear_r2,
                              lin_npass,
                              1, 'linear',
                              ax=axs[0, 0])

                self.plot_res(lin_cal_dict["x"],
                              lin_cal_dict["linear_relative_residuals"],
                              metabolite, self.linear_r2,
                              lin_npass, "linear",
                              ax=axs[1, 0])

            if quad_ipass > 0:
                del (quad_cal_dict["quad_relative_residuals"][0])
                del (quad_cal_dict["x"][0])
                del (quad_cal_dict["y"][0])
                del (quad_cal_dict["sources"][0])

                quad_npass += 1

                self._logger.info(
                    "Removed {} data point {}  because quadratic regression residual too large\n".format(
                        metabolite, quad_cal_dict["sources"][0])
                )

                self.plot_reg(quad_cal_dict,
                              metabolite,
                              self.quad_r2,
                              quad_npass,
                              2, 'quadratic',
                              ax=axs[0, 1])

                self.plot_res(quad_cal_dict["x"],
                              quad_cal_dict["quad_relative_residuals"],
                              metabolite, self.quad_r2,
                              quad_npass, "quadratic",
                              ax=axs[1, 1])

            fig.tight_layout()

            pdf.savefig(fig, bbox_inches='tight')

            if lin_ipass > 0:
                lin_ipass -= 1
            if quad_ipass > 0:
                quad_ipass -= 1

    def plot_reg(self, cal_dict, metabolite, r2, npass, pol_order, reg, ax):
        """
        Function to plot the polynomial regression curve

        :param cal_dict: Dictionary containing calibration data
        :type cal_dict: dict
        :param metabolite: Metabolite to process
        :type metabolite: str
        :param r2: Determination coefficient for given metabolite and regression
        :type r2: float
        :param npass: Counter for number of plots created
        :type npass: int
        :param pol_order: Order of the polynomial
        :type pol_order: int
        :param reg: name of the regression
        :type reg: str
        :param ax: axis figure in which to put plot
        :type ax: class: 'matplotlib.axes.Axes'

        :return: Axes object containing the regression plot
        :rtype: class: 'matplotlib.axes.Axes'
        """

        regression_plot = sns.regplot(data=cal_dict,
                                      x=cal_dict['x'],
                                      y=cal_dict['y'],
                                      ci=None,
                                      line_kws={"color": 'red'},
                                      order=pol_order,
                                      ax=ax
                                      )

        regression_plot.set_xlabel("Concentration (µM)")
        regression_plot.set_ylabel("M0/Mn")
        regression_plot.set_title(f'{metabolite} {reg} regression plot\nR²={r2}\nPass no={npass}')

        return regression_plot

    def plot_res(self, x, residuals, metabolite, r2, npass, reg, ax):
        """
        Plot calculated residuals and refresh if residual too high or too low

        :param x: x values on which to plot residuals
        :type x: float
        :param residuals: residual values
        :type residuals: float
        :param metabolite: metabolite to process
        :type metabolite: str
        :param r2: determination coefficient for given plot
        :type r2: float
        :param npass: counter for number of plots created
        :type npass: int
        :param reg: name of the regression
        :type reg: str
        :param ax: axis figure in which to put plot
        :type ax: class: 'matplotlib.axes.Axes'

        :return: Axes object containing the residual plot
        :rtype: class: 'matplotlib.axes.Axes'
        """

        residual_plot = sns.scatterplot(x=x,
                                        y=residuals,
                                        ax=ax)

        # Get limits for plot y size depending on value of minimum and maximum relatiive residuals
        if min(residuals) < -30:
            min_ylim = min(residuals) * 1.25
        else:
            min_ylim = -30

        if max(residuals) > 30:
            max_ylim = max(residuals) * 1.25
        else:
            max_ylim = 30

        residual_plot.set_ylim(min_ylim, max_ylim)
        residual_plot.axhline(y=-20,
                              color='r')
        residual_plot.axhline(y=20,
                              color='r')
        residual_plot.axhline(y=0,
                              linestyle='--',
                              color='k')

        residual_plot.set_xlabel("X")
        residual_plot.set_ylabel("Y")
        residual_plot.set_title(f'{metabolite} {reg} residual plot\nR²={r2}\nPass no= {npass}\n')

        return residual_plot

    def predict_x_value(self, y_value):
        """
        Function to calculate concentration from value using the roots of the polynomial

        :param y_value: y value for which we predict x
        :type y_value: float

        :return quad_x_pred: x value returned by quadratic regression
        :rtype:float
        :return lin_x_pred: x value returned by linear regression
        :rtype lin_x_pred: float
        """

        # To get roots from polynomial we must have equation in form ax² + bx + c - y = 0
        # For linear predictions we use x = (y - b)/a

        if not np.isfinite(y_value):
            self._logger.error(f'Error: exp. value ({y_value}) must be a number.')

        else:

            nul_eq = self.quad_equation - y_value
            self._logger.debug(f'Roots are equal to: {nul_eq.roots}\n')
            roots = nul_eq.roots  # Remember we are solving for x values, not y

            self._logger.debug(f'minimum: {min(self.tmp_cal_dict["x"])}\n maximum: {max(self.tmp_cal_dict["x"])}\n')

            for ind, val in enumerate(roots):

                if min(self.tmp_cal_dict["x"]) <= val <= max(self.tmp_cal_dict["x"]):

                    quad_x_pred = val

                elif ind == (len(roots) - 1):

                    quad_x_pred = f'Roots out of range ({roots})'

            lin_x_pred = (y_value - self.linear_convert[1]) / self.linear_convert[0]

            if lin_x_pred < 0:
                lin_x_pred = 'Negative concentration (Out of range)'

            return quad_x_pred, lin_x_pred

    def calculate_predictions(self):
        """Function to coordinate predictions and check that datas are conform"""

        self._logger.debug("Predicting concentrations\n")

        for val in self.tmp_sample_dict["y_to_pred"]:
            self._logger.debug(f"Trying to predict {val} of type {type(val)} \n")

            if isinstance(val, str):

                if val == "Under range":
                    self.tmp_sample_dict["Linear result"].append("Under range")
                    self.tmp_sample_dict["Quadratic result"].append("Under range")

                elif val == "Over range":
                    self.tmp_sample_dict["Linear result"].append("Under range")
                    self.tmp_sample_dict["Quadratic result"].append("Over range")

                else:
                    self._logger.warning('Recieved non admissable value to predict: {} \n'.format(val))
                    continue

            elif np.isnan(val):
                self._logger.warning('NaN detected in values to predict')
                self.tmp_sample_dict["Linear result"].append(np.nan)
                self.tmp_sample_dict["Quadratic result"].append(np.nan)

            else:

                try:
                    quad_x_pred, lin_x_pred = self.predict_x_value(val)
                    self._logger.debug(
                        f'The calculated values are: \nQuadratic = {quad_x_pred}\nLinear = {lin_x_pred}\n')
                    self.tmp_sample_dict["Linear result"].append(lin_x_pred)
                    self.tmp_sample_dict["Quadratic result"].append(quad_x_pred)

                except Exception as e:
                    self._logger.error(
                        f'There was a problem while predicting the value from "{val}".\nError: {e}')

    def rebuild_dataframes(self, metabolite):
        """
        Function to assemble datas after calculations into final dataframe
        with a natural sort on index level "sources"

        :param metabolite: metabolite to process
        :type metabolite: str

        :return data_df: dataframe ccntaining calculated concentrations
        :rtype: class: 'Pandas.Dataframe'
        :return cal_df: dataframe containing data used for calibration curves
        :rtype cal_df: class: 'Pandas.Dataframe'
        """

        self._logger.debug("Rebuilding dataframes\n")

        data_df = pd.DataFrame.from_dict({"sources": self.tmp_sample_dict["sources"],
                                          "Linear Concentrations (µM)": self.tmp_sample_dict["Linear result"],
                                          "Quadratic Concentrations (µM)": self.tmp_sample_dict["Quadratic result"]})
        cal_df = pd.DataFrame.from_dict(self.tmp_cal_dict)

        # We sort values naturally by sources before putting in dfs
        data_df = data_df.sort_values(by="sources",
                                      key=lambda x: np.argsort(index_natsorted(x)))

        cal_df = cal_df.sort_values(by="sources",
                                    key=lambda x: np.argsort(index_natsorted(x)))

        # We need to insert the metabolite name before returning the dfs
        data_df = data_df.assign(metabolite=metabolite)
        cal_df = cal_df.assign(metabolite=metabolite)

        data_df.set_index(['metabolite', 'sources'], inplace=True)
        cal_df.set_index(['metabolite', 'sources'], inplace=True)

        return data_df, cal_df

    def main(self):
        """
        Main function to create a dataframe with quadratic polynomial
        predictions of concentrations from MS experiments
        """

        self._logger.info("Starting to process calculations...\n")

        list_of_cal_dfs, list_of_data_dfs = [], []

        with PdfPages(f'{self.run_name}.pdf') as pdf:

            for metabolite in self.cal_data["compound"].unique():

                self.get_data(metabolite)

                try:

                    self.nan_clean(metabolite)

                    self.build_polynome(metabolite)

                    self.build_plots(metabolite, pdf)

                    self.calculate_predictions()

                    data_df, cal_df = self.rebuild_dataframes(metabolite)

                    list_of_cal_dfs.append(cal_df)
                    list_of_data_dfs.append(data_df)

                except Exception as err:

                    self._logger.exception(f"There was a problem while calculating for {metabolite}\n Error: {err} \n ")

                    continue

        self.final_data_df, self.final_cal_df = pd.concat(list_of_data_dfs), pd.concat(list_of_cal_dfs)

        self.final_data_df.to_excel(r'Calculated datas.xlsx', index=True)
        self.final_cal_df.to_excel(r'Calibration datas.xlsx', index=True)

        self._logger.info("Done!")

        logging.shutdown()
