import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import numpy as np

ztp_keep1 = [i for i in range(20, 130 + 1)]
ztp_keep2 = [i for i in range(159, 162 + 1)]
ztp_keep = ztp_keep1 + ztp_keep2
flu_keep1 = [i for i in range(20, 82 + 1)]
flu_keep2 = [i for i in range(113, 116 + 1)]
flu_keep = flu_keep1 + flu_keep2
pg_keep1 = [i for i in range(20, 137 + 1)]
pg_keep2 = [i for i in range(169, 172 + 1)]
pg_keep = pg_keep1 + pg_keep2


class Error(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


parser = argparse.ArgumentParser(description="""This script makes a hexbin plot of two cotranscriptional probing 
                                datasets. Each nucleotide in each transcript length is plotted with its x and y values 
                                represented by the reactivities in the first and second datasets supplied, respectively.
                                Depth of overlapping data are shown by color related to the supplied colorbar. A
                                Pearson's Correlation Coefficient (r) is calculated from the whole dataset comparison
                                and plotted on the graph.""")
required = parser.add_argument_group(title='Required Optional')
# Argument for input files. Requires 2 filepaths, one for each replicate to be compared. Type should be string of the
# path to the .csv reactivity matrix of user's choice generated from the 'generate_cotrans_heatmap' script.
required.add_argument('--input_matrix', '-i',
                    nargs=2,
                    type=str,
                    required=True,
                    help="""Input requires two .csv files of replicates to compare. These should be reactivity matrix
                    .csv files generated from the 'generate_cotrans_heatmap' script. Variables that can be controlled
                    include: whether the riboswitch being assessed currently has ranges of transcripts stored
                    in the script to perform calculations on (currently this script has the Clostridium beijerinckii 
                    ZTP (21-130 and 159-162), Bacillus cereus fluoride (21-82 and 113-116), and Clostridiales bacterium 
                    oral taxon 876 str. F0540 ppGpp (21-137, 169-172) riboswitches) (whether these are stored or not, 
                    the user also has the option to define their own ranges of transcript lengths to assess) and general 
                    graph parameters (x/y-axis limits, graph title, figure size, filetype of saved graph etc.)""")
# Argument for text file that contains user-defined variables. Only takes one string input, which is the file path to
# the .txt file as a string.
required.add_argument('--variables', '-v',
                    nargs=1,
                    type=str,
                    required=True,
                    help="""Input should be a .txt file (supplied: plot_cotrans_correlation_variables) that contains the
                    variables to be controlled by the user to generate the hexbin correlation plots.""")
# Argument for output folder. Takes only one string input, which is the file path to the output folder.
required.add_argument('--output_folder', '-o',
                    nargs=1,
                    type=str,
                    required=True,
                    help="""Provide a folder where all generated data matrices should be saved to.""")

args = parser.parse_args()
# Variables 'input_matrix', 'output_folder', and 'variables' reference lists which contain filepaths to each respective
# file(s)/folder.
input_matrix = args.input_matrix
output_folder = args.output_folder
variables = args.variables

# Ensures that all supplied files and the output folder exist and can be opened.
try:
    open(args.input_matrix[0])
    open(args.input_matrix[1])
    open(args.variables[0])
except IOError:
    raise Error("Include a working path to all required files.")

if os.path.isdir(output_folder[0]) == False:
    raise Error('Provide a valid output folder')

# Initialize an empty dictionary (variables_dict) to store variables from the plot_cotrans_correlation_variables.txt
variables_dict = {}
with open(variables[0], 'r') as f:
    for line in f:
        if line.startswith('>'):
            strip_line = line.strip('> ')
            name, value = strip_line.replace('\n', '').split("=")
            if ',' in value:
                value_list = value.split(',')
                variables_dict[name] = value_list
                continue
            else:
                if value == '':
                    variables_dict[name] = None  # Stores value as None if user leaves the value for the variable blank.
                else:
                    variables_dict[
                        name] = value
                continue


# Initialize SetVariables class to store general user-defined variables in an easy to call method.
class SetVariables:
    def __init__(self, matrix):
        self.matrix = matrix
        self.set_len = []
        self.ribo_stored = True
        self.user_txpt_len()
        self.name = None
        self.txpt_len = None
        self.input_name()
        self.df = None
        self.matrix_df()
        self.stacked_df = None
        self.min_val = None
        self.max_val = None
        self.get_min_max_val()
        self.gridsize = self.single_input('gridsize', 100)
        self.filename = None
        self.set_filename()

    # Create a list of transcript lengths to plot nucleotide reactivities of on the hexbin plot. This function is only
    # run if the user sets 'ribo_stored' to False (in the plot_cotrans_correlation_variables.txt file) and provides both
    # the number of transcript ranges to plot, and those ranges.
    def user_txpt_len(self):
        num_range = 0
        txpt_keep = []
        if variables_dict['ribo_stored'] is None:
            variables_dict['ribo_stored'] = 'True'
        if variables_dict['ribo_stored'].lower() == 'false':
            self.ribo_stored = False
            try:
                int(variables_dict['num_range'])
                num_range = int(variables_dict['num_range'])
            except ValueError:
                raise Error("Must be a single value of int type.")
            except TypeError:
                raise Error("If 'ribo_stored' is set to False, you must provide an integer value for 'num_range'")
        if self.ribo_stored is False:
            try:
                for i in variables_dict['txpt_keep']:
                    int(i)
                for i in variables_dict['txpt_keep']:
                    txpt_keep.append(int(i))
            except ValueError:
                raise Error("Must be values of int type.")
            except TypeError:
                raise Error("If 'ribo_stored' is set to False, you must provide integer values for 'txpt_keep'")
            for i in range(1, num_range + 1):
                min_num1 = min(txpt_keep)
                txpt_keep.remove(min_num1)
                min_num2 = min(txpt_keep)
                txpt_keep.remove(min_num2)
                for j in range(min_num1, min_num2 + 1):
                    self.set_len.append(j)

    def input_name(self):
        name_split = os.path.normpath(self.matrix).split(os.sep)[-1].split('_')
        self.name = '_'.join(name_split[0:7])
        for i in name_split:
            if self.ribo_stored is True:
                if 'crcbf' in i.lower():
                    self.txpt_len = flu_keep
                elif 'g4p' in i.lower():
                    self.txpt_len = pg_keep
                elif 'ztp' in i.lower():
                    self.txpt_len = ztp_keep
            else:
                self.txpt_len = self.set_len

    # Remove the transcript lengths not specified by the script or user.
    def matrix_df(self):
        self.df = pd.read_csv(self.matrix)
        self.df.rename_axis('Nucleotide', axis=1, inplace=True)
        self.df.set_index('Transcript Length', inplace=True)
        index_list = self.df.index.to_list()
        drop_list = []
        for i in index_list:
            if i not in self.txpt_len:
                drop_list.append(i)
        self.df.drop(labels=drop_list, axis=0, inplace=True)

    # Calculate and store the min/max reactivity value from each dataset to determine axes for the graph.
    def get_min_max_val(self):
        stacked = self.df.stack()
        self.stacked_df = stacked.to_frame(name=self.name)
        self.min_val = self.stacked_df[self.name].min()
        self.max_val = self.stacked_df[self.name].max()

    # Reads the saved value for the 'filename' variable in 'plot_cotrans_correlation.txt'. Checks that the name
    # contains only alphanumeric, -, or _ characters.
    def set_filename(self):
        save_file_regex = '^[A-Za-z0-9_-]+$'
        value = False
        # If the user does not include a value for 'filename', filename is set to contain the RNA name, probe method,
        # chemical probe, ligand name/concentration, and appends the gridsize to the end.
        if variables_dict['filename'] is None:
            self.filename = self.name + "_" + str(self.gridsize)
            pass
        # If the value for 'filename' contained commas (,), the entry was stored as a list in the script and is now
        # concatenated with '_'.
        elif type(variables_dict['filename']) is list:
            value = '_'.join(variables_dict['filename'])
        else:
            value = variables_dict['filename']
        # If the user did input a filename, regex is used to check the filename is suitable and appends the RNA name,
        # probe method, chemical probe, ligand name/concentration, and 'remaining labels' to the end.
        if value is False:
            pass
        else:
            no_space = value.replace(' ', '_')  # Replace spaces with underscores.
            regex_save = re.compile(save_file_regex)
            if regex_save.search(no_space):  # If only alphanumeric characters are used, filename is stored.
                self.filename = no_space + "_" + self.name + "_" + str(self.gridsize)
            else:
                raise Error('Choose a file name with only accepted characters')

    def single_input(self, key, default):
        if variables_dict[key] == None:
            return default
        else:
            try:
                int_val = int(variables_dict[key])
                return int_val
            except ValueError:
                raise Error('Input must be a single int.')


# Concatenate the two replicates into a single dataframe. The individual dataframes are stacked so a single column
# represents all reactivities from a single replicate.
def concat_matrices(matrix_1, matrix_2):
    concat_df = pd.concat([matrix_1.stacked_df, matrix_2.stacked_df], axis=1)
    con_col = concat_df.columns.values.tolist()
    num_con_col = []
    # Numbers the replicates if the files for the two datasets have the same value for self.name.
    for n, v in enumerate(con_col):
        totalcount = con_col.count(v)
        count = con_col[:n].count(v)
        num_con_col.append(v + "_" + str(count + 1) if totalcount > 1 else v)
    concat_df.columns = num_con_col  # Replace the column names with numbered versions.
    concat_df.replace(0, np.nan, inplace=True)  # Replace all 0's with NaN.
    concat_df.loc[concat_df.isnull().any(axis=1), :] = np.nan  # Locate all rows that have a NaN and replace the other
                                                               # value with NaN (if not already).
    concat_df.dropna(inplace=True)
    concat_df.to_csv(os.path.join(output_folder[0], matrix_1.filename + '.csv'))
    return concat_df

# Class to set graphing variables for generating the hexbin plot.
class GraphVariables:
    # noinspection PyTypeChecker
    def __init__(self, concat_df):
        self.concat_df = concat_df
        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        self.set_min_max()
        self.x_minor, self.x_major, self.x_major_label = self.set_tick_labels(self.x_min, self.x_max)
        self.y_minor, self.y_major, self.y_major_label = self.set_tick_labels(self.y_min, self.y_max)
        self.title = None
        self.set_title()
        self.palette = 'viridis'
        self.set_palette()
        self.filetype = 'png'
        self.set_filetype()
        self.vmax = self.single_input('vmax', 500)  # Default = 500
        self.x_tick_size = self.single_input('x_tick_size', 10)  # Default = 10
        self.y_tick_size = self.single_input('y_tick_size', 10)  # Default = 10
        self.x_label_size = self.single_input('x_label_size', 10)  # Default = 10
        self.y_label_size = self.single_input('y_label_size', 10)  # Default = 10
        self.title_size = self.single_input('title_size', 15)  # Default = 15
        self.fig_width, self.fig_length = self.double_input('fig_size', 10, 10)  # Default = 10, 10
        self.dpi = self.single_input('dpi', 'figure') # Default is determined based on filetype.

    # Use either the min/max values determined from the data, or values chosen by the user to set the bounds of the
    # x-/y-axes for the hexbin plot.
    def set_min_max(self):
        calc_min = min(matrix_variable[0].min_val, matrix_variable[1].min_val)
        calc_max = max(matrix_variable[0].max_val, matrix_variable[1].max_val)
        # These are to increase the limits of the min/max values to include some white space on the graph so that the
        # data points represented by the min/max values are not right on the edge of the graph.
        if calc_min < 0:
            if calc_min < -1:  # If the lowest reactivity is below -1, subtract 0.1 and round to the nearest 1st
                               # decimal place.
                calc_min_ax = round(calc_min - 0.1, 1)
            elif calc_min > -0.1:  # If the lowest reactivity is between -0.1 and 0, subtract 0.001 and round to the
                                 # nearest 5th decimal place.
                calc_min_ax = round(calc_min - 0.001, 4)
            else:  # If the reactivity is between more negative than -0.1, subtract 0.01 and round to the nearest 3rd
                   # decimal place.
                calc_min_ax = round(calc_min - 0.01, 2)
        elif calc_min == 0:
            calc_min_ax = 0
        else:
            if calc_min > 1:  # If the lowest reactivity is greater than 1, subtract 0.1 and round to the nearest 1st
                             # decimal place.
                calc_min_ax = round(calc_min - 0.1, 1)
            elif calc_min < 0.1:  # If the lowest reactivity is between 0 and 0.1, subtract 0.001 and round to the
                                  # nearest 5th decimal place.
                calc_min_ax = round(calc_min - 0.001, 4)
            else:  # If the lowest reactivity is between 0.1 and 1, subtract 0.01 and round to the nearest 3rd decimal
                   # place.
                calc_min_ax = round(calc_min - 0.01, 2)
        if calc_max > 1:  # If the highest reactivity is greater than 1, add 0.1 and round to the nearest 1st decimal
                          # place.
            calc_max_ax = round(calc_max + 0.1, 1)
        elif calc_max == 0:
            calc_max_ax = 0
        elif calc_max > 0.1:  # If the highest reactivity is between 0.1 and 1, add 0.01 and round to the nearest 3rd
                              # decimal place.
            calc_max_ax = round(calc_max + 0.001, 2)
        else:  # If the highest reactivity is between 0 and 0.1, add 0.01 and round to the nearest 3rd decimal place.
            calc_max_ax = round(calc_max + 0.001, 4)
        self.x_min, self.x_max = self.double_input('x-axis', calc_min_ax, calc_max_ax)
        self.y_min, self.y_max = self.double_input('y-axis', calc_min_ax, calc_max_ax)

    # Calculate and set tick labels/positions.
    def set_tick_labels(self, min_v, max_v):
        if abs(min_v) == abs(max_v):  # If the absolute value of min and max are the same.
            steps = round(float(max_v) / 2, 2)  # Increments are max divided by 2.
        else:  # Otherwise, subtract min value from max value and divide by 3.
            steps = round((float(max_v) - float(min_v)) / 3, 2)
        major_num = [round(i, 2) for i in np.arange(min_v, max_v + steps, steps)]
        minor_num = [round(i, 2) + steps / 2 for i in major_num]
        del minor_num[-1]  # Ensures max edge of graph has a major tick mark.
        major_label = ['{:.2f}'.format(round(i, 2)) for i in major_num]
        return minor_num, major_num, major_label

    # Checks if user entered a title. If not, the title is set to the filename generated in the SetVariables class.
    def set_title(self):
        if variables_dict['title'] is None:
            self.title = matrix_variable[0].filename
        # If the title contained ',', it was converted to a list, so the list is joined back into a single string to
        # set as the graph title.
        elif type(variables_dict['title']) is list:
            self.title = ','.join(variables_dict['title'])
        else:
            self.title = variables_dict['title']

    def set_palette(self):
        try:  # Check if the palette the user set is a valid seaborn palette.
            sns.color_palette(variables_dict['palette'])
            self.palette = variables_dict['palette']
        except ValueError:
            raise Error("Must be a valid seaborn palette")

    def set_filetype(self):
        if variables_dict['filetype'] is not None:
            filetype = variables_dict['filetype']
            if filetype in plt.gcf().canvas.get_supported_filetypes():  # Check that the supplied filetype is valid.
                self.filetype = filetype
            else:
                raise Error('Choose an accepted filetype')

    def single_input(self, key, default):
        if variables_dict[key] == None:
            return default
        else:
            try:
                float_val = float(variables_dict[key])
                return float_val
            except ValueError:
                raise Error('Input must be a single int number.')

    def double_input(self, key, default1, default2):
        if variables_dict[key] == None:
            return default1, default2
        elif key == 'x-axis' or 'y-axis':
            if len(variables_dict[key]) == 2:
                try:
                    temp_list = [float(i) for i in variables_dict[key]]
                    min_val = min(temp_list)
                    max_val = max(temp_list)
                    return min_val, max_val
                except ValueError:
                    raise Error('Inputs must be of int or float type.') from None
            else:
                raise Error('Must supply 2 int or float values separated by a comma.')
        elif key == 'fig_size':
            if len(variables_dict[key]) == 2:
                try:
                    temp_list = [float(i) for i in variables_dict[key]]
                    val1 = temp_list[0]
                    val2 = temp_list[1]
                    return val1, val2
                except ValueError:
                    raise Error('Inputs must be of int or float type.') from None
            else:
                raise Error('Must supply 2 int or float values separated by a comma.')

# Draw the hexbin plot using the defined parameters.
def draw_graph(figsize, df, gridsize, cmap, vmax, xlim, ylim, title, title_size, x_major, x_major_label, x_tick_size,
               x_minor, y_major, y_major_label, y_minor, y_tick_size, col1, col2, x_label_size, y_label_size, filename,
               filetype, dpi):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect(aspect='equal')
    # Draw horizontal and vertical red dashed lines at 0 position to clearly show where x and y = 0.
    ax.axhline(0, color='red', lw=.3, linestyle='dashed', alpha=0.3, zorder=0)
    ax.axvline(0, color='red', lw=.3, linestyle='dashed', alpha=0.3, zorder=1)
    plot = plt.hexbin(df.iloc[:,0], df.iloc[:,1], gridsize=gridsize, cmap=cmap, mincnt=1, edgecolors='None', vmax=vmax)
    cbar_major_ticks = [i for i in np.arange(0, vmax+1, 100)]
    cbar_minor_ticks = [i + 50 for i in cbar_major_ticks[:-1]]
    cbar_major_ticks[0] = 1
    cbar = plt.colorbar(label='# of Instances', ax=ax, ticks=cbar_major_ticks)
    cbar.ax.yaxis.set_ticks(cbar_minor_ticks, minor=True)
    cbar.ax.tick_params(labelsize=6)
    ax.set(xlim=xlim, ylim=ylim)
    ax.annotate("r = {:.3f}".format(np.corrcoef(df.iloc[:, 0], df.iloc[:, 1])[0, 1]),
                xy=(.05, .95), xycoords='axes fraction')
    ax.set_title(title, fontsize=title_size)
    ax.set_xticks(x_major, x_major_label)
    ax.set_xticks(x_minor, minor=True)
    ax.set_xticklabels(x_major_label, size=x_tick_size)
    ax.set_yticks(y_major, y_major_label)
    ax.set_yticks(y_minor, minor=True)
    ax.set_yticklabels(y_major_label, size=y_tick_size)
    plt.xlabel(col1, fontsize=x_label_size)
    plt.ylabel(col2, fontsize=y_label_size)
    plt.savefig(os.path.join(output_folder[0], filename + '.{}'.format(filetype)),
                format=filetype, dpi=dpi, bbox_inches="tight")
    plt.close()
    plt.figure().clear()

matrix_variable = [SetVariables(i) for i in input_matrix]
concat_matx = concat_matrices(matrix_variable[0], matrix_variable[1])
graph_variables = GraphVariables(concat_matx)
draw_graph(figsize=[graph_variables.fig_width, graph_variables.fig_length],
           df=concat_matx,
           gridsize=matrix_variable[0].gridsize,
           cmap=graph_variables.palette,
           vmax=graph_variables.vmax,
           xlim=[graph_variables.x_min, graph_variables.x_max],
           ylim=[graph_variables.y_min, graph_variables.y_max],
           title=graph_variables.title,
           title_size=graph_variables.title_size,
           x_major=graph_variables.x_major,
           x_major_label=graph_variables.x_major_label,
           x_minor=graph_variables.x_minor,
           y_minor=graph_variables.y_minor,
           x_tick_size=graph_variables.x_tick_size,
           y_major=graph_variables.y_major,
           y_major_label=graph_variables.y_major_label,
           y_tick_size=graph_variables.y_tick_size,
           col1=matrix_variable[0].name,
           col2=matrix_variable[1].name,
           x_label_size=graph_variables.x_label_size,
           y_label_size=graph_variables.y_label_size,
           filename=matrix_variable[0].filename,
           filetype=graph_variables.filetype,
           dpi=graph_variables.dpi)

