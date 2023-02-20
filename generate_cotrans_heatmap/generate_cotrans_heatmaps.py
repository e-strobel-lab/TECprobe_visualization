import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import argparse
import seaborn as sns
import re


# Class to call errors with specific messages to alert the user as to why the error was called.
class Error(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


parser = argparse.ArgumentParser(
                                description='Plot a heatmap of reactivities for either raw or normalized reactivities',
                                )
required = parser.add_argument_group(title='Required Optional')
# Argument for input file(s). Requires at least 1 filepath. Type should be string of the path to full, compiled data
# generated into .csv file from 'compile_sm2_output' script.
required.add_argument('--input_matrix',
                      '-i',
                      nargs='*',
                      type=str,
                      required=True,
                      help="""Input file(s) should be the full table of compiled SHAPEmapper2 data exported from the 
                      'compile_sm2_output' script. These .csv files end with 'full_table'"""
                      )
# Argument for text file that contains user-defined variables. Only takes one string input, which is the file path to
# the .txt file as a string.
required.add_argument('--variables',
                      '-v',
                      nargs=1,
                      type=str,
                      required=True,
                      help="""Input should be a .txt file (supplied: generate_cotrans_heatmap.txt) that contains the 
                      variables to be controlled by the user to generate custom heatmap matrices. Variables that can be 
                      controlled include: the reactivity data used to generate the heatmap, vmax values for drawing 
                      heatmaps, pallets to draw heat maps, the filetype used to save heatmaps, and various aspects of 
                      the heatmap graph (axes cutoffs, text size, etc.)."""
                      )
# Argument for output folder. Takes only one string input, which is the file path to the output folder.
required.add_argument('--output_folder',
                      '-o',
                      nargs=1,
                      type=str,
                      required=True,
                      help='Provide a folder where all generated data matrices and graphs should be saved to.'
                      )

args = parser.parse_args()
# Variables 'input_matrix', 'output_folder', and 'variables' reference lists which contain filepaths to each respective
# file(s)/folder(s).
input_matrix = args.input_matrix
output_folder = args.output_folder
variables = args.variables

# Ensures that all supplied folders and the variables .txt file exist and can be opened.
for i in range(0, len(input_matrix)):
    try:
        open_input = open(args.input_matrix[i])
        open_var = open(args.variables[0])
    except IOError:
        raise Error("Include the path to all files/folders.")

if os.path.isdir(output_folder[0]) is False:
    raise Error('Provide a valid output folder')

# Initialize an empty dictionary (variables_dict) to store variables from the generate_cotrans_heatmap_variables.txt
# file.
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
                    variables_dict[name] = None  # Stores value as None if user leaves the value for the variable blank
                else:
                    variables_dict[name] = value
                continue


# Initialize SetVariables class to store general user-defined variables in an easy to call method.
class SetVariables:
    def __init__(self, input_matrix):
        self.matrix_name = input_matrix
        self.matrix = pd.read_csv(self.matrix_name, index_col='Transcript Length', low_memory=False)
        self.rna_name = None
        self.probe_method = None
        self.chemical_probe = None
        self.ligand_conc = None
        self.ligand_name = None
        self.ligand = None
        self.remain = None
        self.concat_names = None
        self.set_names()
        self.filename = None
        self.set_filename()
        self.col = 'Reactivity_profile'  # Default = Raw reactivity
        self.set_react_col()

# Based on the filename (generated in the compile_sm2_output script), splits it to store the RNA name, the probe method
# (CoTxn or Equil), chemical probe, ligand concentration/name, and all remaining labels. Labels are generated based on
# index position.
    def set_names(self):
        filename_split = (os.path.normpath(self.matrix_name).split(os.sep)[-1]).split('_')
        self.rna_name = filename_split[0]
        self.probe_method = filename_split[1]
        self.chemical_probe = filename_split[2]
        self.ligand_conc = filename_split[3]
        self.ligand_name = filename_split[4]
        self.ligand = '_'.join([self.ligand_conc, self.ligand_name])
        self.remain = '_'.join(filename_split[5:-2])
        if self.remain == '':
            self.concat_names = '_'.join([self.rna_name, self.probe_method, self.chemical_probe, self.ligand])
        else:
            self.concat_names = '_'.join([self.rna_name, self.probe_method, self.chemical_probe, self.ligand,
                                          self.remain])

# Reads the saved value for the 'filename' variable in 'generate_cotrans_heatmap.txt'. Checks that the name contains
# only alphanumeric, -, or _ characters.
    def set_filename(self):
        save_file_regex = '^[A-Za-z0-9_-]+$'
        value = False
        # If the user does not include a value for 'filename', filename is set to contain the RNA name, probe method,
        # chemical probe, ligand name/concentration, and any 'remaining' labels from filename.
        if variables_dict['filename'] is None:
            self.filename = self.concat_names
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
                self.filename = no_space + '_' + self.concat_names
            else:
                raise Error('Choose a file name with only accepted characters')

    # Defines what columns are used for drawing the reactivity matrix. Default = 'Reactivity_profile'.
    def set_react_col(self):
        if variables_dict['col'] is None or variables_dict['col'].lower() == 'reactivity_profile':
            self.col = 'Reactivity_profile'
            pass
        elif variables_dict['col'].lower() == 'norm_calc_profile':
            self.col = 'Norm_calc_profile'  # Whole dataset normalized reactivities.
        elif variables_dict['col'].lower() == 'norm_profile':
            self.col = 'Norm_profile'  # SHAPEmapper2 normalized reactivities (normalized individually at each length)
                                       # where whole dataset normalization was not applied.
        elif variables_dict['col'].lower() == 'indiv_norm_profile':
            self.col = 'Indiv_Norm_profile'  # SHAPEmapper2 normalized reactivities (normalized individually at each
                                            # length). These are the same as 'Norm_profile', but have were renamed if
                                            # the user chose to calculate a whole dataset normalization factor. These
                                            # columns were renamed to avoid confusion with the other
                                            # ('Norm_calc_profile') normalized values.
        else:
            raise Error("""Column name is entered incorrectly. Make sure it is spelled correctly with proper "_" 
                        positioning, or leave blank to automatically plot "Reactivity_profile" values.""")


# Function to create (and save) 3 separate data matrices from 3 different types of columns of data from the full data
# table (input matrix). One each for untreated effective read depth, modified effective read depth, and the reactivity
# column chosen in the 'generate_cotrans_heatmap.txt' file and defined in the SetVariables class.
def split_matrix(matrix):
    df_list = []
    txpt_len = [i for i in range(1, matrix.matrix.shape[0] + 1)]
    matrix_df = pd.read_csv(matrix.matrix_name, index_col='Transcript Length', low_memory=False)
    col_matrix_list = ['Modified_effective_depth', 'Untreated_effective_depth']
    col_matrix_list.append(var[0].col)  # Adds user-defined reactivity column to list of columns to draw matrices of.
    for name in col_matrix_list:
        tl_col_names = {}
        indiv_df = matrix_df.filter(regex=name, axis=1)  # Filter all columns that contain one of the column names.
        if indiv_df.empty:
            raise Error("""Dataframe is empty. Check that you provided a valid reactivity column and have not 
            manipulated column names in the datatable. """)
        for col in indiv_df.columns:
            tl_col_names[col] = int(col[-3:])
        indiv_df.insert(0, 'Nucleotide', txpt_len)  # New column that represents all nucleotide positions.
        indiv_df2 = indiv_df.copy()
        indiv_df2.rename(columns=tl_col_names, inplace=True)  # Rename columns to only transcript lengths.
        indiv_df2.set_index('Nucleotide', inplace=True)  # Set nucleotide position as the index.
        indiv_df2 = indiv_df2.transpose()  # Switch axes so nucleotide positions are columns and transcript lengths
                                           # are rows.
        indiv_df2.index.rename('Transcript Length', inplace=True)  # Set the transcript lengths as the index.
        indiv_df2.sort_index(inplace=True)  # Sort transcript lengths numerically.
        indiv_df2.to_csv(os.path.join(output_folder[0], "{}_{}.csv").format(matrix.filename, name),
                         index_label='Transcript Length')  # Save dataframe to output folder.
        df_list.append(indiv_df2)  # Save all 3 dataframes as a list that is returned by the function.
    return df_list


# Class to set graphing parameters for matrices
class GraphVariables:
    def __init__(self, n, matrices):
        self.n = n  # Index position of dataset
        self.matrices = matrices
        self.full_length = True
        self.mod_depth_matrix_log = None
        self.untreated_depth_matrix_log = None
        self.mod_depth_matrix = self.matrices[0]
        self.untreated_depth_matrix = self.matrices[1]
        self.react_matrix = self.matrices[2]
        self.log_df()
        self.full_x = None
        self.xmin = None
        self.xmax = None
        self.x_window = None
        self.full_y = None
        self.ymin = None
        self.ymax = None
        self.y_window = None
        self.set_min_max()
        self.x_major = None
        self.x_pos = None
        self.x_minor = None
        self.y_major = None
        self.y_pos = None
        self.y_minor = None
        self.set_tick_labels()
        self.mod_depth_matrix_log_subset = None
        self.untreated_depth_matrix_log_subset = None
        self.react_matrix_subset = None
        self.subset_df()
        self.vmax_react_calc = None
        self.vmax_react()
        self.title = None
        self.set_title()
        self.palette_react = None
        self.palette_read = None
        self.set_palette()
        self.filetype = 'png'
        self.set_filetype()
        self.vmax_react = self.single_input('vmax_react', self.vmax_react_calc)  # Default = Reactivity at 95th
                                                                                 # percentile for raw reactivities
                                                                                 # and 2 for normalized reactivities.
        self.vmax_read = None
        self.set_vmax_read()
        self.x_tick_size = self.single_input('x_tick_size', 10)  # Default = 10
        self.y_tick_size = self.single_input('y_tick_size', 10)  # Default = 10
        self.x_label_size = self.single_input('x_label_size', 15)  # Default = 15
        self.y_label_size = self.single_input('y_label_size', 15)  # Default = 15
        self.title_size = self.single_input('title_size', 20)  # Default = 20
        self.colorbar_fontsize = self.single_input('colorbar_fontsize', 15)  # Default = 15
        self.fig_width, self.fig_length = self.double_input('fig_size', 10, 10)  # Default = 10, 10
        self.dpi = self.single_input('dpi', 'figure')  # Default is determined based on filetype.

# Take the log of values in effective read depth dataframes for simpler heatmap read out. These are also saved to the
# output folder.
    def log_df(self):
        self.mod_depth_matrix_log = self.mod_depth_matrix.copy()
        self.mod_depth_matrix_log.replace(0, 1, inplace=True)  # Replace value of 0 with 1 to avoid errors with log(0).
        self.untreated_depth_matrix_log = self.untreated_depth_matrix.copy()
        self.untreated_depth_matrix_log.replace(0, 1, inplace=True)  # Replace value of 0 with 1 to avoid errors with
                                                                     # log(0).
        for col in self.mod_depth_matrix_log.columns:
            self.mod_depth_matrix_log.loc[:, col] = np.log10(self.mod_depth_matrix_log.loc[:, col])
        for col in self.untreated_depth_matrix_log.columns:
            self.untreated_depth_matrix_log.loc[:, col] = np.log10(self.untreated_depth_matrix_log.loc[:, col])
        self.mod_depth_matrix_log.to_csv(os.path.join(output_folder[0],
                                                      var[self.n].filename + '_' + 'mod_depth_log.csv'))
        self.untreated_depth_matrix_log.to_csv(os.path.join(output_folder[0],
                                                            var[self.n].filename + '_' + 'untreated_depth_log.csv'))

    def set_min_max(self):
        self.full_x = [int(i) for i in self.mod_depth_matrix.columns.values]  # List of all nucleotide positions.
        if min(self.full_x) == 1:  # Inserts a 0 at index 0 for graphing purposes.
            self.full_x.insert(0, 0)
        self.full_y = [int(i) for i in self.mod_depth_matrix.index.values]  # List of all transcript lengths.
        # For both self.double_input() function calls, the function looks at the dictionary of user-inputted variables
        # (variables_dict) and sees if any values were entered for x-axis or y-axis cutoffs. If there are,
        # self.full_length is changed to False (in the double_input() function) and the min/max values for x & y axes
        # are returned. If no x/y-axis cutoffs were provided, the default values are listed
        # x-axis (0, max nucleotide position) and y-axis (20, max transcript length).
        self.xmin, self.xmax = self.double_input('x-axis', self.full_x[0], self.full_x[-1])
        self.ymin, self.ymax = self.double_input('y-axis', self.full_y[0], self.full_y[-1])
        # Ensures the value of self.xmax does not exceed the value of self.ymax. The graphing area will be blank for
        # those x values that exceed the value of self.ymax.
        if self.xmax > self.ymax:
            self.xmax = self.ymax
        # Ensures the value of self.ymin is not lower than the value of self.xmin. The graphing area will be blank for
        # those y values that fall below the value of self.xmin.
        if self.ymin < self.xmin:
            self.ymin = self.xmin
        # Creates a list of columns/nucleotide positions (self.x_window) and rows/transcript lengths (self.y_window).
        # that will be plotted.
        self.x_window = [i for i in range(self.xmin, self.xmax + 1)]
        self.y_window = [i for i in range(self.ymin, self.ymax + 1)]

    def set_tick_labels(self):
        self.x_major = [i for i in range(self.xmin, round(self.xmax, -1), 10)]  # Labels for x-major tick marks,
                                                                                # multiple of 10.
        # Generate positions for x-major tick marks in the middle of the cell they represent.
        if self.xmin == 0:
            self.x_pos = [i - 0.5 for i in self.x_major]
        else:
            self.x_pos = [i - self.xmin + 0.5 for i in self.x_major]
        if (self.xmax - self.xmin) // 10 == len(self.x_major):  # Add additional x-major label to end of graph.
            self.x_major.append(self.x_major[-1] + 10)
            self.x_pos.append(self.x_pos[-1] + 10)
        self.x_minor = [i + 5 for i in self.x_pos]  # Positions for x-minor tick marks (halfway between x-major marks)
        # Makes a list of y-axis positions (in multiples of 10) to label major tick marks.
        self.y_major = [i for i in range(self.ymin, round(self.ymax, -1), 10)]  # Labels for y-major tick marks,
                                                                                # multiple of 10.
        self.y_pos = [i - self.ymin + 0.5 for i in self.y_major]  # Positions for y-major tick marks, in middle of
                                                                  # cell they represent.
        self.y_minor = [i + 5 for i in self.y_pos]  # Positions for y-minor tick marks (halfway between y-major marks).
        if (self.ymax - self.ymin) // 10 == len(self.y_major):  # Add additional y-major label to end/bottom of graph.
            self.y_major.append(self.y_major[-1] + 10)
            self.y_pos.append(self.y_pos[-1] + 10)

    def subset_df(self):  # Filter out columns/rows if not plotting all data.
        if self.full_length is True:  # Set 'subset' versions of dataframes to the all nucleotides/transcript lengths.
            self.mod_depth_matrix_log_subset = self.mod_depth_matrix_log.copy()
            self.untreated_depth_matrix_log_subset = self.untreated_depth_matrix_log.copy()
            self.react_matrix_subset = self.react_matrix.copy()
        else:  # Create subsets of the dataframes using the windows provided.
            if self.x_window[0] == 0:
                self.x_window.remove(0)
            self.mod_depth_matrix_log_subset = self.mod_depth_matrix_log.loc[self.y_window, self.x_window]
            self.untreated_depth_matrix_log_subset = self.untreated_depth_matrix_log.loc[self.y_window, self.x_window]
            self.react_matrix_subset = self.react_matrix.loc[self.y_window, self.x_window]

    def vmax_react(self):  # Calculate vmax for reactivity heatmap.
        if var[0].col == 'Reactivity_profile':
            stack = self.react_matrix_subset.stack()
            self.vmax_react_calc = round(stack.quantile(q=0.95), 2)
        else:
            self.vmax_react_calc = 2  # Default vmax for normalized reactivity matrices = 2.

    # Vmax for effective read depths is determined by the largest value from both of the log dataframes. The
    # highest value between both dataframes is then set as the Vmax value for both effective read depth heatmaps (mod
    # and untreated)
    def set_vmax_read(self):
        max_values = []
        max_values.append(max(self.mod_depth_matrix_log.max(axis=1)))
        max_values.append(max(self.untreated_depth_matrix_log.max(axis=1)))
        max_log_read = round(max(max_values) + 0.5)
        self.vmax_read = self.single_input('vmax_read', max_log_read)

# Checks if the user entered a title. If not, the graph title is set to the RNA name + probe method + chemical probe +
# ligand name/concentration.
    def set_title(self):
        if variables_dict['title'] is None:
            self.title = var[self.n].concat_names
        # If the title contained ',', it was converted to a list, so the list is joined back into a single string to
        # set as the graph title.
        elif type(variables_dict['title']) is list:
            self.title = ','.join(variables_dict['title'])
        else:
            self.title = variables_dict['title']

    def set_palette(self):
        if variables_dict['palette_react'] is None:
            self.palette_react = 'viridis'
            variables_dict['palette_react'] = 'viridis'
        try:  # Check if the palette the user set is a valid seaborn palette.
            sns.color_palette(variables_dict['palette_react'])
            self.palette_react = variables_dict['palette_react']
        except ValueError:
            raise Error("Must be a valid seaborn palette")
        if variables_dict['palette_read'] is None:
            self.palette_read = 'magma'
            variables_dict['palette_read'] = 'magma'
        try:
            sns.color_palette(variables_dict['palette_read'])
            self.palette_read = variables_dict['palette_read']
        except ValueError:
            raise Error("Must be a valid seaborn palette")

    def set_filetype(self):
        if variables_dict['filetype'] is not None:
            filetype = variables_dict['filetype']
            if filetype in plt.gcf().canvas.get_supported_filetypes():  # Check that the supplied filetype is valid.
                self.filetype = filetype
            else:
                raise Error('Choose an accepted filetype')

    def single_input(self, key, default):  # Function to check if user has defined an option, if not sets to default.
        if variables_dict[key] is None:
            return default
        else:
            try:
                float_val = float(variables_dict[key])
                return float_val
            except ValueError:
                raise Error('Input must be a single int or float number.')

    def double_input(self, key, default1, default2):  # Function that checks if user has defined options, if not sets
                                                      # to deault1 and default2.
        if variables_dict[key] is None:
            return default1, default2
        elif key == 'x-axis' or key == 'y-axis':
            if len(variables_dict[key]) == 2:  # If two inputs were entered, checks that they can be converted to
                                               # integers.
                try:
                    temp_list = [int(i) for i in variables_dict[key]]
                    if min(temp_list) > default1:
                        min_val = min(temp_list)
                    else:
                        min_val = default1
                    if min_val == 1:  # Changes min value of 1 to 0 for easier graphing purposes.
                        min_val = 0
                    max_val = max(temp_list)
                    self.full_length = False
                    return min_val, max_val
                except ValueError:
                    raise Error('Inputs must be of int type.') from None
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


# Draw the heatmaps using all the parameters defined in the previous class.
def draw_graph(figsize, react_df, untreated_df, mod_df, react_vmax, read_vmax, react_palette, read_palette,
               colorbar_size, title, title_size, x_tick_pos, x_tick_labels, x_minor, x_tick_size, y_tick_pos,
               y_tick_labels, y_minor, y_tick_size, x_label_size, y_label_size, filename, filetype, dpi):
    for i in range(0, 3):
        fig, ax = plt.subplots(figsize=figsize)
        if i == 0:  # Reactivity heatmap
            sns.heatmap(react_df, vmin=0, vmax=react_vmax, cmap=react_palette, ax=ax, square=True,
                        cbar_kws={'shrink': 0.6})
            ax.collections[0].colorbar.set_label('Reactivity', labelpad=20, size=colorbar_size)
            ax.set_title(title + '_reactivity', fontsize=title_size, pad=20)
        if i == 1:  # Untreated effective read depth heatmap
            sns.heatmap(untreated_df, vmin=0, vmax=read_vmax, cmap=read_palette, ax=ax, square=True,
                        cbar_kws={'shrink': 0.6})
            ax.collections[0].colorbar.set_label('Log Read Depth', labelpad=20, size=colorbar_size)
            ax.set_title(title + '_untreated_effective_read_depth', fontsize=title_size, pad=20)
        if i == 2:  # Modified effective ready depth heatmap
            sns.heatmap(mod_df, vmin=0, vmax=read_vmax, cmap=read_palette, ax=ax, square=True,
                        cbar_kws={'shrink': 0.6})
            ax.collections[0].colorbar.set_label('Log Read Depth', labelpad=20, size=colorbar_size)
            ax.set_title(title + '_mod_effective_read_depth', fontsize=title_size, pad=20)
        ax.set_xticks(x_tick_pos, x_tick_labels)
        ax.set_xticks(x_minor, minor=True)
        ax.set_xticklabels(x_tick_labels, size=x_tick_size)
        ax.set_yticks(y_tick_pos, y_tick_labels)
        ax.set_yticks(y_minor, minor=True)
        ax.set_yticklabels(y_tick_labels, size=y_tick_size)
        ax.tick_params(axis='y', which='both', direction='inout')
        plt.xlabel('Nucleotide Position', size=x_label_size)
        plt.ylabel('Transcript Length', size=y_label_size)
        img = ax.get_figure()
        if i == 0:
            img.savefig(os.path.join(output_folder[0],
                        filename + "_Vmax{}_reactivity.{}").format(str(float(react_vmax)), str(filetype)),
                        format=filetype, dpi=dpi)
        if i == 1:
            img.savefig(os.path.join(output_folder[0],
                        filename + "_Vmax{}_untreated_effective_read_depth.{}".format(str(float(read_vmax)),
                                                                                      str(filetype))),
                        format=filetype, dpi=dpi)
        if i == 2:
            img.savefig(os.path.join(output_folder[0],
                        filename + "_Vmax{}_mod_effective_read_depth.{}".format(str(float(read_vmax)),
                                                                                str(filetype))),
                        format=filetype, dpi=dpi)
        plt.close()
        plt.figure().clear()


var = [SetVariables(i) for i in input_matrix]
matrices = [split_matrix(i) for i in var]
graph_var = [GraphVariables(n, i) for n, i in enumerate(matrices)]
for n, i in enumerate(graph_var):
    draw_graph(figsize=[i.fig_width, i.fig_length],
               react_df=i.react_matrix_subset,
               untreated_df=i.untreated_depth_matrix_log_subset,
               mod_df=i.mod_depth_matrix_log_subset,
               react_vmax=i.vmax_react,
               read_vmax=i.vmax_read,
               react_palette=i.palette_react,
               read_palette=i.palette_read,
               colorbar_size=i.colorbar_fontsize,
               title=i.title,
               title_size=i.title_size,
               x_tick_pos=i.x_pos,
               x_tick_labels=i.x_major,
               x_minor=i.x_minor,
               x_tick_size=i.x_tick_size,
               y_tick_pos=i.y_pos,
               y_tick_labels=i.y_major,
               y_minor=i.y_minor,
               y_tick_size=i.y_tick_size,
               x_label_size=i.x_label_size,
               y_label_size=i.y_label_size,
               filename=var[n].filename,
               filetype=i.filetype,
               dpi=i.dpi)
