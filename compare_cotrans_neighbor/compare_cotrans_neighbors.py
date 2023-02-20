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


parser = argparse.ArgumentParser(description="""This script calculates the Pearson Correlation Coefficient (r) for
                                            Reactivity_profile and Untreated_rate values of neighboring (n and n+1)
                                            transcripts. The  resulting r values are plotted as a violin plot (one plot 
                                            each for Reactivity_profile and Untreated_rate) with the number of 
                                            comparisons (counts) and data median reported.""")
required = parser.add_argument_group(title='Required Optional')
# Argument for input file(s). Requires at least 1 filepath. Type should be string of the path to full, compiled data
# generated into .csv file from 'compile_sm2_output' script.
required.add_argument('--input_matrix', '-i',
                    nargs='*',
                    type=str,
                    required=True,
                    help="""Input file(s) should be the full table of compiled SHAPEmapper2 data exported from the 
                      'compile_sm2_output' script. These .csv files end with 'full_table'. Files should also be from
                      the same RNA/riboswitch, or at least the same transcript lengths.""")
# Argument for text file that contains user-defined variables. Only takes one string input, which is the file path to
# the .txt file as a string.
required.add_argument('--variables', '-v',
                    nargs=1,
                    type=str,
                    required=True,
                    help="""Input should be a .txt file (supplied: compare_cotrans_neighbors.txt) that contains the 
                    variables to be controlled by the user to generate the violin plots. Variables that can be 
                    controlled include: whether the riboswitch being assessed currently has ranges of transcripts stored
                    in the script to perform calculations on (currently this script has the Clostridium beijerinckii 
                    ZTP (21-130 and 159-162), Bacillus cereus fluoride (21-82 and 113-116), and Clostridiales bacterium 
                    oral taxon 876 str. F0540 ppGpp (21-137, 169-172) riboswitches) and general graph parameters 
                    (y-axis limits, graph title, figure size, filetype of saved graph, etc.)
                    """)
# Argument for output folder. Takes only one string input, which is the file path to the output folder.
required.add_argument('--output_folder', '-o',
                    nargs=1,
                    type=str,
                    required=True,
                    help='Provide a folder where generated graphs should be saved to.')

args = parser.parse_args()
# Variables 'input_matrix', 'output_folder', and 'variables' reference lists which contain filepaths to each respective
# file(s)/folder.
input_matrix = args.input_matrix
output_folder = args.output_folder
variables = args.variables

# Ensures that all supplied files and the output folder exist and can be opened.
for i in range(0, len(input_matrix)):
    try:
        open_input = open(args.input_matrix[i])
        open_var = open(args.variables[0])
    except IOError:
        raise Error("Include a working path to all required files.")

if os.path.isdir(output_folder[0]) is False:
    raise Error('Provide a valid output folder')

# Initialize an empty dictionary (variables_dict) to store variables from the compare_cotrans_neighbors_variables.txt
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
                    variables_dict[name] = None  # Stores value as None if user leaves the value for the variable blank.
                else:
                    variables_dict[name] = value
                continue


# Initialize SetVariables class to store general user-defined variables in an easy to call method.
class SingleVariable:
    def __init__(self):
        self.title = None
        self.set_title()
        self.filename_react = None
        self.filename_bkg_mut = None
        self.set_filename()
        self.filetype = 'png'
        self.set_filetype()
        self.title_size = self.single_input('title_size', 10)  # Default = 10
        self.fig_width, self.fig_length = self.double_input('fig_size', 10, 10)  # Default = 10, 10
        self.dpi = self.single_input('dpi', 'figure')  # Default is determined based on filetype.
        self.ymin, self.ymax = self.double_input('ylim', 0, 1)
        self.filename = None
        self.set_filename()

    # Checks if the user entered a title, if not the graph title is left empty.
    def set_title(self):
        if variables_dict['title'] is None:
            self.title = ''
            pass
        # If the title contained ',', it was converted to a list, so the list is joined back into a single string to
        # set as the graph title.
        elif type(variables_dict['title']) is list:
            self.title = ','.join(variables_dict['title']) + "_"
        else:
            self.title = variables_dict['title'] + "_"

    # Reads the saved value for the 'filename' variable in 'compare_cotrans_neighbors.txt'. Checks that the name
    # contains only alphanumeric, -, or _ characters.
    def set_filename(self):
        save_file_regex = '^[A-Za-z0-9_-]+$'
        value = False
        # If the user does not include a value for 'filename', filename is left empty and will be named later based
        # on information from a single input filename.
        if variables_dict['filename'] is None:
            self.filename = ''
        # If the value for 'filename' contained commas (,), the entry was stored as a list in the script and is now
        # concatenated with '_'.
        elif type(variables_dict['filename']) is list:
            value = '_'.join(variables_dict['filename'])
        else:
            value = variables_dict['filename']
        # If the user did input a filename, regex is used to check the filename is suitable.
        if value is False:
            pass
        else:
            no_space = value.replace(' ', '_')  # Replace spaces with underscores
            regex_save = re.compile(save_file_regex)
            if regex_save.search(no_space):  # If only alphanumeric characters are used, filename is stored.
                self.filename = no_space + "_"
            else:
                raise Error('Choose a file name with only accepted characters')

    def set_filetype(self):
        if variables_dict['filetype'] is not None:
            filetype = variables_dict['filetype']
            if filetype in plt.gcf().canvas.get_supported_filetypes():  # Check that the supplied filetype is valid.
                self.filetype = filetype
            else:
                raise Error('Choose an accepted filetype')

    def single_input(self, key, default):
        if variables_dict[key] is None:
            return default
        else:
            try:
                int_val = int(variables_dict[key])
                return int_val
            except ValueError:
                raise Error('Input must be a single int number.')

    def double_input(self, key, default1, default2):  # Function that checks if user has defined options, if not sets
                                                      # to deault1 and default2
        if variables_dict[key] is None:
            return default1, default2
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
        elif key == 'ylim':
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
        else:
            raise Error('Must supply 2 int or float values separated by a comma.')


class SetVariables:
    def __init__(self, matrix):
        self.split_file = matrix.split('/')
        self.matrix = pd.read_csv(matrix, low_memory=False)
        self.read_df = None
        self.react_df = None
        self.bkg_mut_df = None
        self.df_manip()
        self.set_len = []
        self.ribo_stored = True
        self.user_txpt_len()
        self.riboswitch = None
        self.name_react = None
        self.name_bkg_mut = None
        self.riboswitch = None
        self.txpt_len = None
        self.set_name()

    # Filter the full data table(s) into individual matrices for the 3 columns of interest: Modified_effective_depth,
    # Reactivity_profile, and Untreated_rate
    def df_manip(self):
        cols = ['Modified_effective_depth', 'Reactivity_profile', 'Untreated_rate']
        txpt_len = [i for i in range(1, self.matrix.shape[0] + 1)]
        df_list = []
        for i in cols:
            col_dict = {}
            temp_df = self.matrix.filter(regex=i, axis=1)
            for col in temp_df.columns:
                col_dict[col] = int(col[-3:])
            temp_df.insert(0, 'Nucleotide', txpt_len)  # New column that represents all nucleotide positions.
            temp_df2 = temp_df.copy()
            temp_df2.rename(columns=col_dict, inplace=True)
            temp_df2.set_index('Nucleotide', inplace=True)
            temp_df2 = temp_df2.transpose()
            temp_df2.index.rename('Transcript Length', inplace=True)
            temp_df2.sort_index(inplace=True)
            df_list.append(temp_df2)
        self.read_df = df_list[0]
        self.react_df = df_list[1]
        self.bkg_mut_df = df_list[2]

    # Create a list of transcript lengths to perform Pearson's Correlation Coefficient calculation on if ranges were
    # provided by the user.
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

    # Set names for reactivity and background mutation data based on information from filename. Will be used in the
    # title of the graph and the filename to save the graph.
    def set_name(self):
        name_split = self.split_file[-1].split('_')
        self.name_react = '_'.join(name_split[:-2]) + '_react'
        self.name_bkg_mut = '_'.join(name_split[:-2]) + '_bkg_mut'
        for i in name_split:
            if self.ribo_stored is True:
                if 'crcbf' in i.lower():
                    self.riboswitch = 'Flu'
                    self.txpt_len = flu_keep
                elif 'g4p' in i.lower():
                    self.riboswitch = 'ppGpp'
                    self.txpt_len = pg_keep
                elif 'ztp' in i.lower():
                    self.riboswitch = 'ZTP'
                    self.txpt_len = ztp_keep
            else:
                self.riboswitch = name_split[0]
                self.txpt_len = self.set_len


react_dict = {}
bkg_mut_dict = {}


def gather_r(matrix):
    all_react = []
    all_bkg_mut = []
    for i in matrix.read_df.index.tolist():
        if i == matrix.read_df.index.tolist()[0]:  # Require a transcript length at n-1 for comparison.
            pass
        elif i not in matrix.txpt_len or i-1 not in matrix.txpt_len:
            pass
        else:  # Calculate Pearson Correlation Coefficient (r) for all nucleotide positions, except the last position
               # in transcript n, as it does not exist in transcript n-1.
            pearson_react = np.corrcoef(matrix.react_df.loc[i, : i - 1],
                                        matrix.react_df.loc[i - 1, : i - 1])[0, 1]
            pearson_bkg_mut = np.corrcoef(matrix.bkg_mut_df.loc[i, : i - 1],
                                          matrix.bkg_mut_df.loc[i - 1, : i - 1])[0, 1]
            all_react.append(pearson_react)
            all_bkg_mut.append(pearson_bkg_mut)
    #  Append to dictionary with the identifier of the dataset as key and a list r values for the value.
    react_dict[matrix.name_react] = all_react
    bkg_mut_dict[matrix.name_bkg_mut] = all_bkg_mut


def r_df(dictionary, figsize, ymin, ymax, title, title_size, riboswitch, filename, filetype, dpi):
    medians = []
    counts = []
    df = pd.DataFrame.from_dict(dictionary)
    for key, value in dictionary.items():  # Calculate median of r values and number of comparisons performed.
        medians.append('median: ' + str(round(df[key].median(), 2)))
        counts.append('n: ' + str(df[key].count()))
    melt_df = pd.melt(df, value_name='r', var_name='Dataset')
    fig, ax = plt.subplots(figsize=figsize)
    sns.violinplot(x='Dataset', y='r', data=melt_df, hue='Dataset', cut=0, ax=ax, linewidth=0.4,
                   palette=['#A967AA', '#5BB75B'])
    plt.ylim(ymin, ymax)
    pos = range(len(counts))
    for tick, label in zip(pos, ax.get_xticklabels()):
        ax.text(pos[tick], 1, counts[tick] + '\n' + medians[tick],
                horizontalalignment='center',
                size='xx-small',
                color='black',
                weight='semibold')
    ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    ax.set_title(str(title) + riboswitch, fontsize=title_size, pad=20)
    fig.tight_layout()
    img = ax.get_figure()
    img.savefig(os.path.join(output_folder[0], filename + riboswitch + ".{}".format(str(filetype))), format=filetype,
                dpi=dpi)
    plt.figure().clear()


single_var = SingleVariable()
var = [SetVariables(i) for i in input_matrix]
rs = [gather_r(i) for i in var]
swarm_react = r_df(dictionary=react_dict,
                   figsize=[single_var.fig_width, single_var.fig_length],
                   ymin=single_var.ymin,
                   ymax=single_var.ymax,
                   title=single_var.title,
                   title_size=single_var.title_size,
                   riboswitch=var[0].name_react,
                   filename=single_var.filename,
                   filetype=single_var.filetype,
                   dpi=single_var.dpi)
swarm_bkg_mut = r_df(dictionary=bkg_mut_dict,
                     figsize=[single_var.fig_width, single_var.fig_length],
                     ymin=single_var.ymin,
                     ymax=single_var.ymax,
                     title=single_var.title,
                     title_size=single_var.title_size,
                     riboswitch=var[0].name_bkg_mut,
                     filename=single_var.filename,
                     filetype=single_var.filetype,
                     dpi=single_var.dpi)
