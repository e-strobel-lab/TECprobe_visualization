import pandas as pd
import os
import numpy as np
import re
import argparse


# Class to call errors with specific messages to alert the user as to why the error was called.
class Error(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


parser = argparse.ArgumentParser(
    description="""Compile SHAPEmapper2 output nt profile files for each transcript length into a single .csv file"""
    )
required = parser.add_argument_group(title='Required Optional')

# Argument for input folders, requires at least 1 input folder, file path provided as a string.
required.add_argument('--input_folder',
                      '-i',
                      nargs='*',
                      type=str,
                      required=True,
                      help="""Each input folder should contain a single dataset. Due to how TECprobe-ML operates, each
                      transcript length is processed through SHAPEmapper2 separately, thus each transcript length will
                      have its own folder output from SHAPEmapper2. A folder which contains all transcript length
                      folders from a single dataset is the required format of the input folder for this script.
                      Multiple datasets can be processed at once if each folder path is supplied. The script will not
                      run as intended if a single folder with multiple datasets is provided."""
                      )

# Argument for text file that contains user-defined variables. Only takes one string input, which is the file path to
# the .txt file.
required.add_argument('--variables',
                      '-v',
                      nargs=1,
                      type=str,
                      required=True,
                      help="""Input should be a .txt file (supplied) that contains the variables to be controlled by 
                              the user. These variables allow the user to control: whether whole dataset normalization
                              is applied, the cutoff values of background mutation rate and effective read depth for 
                              inclusion of a reactivity in normalization calculation, and what the reverse primer
                              sequence of the DNA template is (used for determining which transcripts are not enriched
                              """
                      )

# Argument for output folder. Takes only one string input, which is the file path to the output folder.
required.add_argument('--output_folder',
                      '-o',
                      nargs=1,
                      type=str,
                      required=True,
                      help='Provide a folder path where all output compiled tables should be saved to.'
                      )

args = parser.parse_args()
# Variables 'input_folder', 'output_folder', and 'variables' reference lists that contain filepaths to each respective
# file(s)/folder(s)
input_folder = args.input_folder
output_folder = args.output_folder
variables = args.variables

# Ensures that all supplied folders and the variables .txt file exist and can be opened.
for i in range(0, len(input_folder)):
    if not os.path.isdir(input_folder[i]):
        raise Error('Provide a valid directory for folder(s) containing SHAPEmapper2 output files.')
try:
    open(args.variables[0])
except IOError:
    raise Error("Include the path to the compile_sm2_output_variables.txt document.")

if not os.path.isdir(output_folder[0]):
    raise Error('Provide a valid output folder')

# Initialize an empty dictionary (variables_dict) to store variables from the compile_sm2_output_variables.txt file.
variables_dict = {}
with open(variables[0], 'r') as f:
    for line in f:
        if line.startswith('>'):
            strip_line = line.strip('>')
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


# Initialize SetVariables class to store general user-defined variables in an easy to call method. Also defines
# 'self.col', which references all columns of interest in the code.
class SetVariables:
    def __init__(self):
        self.col = ['Reactivity_profile', 'Norm_profile', 'Modified_effective_depth',
                    'Untreated_effective_depth', 'Untreated_rate']
        self.norm = None  # Defines whether the whole dataset normalization factor will be applied.
        self.normalization()
        self.bkg_mut = None  # Defines the background mutation cutoff.
        self.background_mutation()
        self.mod_eff_depth = None  # Defines the modified effective read depth cutoff.
        self.effective_depth()

    # Defines whether normalization will be recalculated for the whole data set.
    def normalization(self):
        if variables_dict['normalization'] is None or variables_dict['normalization'].lower() == 'true':
            self.norm = True
        elif variables_dict['normalization'].lower() == 'false':
            self.norm = False
        else:
            raise Error('Normalization must be set to "True" or "False", or left blank.')

    # Defines the background mutation cutoff to be applied when selecting reactivities to be used for calculating the
    # normalization factor
    def background_mutation(self):
        if variables_dict['bkg_mut'] is None:
            self.bkg_mut = 0.05
        elif variables_dict['bkg_mut'] is not None:
            try:
                self.bkg_mut = float(variables_dict['bkg_mut'])
            except ValueError:
                raise Error('Number must be of int or float type.')

    # Defines the modified effective read depth cutoff to be applied when selecting reactivities to be used for
    # calculating the normalization factor.
    def effective_depth(self):
        if variables_dict['mod_eff_depth'] is None:
            self.mod_eff_depth = 5000
        elif variables_dict['mod_eff_depth'] is not None:
            try:
                self.mod_eff_depth = int(variables_dict['mod_eff_depth'])
            except ValueError:
                raise Error('Number must be of int type.')


# Function that travels the path of each input folder to find the relevant SHAPEmapper2 output file (nt_profile.txt).
# A dictionary is created where key is the transcript length and the value is the path to the nt_profile.txt file for
# that trasncript length, which contains SHAPEmapper2 calculations of reactivities, read depth, etc. This dictionary is
# returned by the function.
def file_dictionary(inputs):
    file_dict = {}
    for (root, dirs, files) in os.walk(inputs):
        for entry in files:
            if entry.endswith("profile.txt"):
                txp_len = entry[-18:-15]
                if txp_len.startswith('_'):
                    txp_len = txp_len.replace("_", '0')
                file_dict[txp_len] = os.path.join(root, entry)
    return file_dict


# Creates a class to store information specific to each set of data. Information stored is: 1. Whether the probing
# was done cotranscriptionally or at equilibrium, 2. The name of the RNA, 3. The name and concentration of the
# ligand (if applicable), and 4. The chemical probe.
class DatasetName:
    def __init__(self, dictionary):
        self.dictionary = dictionary
        self.name_list = []  # Contains name/labels of interest when filepath name is split
        self.probe_method = None
        self.chemical_probe = None  # Defines the chemical probe
        self.cotx_equil()
        self.rna_name = self.name_list[0]
        self.ligand_name = None  # Defines the ligand name
        self.ligand_conc = None  # Defines the ligand concentration
        self.ligand_name_conc = None  # Concatenates ligand name and concentration
        self.remain = None  # Defines any remaining part of the folder name for file naming purposes
        self.ligand()

# Stores one of the transcript length filepaths to parse to determine whether filename contains 'CoTxn'
# (cotranscriptional) or 'Equil' (equilibrium). The index of this is then used to determine the RNA name and the
# chemical probe used.
    def cotx_equil(self):
        for key, value in self.dictionary.items():
            while len(self.name_list) == 0:
                filename = os.path.normpath(value).split(os.sep)
                self.name_list = filename[-2].split('_')
                if 'CoTxn' in self.name_list:
                    cotx_equl_idx = self.name_list.index('CoTxn')  # Return index position of 'CoTxn'
                    self.probe_method = 'CoTxn'
                elif 'Equil' in self.name_list:
                    cotx_equl_idx = self.name_list.index('Equil')  # Return index position of 'Equil'
                    self.probe_method = 'Equil'
                else:
                    raise Error('Filename of dataset incorrectly formatted')
                if cotx_equl_idx > 1:  # Join all items before the CoTxn/Equil tag into a single string
                    self.name_list[0] = '_'.join(self.name_list[0:cotx_equl_idx])
                    del self.name_list[1:cotx_equl_idx]
                    cotx_equl_idx = 1
                else:
                    pass
                # Chemical probe always listed after Cotxn/Equil
                self.chemical_probe = self.name_list[cotx_equl_idx + 1]

# Uses regex to find the ligand concentration by searching for 3 digits in a row followed by any letter and then a
# capital M. If a match is found, the index of the ligand concentration is used to find the ligand name and these two
# values are stored separately (ligand_conc and ligand_name, respectively), and also joined together to store as a
# single variable (ligand_name_conc). If no match is found, these variables are all stored as "No Ligand".
    def ligand(self):
        lig_conc_regex = '[0-9][0-9][0-9][a-z][M]'
        r = re.compile(lig_conc_regex)
        if any(match := r.match(i) for i in self.name_list):
            ligand_conc_unmod = match.group(0)
            ligand_conc_split = re.split('([a-z])', ligand_conc_unmod)
            ligand_conc_split[0] = str(int(ligand_conc_split[0]))
            self.ligand_conc = ''.join(ligand_conc_split)
            ligand_conc_idx = self.name_list.index(ligand_conc_unmod)
            self.ligand_name = self.name_list[ligand_conc_idx - 1]
            self.ligand_name_conc = self.ligand_conc + '_' + self.ligand_name
            self.remain = '_'.join(self.name_list[ligand_conc_idx + 1:-2])
        else:
            self.ligand_conc = 'No_ligand'
            self.ligand_name = 'No_ligand'
            self.ligand_name_conc = 'No_ligand'


# Function to parse through the dictionary returned from file_dictionary() and returns a dataframe for each dataset.
def csv_to_dataframe(csv_dict):
    df_list = []
    for key, value in csv_dict.items():
        drop_list = []
        df = pd.read_csv(value, delimiter='\t')
        # Drops all rows that contain lowercase nucleotides in the 'Sequence' column (which indicates
        # non-riboswitch/RNA of interest sequence). In our files, this represents the SC1 hairpin, which is 5' to
        # the riboswitch and is used for facilitating sample processing.
        for idx, row in df.iterrows():
            if row['Sequence'].islower():
                drop_list.append(idx)
            else:
                continue
        df.drop(labels=drop_list, axis=0, inplace=True)
        last_row = df.index[-1]
        # Some final nucleotide positions have a value of 0, while some have a value of NaN. NaN will not plot for the
        # heatmap, and will make the diagonal appear jagged. This code looks for NaN values at the last nucleotide
        # position of each transcript length and changes to 0 (for Reactivity_profile, Norm_profile, Modified_rate,
        # Untreated_rate, and Denatured_rate columns) or 1 (for Modified_effective_depth, Untreated_effective_depth, and
        # Denatured_effective_depth).
        for i in var.col:
            if i in df.columns:
                if 'profile' in i.lower() or 'rate' in i.lower():
                    if np.isnan(df.at[last_row, i]):
                        df.at[last_row, i] = 0
                elif 'effective' in i.lower():
                    if np.isnan(df.at[last_row, i]):
                        df.at[last_row, i] = 1
            if i in df.columns:
                df[i].fillna(-9999, inplace=True)  # Fill all remaining NaN with -9999. Used as marker for user.
        # Appends the transcript length to the end of each column name
        rename_col = {}
        for col in df.columns:
            rename_col[col] = '{}_{}'.format(col, key)
        df.rename(columns=rename_col, inplace=True)
        df_list.append(df)
    # Concatenates the list of dataframes to return a consolidated dataframe containing all of the data columns for a
    # single experiment in a single table.
    df_react = pd.concat(df_list, axis=1)
    df_react.index = np.arange(1, len(df_react) + 1)
    return df_react


# Defines the reverse primer used to amplify the DNA template. This information is used to drop certain transcript
# lengths from the dataframe used to calculate the normalization factor. Certain lengths are dropped because they are
# not enriched in the dataset due to a lack of internal biotin modifications in the primer and proximity to the
# terminal biotin modification on the 5' end of the reverse primer.
class RevPrimer:
    def __init__(self, matrix):
        self.df = matrix
        self.rev_primer = None  # Defines the reverse primer sequence
        self.len_rev_primer = None  # Defines the length of self.rev_primer
        self.set_rev_primer()
        self.txpt_len = None  # Defines length of full length RNA
        self.rev_primer_idx_start = None  # Defines the nucleotide position at which self.rev_primer sequence starts
        self.set_rev_primer_index()
        self.df_drop = None  # Stores dataframe with transcript lengths dropped from subset_df() function
        self.subset_df()

    # Function to determine if the user has supplied a reverse primer sequence to remove certain transcript lengths.
    # If this sequence is supplied and does not match our HP4 reverse primer sequence, the sequence is checked to
    # make sure it only contains RNA nucleotides (A, G, U, or C) and then sets it as the rev_primer sequence. The
    # length of the reverse primer sequence is also stored (len_rev_primer).
    def set_rev_primer(self):
        rna_nt = ['G', 'C', 'U', 'A']
        if variables_dict['rev_primer'] is None or \
                variables_dict['rev_primer'].upper() == 'UGAUUCGUCAGGCGAUGUGUGCUGGAAGACAUU':
            self.rev_primer = 'UGAUUCGUCAGGCGAUGUGUGCUGGAAGACAUU'  # HP4 sequence
        else:
            for i in variables_dict['rev_primer']:
                if i.upper() in rna_nt:
                    continue
                else:
                    raise Error(
                        'Ensure sequence is in RNA format and contains only "G", "C", "U", and "A" characters.')
            self.rev_primer = variables_dict['rev_primer'].upper()
        self.len_rev_primer = len(self.rev_primer)

    # Function to create a dictionary (seq_dict) that stores the transcript length as the key and the entire sequence
    # of that RNA in the value. The values are then searched to see all transcript lengths where the full rev_primer
    # sequence is found.
    def set_rev_primer_index(self):
        self.txpt_len = [i for i in range(1, self.df.shape[0] + 1)]
        seq_dict = self.df['Sequence_{}'.format(self.txpt_len[-1])].transpose().to_dict()
        for k, v in seq_dict.items():
            try:
                seq_dict[k] = ''.join([seq_dict[k - 1], v])
            except KeyError:
                continue
        idx_list = []
        for k, v in seq_dict.items():
            if re.findall(self.rev_primer, v):
                idx_list.append(k)
            else:
                continue
        rev_primer_idx_end = idx_list[0]  # Identify first length where full rev_primer sequence found
        self.rev_primer_idx_start = (rev_primer_idx_end - self.len_rev_primer + 1)  # Define transcript length where
                                                                                    # rev_primer sequence starts

    # Function to define transcript lengths that need to be dropped from the dataframe used for normalization factor
    # calculation purposes. A list is made of transcript lengths spanning 10 before the self.rev_primer_idx_start
    # position to the final transcript length. Transcript lengths between 15 and 8 from the 3' end are dropped from
    # this list. The remaining lengths in the list are used to remove the corresponding columns of data from the
    # dataframe.
    def subset_df(self):
        idx_drop = []
        for i in range(self.rev_primer_idx_start - 10, self.txpt_len[-1] + 1):  # Creates a list of all transcript
                                                                                # lengths from 10 nucleotides upstream
                                                                                # of the reverse primer sequence to the
                                                                                # 3' end of the template
            idx_drop.append(i)
        for i in range(idx_drop[-15], idx_drop[-7]):  # Drops lengths 15 to 8 from 3' end;
                                                     # these transcript lengths will be retained in the dataframe for
                                                     # normalization factor calculation purposes. All other lengths
                                                     # in the 'idx_drop' list are the values that will be removed from
                                                     # the dataframe used to calculate the normalization factor
            idx_drop.remove(i)
        col_drop = []
        for i in idx_drop:  # Create a list of all column names to drop based on whether the column contains the
                            # transcript length indicated in the idx_drop list
            for col in self.df.columns:
                if int(col[-3:]) == i:
                    col_drop.append(col)
        self.df_drop = self.df.drop(col_drop, axis=1)  # Drops all columns from col_drop list


# Function to calculate the normalization factor. This is adapted from the original SHAPEmapper2 script
# (https://github.com/Weeks-UNC/shapemapper2) with slight modifications.
def normalization_factor(df_drop, data_d):
    df_drop.replace(0, np.nan, inplace=True)  # Replace 0's with NaN so they do not contribute to
                                              # normalization factor calculation
    for key, value in data_d.items():
        for idx, row in df_drop.iterrows():
            # Check if the "Modified_effective_depth" column exists for that transcript length (some transcript
            # lengths have been dropped from the dataframe (see: subset_df()) in the df_drop dataframe.
            if 'Modified_effective_depth_{}'.format(key) in df_drop.columns:
                # If the transcript length still exists in the dataframe, check whether the
                # 'Modified_effective_depth' or 'Untreated_effective_depth' is below the user-defined cutoff
                # (mod_eff_depth) or 'Untreated_rate' is above the user-defined cutoff (bkg_mut) for each
                # nucleotide position (idx) of a certain transcript length (key). If any of these requirements
                # are met, the 'Reactivity_profile' value for that transcript length (key) at that nucleotide
                # position (idx) is set to NaN, so that it does not contribute to calculation of the normalization
                # factor.
                if row['Modified_effective_depth_{}'.format(key)] < var.mod_eff_depth \
                        or row['Untreated_effective_depth_{}'.format(key)] < var.mod_eff_depth \
                        or row['Untreated_rate_{}'.format(key)] > var.bkg_mut:
                    df_drop.at[idx, 'Reactivity_profile_{}'.format(key)] = np.nan
                else:
                    continue
            else:
                continue
    df_norm_filtered = df_drop.filter(regex='Reactivity_profile')  # Filter the dataframe for only the
                                                                   # 'Reactivity_profile' columns
    num_txpt_len = len(df_norm_filtered.columns)  # Number of transcript lengths present in filtered dataframe
    stacked_df = df_norm_filtered.stack()  # Stacks all columns into a single column
    reactivity_list = sorted(stacked_df.tolist())  # Converts numerically ordered reactivities to a list
    len_list = len(reactivity_list)  # Calculates the total number of reactivities remaining
    q25, q75 = np.nanpercentile(reactivity_list, [25, 75])  # Determine reactivity value at 25th and 75th percentiles
    iqr_limit = 1.5 * abs(q75 - q25)  # Calculate the interquartile range of all reactivities
    top_10 = len_list // 10  # Index position of where top 10% of len_list begins
    top_5 = len_list // 20  # Index position of where top 5% of len_list begins
    ten_limit = reactivity_list[-top_10 - 1]  # Highest reactivity after 'removing' top 10% reactivities
    five_limit = reactivity_list[-top_5 - 1]  # Highest reactivity after 'removing' top 5% reactivities
    # Determines the max value (limit) between iqr_limit and ten_limit (if the number of transcript lengths in the
    # dataframe is > or = 100) or five_limit (if the number of transcript lengths in the dataframe is < 100).
    limit = max(iqr_limit, ten_limit)
    if num_txpt_len < 100:
        limit = max(iqr_limit, five_limit)
    # Creates a list of reactivities (filter_list) that fall below the value of "limit"
    filter_list = []
    for i in range(len_list):
        if reactivity_list[i] < limit:
            filter_list.append(reactivity_list[i])
    a = 0
    # Average together the top 10% of reactivities of 'filter_list' to generate the normalization factor.
    try:
        for i in range(-top_10, 0):
            a += filter_list[i]
        norm_factor = a / top_10
    except IndexError:
        raise Error("Unable to calculate a normalization factor.")
    return norm_factor


# Function to save a .csv file of the full, concatenated dataframe. If the user wanted the new, calculated normalized
# reactivity and std err values, those columns are created and appended before saving.
def matrix_generation(df_react, norm_factor, data_dict, data_vars):
    if var.norm is True:
        norm_col_dict = {}
        for k, v in data_dict.items():
            # Divides all raw reactivity and std err values by the normalization factor and stores these as columns
            # Norm_calc_profile and Norm_calc_stderr, respectively.
            norm_col_dict['Norm_calc_profile_{}'.format(k)] = \
                df_react['Reactivity_profile_{}'.format(k)] / norm_factor
            norm_col_dict['Norm_stderr_{}'.format(k)] = df_react['Std_err_{}'.format(k)] / norm_factor
            # To avoid confusion between the two sets of normalized columns, the normalized columns generated from
            # SHAPEmapper2 are renamed to have 'Indiv' added to the beginning of the column name
            df_react.rename(columns={'Norm_profile_{}'.format(k): 'Indiv_Norm_profile_{}'.format(k),
                                     'Norm_stderr_{}'.format(k): 'Indiv_Norm_stderr_{}'.format(k)}, inplace=True)
        norm_values = pd.DataFrame(norm_col_dict)
        full_df = pd.concat([df_react, norm_values], axis=1)  # Create single dataframe with all info and new,
                                                              # calculated normalized values
        full_df.sort_index(axis=1, inplace=True)  # Sort columns alphabetically
        # Full data matrix saved by naming with the RNA name, chemical probe, and ligand name/concentration.
        full_df.to_csv(os.path.join(output_folder[0],
                                    "{}_{}_{}_{}_Full_Table.csv").format(data_vars.rna_name,
                                                                         data_vars.probe_method,
                                                                         data_vars.chemical_probe,
                                                                         data_vars.ligand_name_conc,
                                                                         data_vars.remain
                                                                         ),
                       index_label='Transcript Length'
                       )
    else:  # Sort the dataframe columns alphabetically and save as .csv if the user doesn't want to add normalized
           # reactivities
        df_react.sort_index(axis=1, inplace=True)
        df_react.to_csv(os.path.join(output_folder[0],
                                     "{}_{}_{}_{}_Full_Table.csv").format(data_vars.rna_name,
                                                                          data_vars.probe_method,
                                                                          data_vars.chemical_probe,
                                                                          data_vars.ligand_name_conc,
                                                                          data_vars.remain
                                                                          ),
                        index_label='Transcript Length'
                        )


var = SetVariables()
data_dict = [file_dictionary(i) for i in input_folder]
data_vars = [DatasetName(i) for i in data_dict]
df_react = [csv_to_dataframe(i) for i in data_dict]
rev = [RevPrimer(i) for i in df_react]
norm_factor = [normalization_factor(i.df_drop, j) for i, j in zip(rev, data_dict)]
matrices = [matrix_generation(w, x, y, z) for w, x, y, z in zip(df_react, norm_factor, data_dict, data_vars)]
