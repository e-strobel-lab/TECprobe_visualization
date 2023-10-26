# TECprobe_visualization

TECprobe_visualation is a suite of tools for visualizing data from Transcription Elongation Complex RNA structure probing (TECprobe) experiments that has been processed through the TECtools and ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2). These visualizations include:

**compile_sm2_output** - Concatenates ShapeMapper2 output files for each transcript length of RNA and normalizes nucleotide reactivity using the whole dataset.

**generate_cotrans_heatmap** - Generates heatmaps of nucleotide reactivity and effective read depth (both modified and untreated samples).

**plot_cotrans_correlation** - Compares raw nucleotide reactivity values from two experimental replicates and plots the comparison as a hexbin plot.

**compare_cotrans_neighbor** - Calculates the Pearson Correlation Coefficient (r) when comparing raw reactivity or background mutation rate values of neighboring transcript lengths at each nucleotide position. These r values are then assembled into violin plots showing the distribution of r values for each parameter comparison.

## Processing ShapeMapper2 data using compile_sm2_output
`compile_sm2_output` compiles `ShapeMapper2` (https://github.com/Weeks-UNC/shapemapper2) output files that end in `profile.txt` into a single .csv file. These .txt files contain all data relevant to reactivity data processing and analysis.

### compile_sm2_output inputs

Required: 

```
-i/--input_folder <input_folder>    Path to folder which contains ShapeMapper2 processed data to be
                                    concatenated. Accepts 1+ folder paths

-v/--variables <variables_txt>      Path to the variables file. This will implement different options
                                    during the run (see below)

-o/--output_folder <output_folder>  Path to the folder to save output file(s)

```

### Basic usage of compile_sm2_output:
  1. `input_folder` filepath(s) should be to a folder that contains subfolders of each individual transcript length, which contains output files from ShapeMapper2 analysis. The specific files `compile_sm2_output` is looking for are the ShapeMapper2 outputs that end with `nt_profile.txt`.
  2. The `variables_txt` file provided to `compile_sm2_output` should be filled in to the users specifications. This file is provided in `TECprobe_visualization/compile_sm2_output` as `compile_sm2_output_variables.txt`. Each variable that can be set contains additional information about that variable, accepted arguments for that variable, as well as the default value (if applicable).
  3. It is important to note that because ShapeMapper2 (https://github.com/Weeks-UNC/shapemapper2) is run on individual transcript lengths, and not the dataset as a whole, the normalization calculation reported by ShapeMapper2 is not calculated in the context of all reactivities for the dataset. The user can choose to have the same normalization calculation performed by ShapeMapper2 run when all data for a given dataset have been compiled. This will add new reactivity (`Norm_calc_profile`) and standard error (`Norm_stderr`) columns  and rename the previous 'normalized' data columns to `Indiv_calc_profile` and `Indiv_stderr`.
  4. `compile_sm2_output` generates a single .csv file for each dataset provided. The .csv file(s) will be saved to the `output_folder` and named similar to how the `nt_profile.txt` files are labeled with 'Full_Table' appended to the end.
  5. `compile_sm2_output` should be run in terminal as follows:

     `python <path to compile_sm2_output.py> -i <input_folder> -v <variables_txt> -o <output_folder>`

## Generating heatmaps using generate_cotrans_heatmap

`generate_cotrans_heatmap` generates heatmaps and corresponding .csv files for reactivity, untreated effective read depth, and modified effective read depth. 

### generate_cotrans_heatmap inputs

Required: 

```
-i/--input_matrix <input_matrix_csv>    Path to Full_Table.csv file(s) generated from running
                                        compile_sm2_output. Accepts 1+ file(s)

-v/--variables <variables_txt>          Path to the variables file. This will implement different options
                                        during the run (see below)

-o/--output_folder <output_folder>      Path to the folder to save output file(s)

```

### Basic usage of generate_cotrans_heatmap

1. `input_matrix_csv` filepath(s) should lead to Full_Table.csv file(s) generated from `compile_sm2_output`.
2. The `variables_txt` file provided to `generate_cotrans_heatmap` should be filled in to the users specifications. This file is provided in `TECprobe_visualization/generate_cotrans_heatmap` as `generate_cotrans_heatmap_variables.txt`. Each variable that can be set contains additional information about that variable, accepted arguments for that variable, as well as the default value (if applicable).
3. Multiple choices exist for reactivity data to assemble into a heatmap. If the user set `Normalization = True` when running `compile_sm2_output`, the user can choose the following columns for generating a reactivity heatmap: `Reactivity_profile` (background-subtracted mutation rate), `Norm_calc_profile` (Reactivities normalized with a single, whole-dataset calculated normalization factor), or `Indiv_Norm_profile` (Reactivities normalized with individual normalization factors for each transcript length). If the user set `Normalization = False` when running `compile_sm2_output`, the user can choose the following columns for generating a reactivity heatmap: `Reactivity_profile` (background-subtracted mutation rate) or `Norm_profile` (The same as `Indiv_Norm_profile` when `Normalization = True`, but this column is not renamed if the user sets `Normalization = False`).
4. Eight total files are generated for each `input_matrix_csv` provided:
   * One drawn heatmap of user-chosen reactivity parameter; filetype determined by user.
   * One .csv file of the reactivity matrix used to draw heatmap; File saved with name of the user-chosen reactivity column appended to the end of the dataset name.
   * One drawn heatmap of untreated effective read depth; filetype determined by user.
   * Two .csv files of the untreated effective read depth matrices: 1. Raw reads and 2. Log reads. Files saved with `Untreated_effective_depth` and `untreated_depth_log`, respectively.
   * One drawn heatmap of modified effective read depth; filetype determined by user.
   * Two .csv files of the modified effective read depth matrices: 1. Raw reads and 2. Log reads. Files saved with `Modified_effective_depth` and `mod_depth_log`, respectively.
5. `generate_cotrans_heatmap` should be run in terminal as follows:

   `python <path to generate_cotrans_heatmap.py> -i <input_matrix_csv> -v <variables_txt> -o <output_folder>`

## Comparing experimental replicates using plot_cotrans_correlation

`plot_cotrans_correlation` plots nucleotide reactivities of two replicates against one another on a hexbin plot to visualize correlation between replicates.

### generate_cotrans_heatmap inputs

Required: 

```
-i/--input_matrix <input_matrix_csv>    Path to reactivity matrix .csv file generated from running
                                        generate_cotrans_heatmap. Requires 2 files

-v/--variables <variables_txt>          Path to the variables file. This will implement different options
                                        during the run (see below)

-o/--output_folder <output_folder>      Path to the folder to save output file(s)

```

### Basic usage of plot_cotrans_correlation

1. `input_matrix_csv` filepaths should lead to .csv reactivity matrix files generated by `generate_cotrans_heatmap`. Ideally these two reactivity matrices are experimental replicates of the same riboswitch and conditions.
2. The `variables_txt` file provided to `plot_cotrans_correlation` should be filled in to the users specifications. This file is provided in `TECprobe_visualization/plot_cotrans_correlation` as `plot_cotrans_correlation_variables.txt`. Each variable that can be set contains additional information about that variable, accepted arguments for that variable, as well as the default value (if applicable).
3. Specific ranges of transcript lengths are used for comparison, as certain transcript lengths are not enriched due to experimental design (for more information, see Szyjka and Strobel (2023)). If the user is analyzing any of the three riboswitches studied in the paper (_Clostridium beijerinckii pfl_ ZTP, _Bacillus cereus crcB_ fluoride, or _Clostridiales bacterium oral taxon_ 876 str. F0540 ppGpp), `ribo_stored` can be set to `True`, otherwise `ribo_stored` should be set to `False`.
4. If `ribo_stored = False`, the user will have to designate the transcript lengths used for replicate comparison by setting the number of different transcript length ranges (`num_range`) and the min/max values of each range (`txpt_keep`). Further detail on these parameters can be read in `plot_cotrans_correlation_variables.txt`. **NOTE: Both `input_matrix_csv` files should be from the same riboswitch or at least compare the same transcript lengths**
5. Two files are generated when running `plot_cotrans_correlation`:
   * A .csv file that reports the nucleotide reactivity at each transcript length and nucleotide position for both replicates.
   * A hexbin plot generated from the data shown in the above .csv file, where the first `input_matrix` is recorded on the x-axis and the second `input_matrix` is recorded on the y-axis. The filetype is chosen by the user.
6. `plot_cotrans_correlation` should be run in terminal as follows:

   `python <path to plot_cotrans_correlation> -i <input_matrix_csv_1> <input_matrix_csv_2> -v <variables_txt> -o <output_folder>`

## Comparing neighboring transcripts using compare_cotrans_neighbor

`compare_cotrans_neighbor` calculates a Pearson Correlation Coefficient (r) when comparing nucleotide reactivity and background mutation rate between neighboring transcript lengths. These r values are then assembled into a violin for each dataset and comparison.

### compare_cotrans_neighbor inputs

```
-i/--input_matrix <input_matrix_csv>    Path to Full_Table.csv file(s) generated from running
                                        compile_sm2_output. Accepts 1+ file(s)

-v/--variables <variables_txt>          Path to the variables file. This will implement different options
                                        during the run (see below)

-o/--output_folder <output_folder>      Path to the folder to save output file(s)

```

### Basic usage of compare_cotrans_neighbor

1. `input_matrix` filepath(s) should lead to Full_Table.csv file(s) generated from `compile_sm2_output`.
2. The `variables_txt` file provided to `compare_cotrans_neighbor` should be filled in to the users specifications. This file is provided in `TECprobe_visualization/compare_cotrans_neighbor` as `compare_cotrans_neighnor_variables.txt`. Each variable that can be set contains additional information about that variable, accepted arguments for that variable, as well as the default value (if applicable).
3. Specific ranges of transcript lengths are used for comparison, as certain transcript lengths are not enriched due to experimental design (for more information, see Szyjka and Strobel (2023)). If the user is analyzing any of the three riboswitches studied in the paper (_Clostridium beijerinckii pfl_ ZTP, _Bacillus cereus crcB_ fluoride, or _Clostridiales bacterium oral taxon_ 876 str. F0540 ppGpp), `ribo_stored` can be set to `True`, otherwise `ribo_stored` should be set to `False`.
4. If `ribo_stored = False`, the user will have to designate the transcript lengths used for replicate comparison by setting the number of different transcript length ranges (`num_range`) and the min/max values of each range (`txpt_keep`). Further detail on these parameters can be read in `compare_cotrans_neighbor_variables.txt`. **NOTE: All `input_matrix_csv` files should be from the same riboswitch or at least compare the same transcript lengths**
5. Four files are generated when running `compare_cotrans_neighbor`:
   * Two violin plots: 1. r values of nucleotide reactivity (background-subtracted mutation rate) and 2. r values of background mutation rate. Each dataset provided will be provided on a single graph.
   * Two .csv files: 1. r values of nucleotide reactivity (background-subtracted mutation rate) saved with `Reactivity.csv` and 2. r values of background mutation rate saved with `BkgrndMutRate.csv`
6. `compare_cotrans_neighbor` should be run in terminal as follows:

   `python <path to compare_cotrans_neighbor.py> -i <input_matrix_csv> -v <variables_txt> -o <output_folder>`
