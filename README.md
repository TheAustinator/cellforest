<p align="center"><img width=80% src="https://raw.githubusercontent.com/theaustinator/cellforest/master/static/cellforest.jpg" alt="cellforest"></p>

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![Python](https://img.shields.io/badge/python-v3.6+-blue.svg)
![Contributions welcome](https://img.shields.io/badge/contributions-welcome-orange.svg)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/python/black)
[![GitHub Issues](https://img.shields.io/github/issues/TheAustinator/cellforest.svg)](https://github.com/TheAustinator/cellforest/issues)
[![Build Status](https://travis-ci.org/TheAustinator/cellforest.svg?branch=master)](https://travis-ci.org/TheAustinator/cellforest)
<!--[![PyPI version](https://badge.fury.io/py/.svg)](https://badge.fury.io/py/cellforest)-->


## A simple, interactive, and customizable single cell workflow manager
<p align="center">
  <a href="#core concepts">Core Concepts</a> •
  <a href="#overview">Overview</a> •
  <a href="#features">Features</a> •
  <a href="#usage">Usage</a> •
  <a href="#upcoming-features">Upcoming Features</a>
</p>


## Core Concepts



### Features

## Usage
**Install**
```
pip install cellforest
```
**Import**
```python
from cellforest import CellBranch
```

### Examples

## Upcoming Features

## Automated quality control (QC) plotting 

### I. Setting up QC plotting

1. All plots are defined in `plot_map` dictionary in the config file (`cellforest/config/default_config.yaml`)

2. A typical `plot_map` dictionary looks like this in a config file:

   ```yaml
   plot_map:
     root:
       _UMIS_PER_CELL_HIST_: ~
       _UMIS_PER_BARCODE_RANK_CURV_: ~
     normalize:
       _GENES_PER_CELL_HIST_:
         plot_method: plot_umis_per_cell_hist
         filename: gene_per_cell_hist
       _UMIS_PER_CELL_HIST_: ~
   ```

   1. At the top level, process name is defined (e.g., `root`, `normalize`, `cluster`) — these would indicate definition of plots per process run
   2. Inside each key is a dictionary of plot names. The format used in CellForest is all-caps bound with underscores (e.g., `_THIS_IS_A_PLOT_NAME_`)
   3. For each plot, you can define `plot_method` or `filename`. `plot_method` will be scraped from scripts in `cellforest/plot` and map to the plotting function
   4. `filename` refers to the name of the image that will be saved upon running the plotting function. Note that you can omit the extension for the filename but you can also include it to force a specific extension for each plot
   5. You can use `~` as value for a plot key to indicate usage of default parameters (we will talk about parametrization in the next section), such as inferred filename (`_THIS_IS_A_PLOT_NAME_` becomes `this_is_plot_name`) as well as inferred plot method (`_THIS_IS_A_PLOT_NAME_` becomes `plot_this_is_a_plot_name()`)

3. Plots are run with the help of the `hook_generate_plots` hook in the config file. Basing from the `plot_map`, after a process successfully runs (e.g., after `branch.process.normalize()`), mapped plotting functions will be triggered

### II. Parametrizing QC plots

1. To parametrize plots as well as the plotting methods that produce the plots, we can use the `plot_kwargs` field in each plot along with the required `plot_kwargs_defaults` in the config

2. Here is an example in a config file:

   ```yaml
   plot_kwargs_defaults:
     stratify:
       default: None
     plot_size:  # width x height in px
       default: [800, 800]
       large: [1600, 1600]
     filename_ext: png
     
   plot_map:
     normalize:
       _GENES_PER_CELL_HIST_:
         plot_kwargs:
           stratify: 
             - sample
             - default
           plot_size: large
   ```

   1. `plot_kwargs_defaults` defines the named values for each keyword argument option (e.g., `plot_size["default"] == [800, 800]`); for each keyword argument, there must a `default` option because this will be used for filename inference (we will talk about it in a bit)

   2. Keyword arguments in individual plots will try to seek for the value in `plot_kwargs_defaults` and if found, will use that value to pass on to the plotting function. If not, the noted value is passed on. For example, in the above case, `default` will map to `None` whereas `sample` will just pass on as `sample`

   3. **Important rule**: “same length or length 1.” In the above case, `stratify` has 2 defined options (is a list of 2 elements) so CellForest will infer that for each of those stratification options, the plot size will be large (`[1600, 1600]`).

   4. If we want to have one of them large, one of them default, we can write in the same order:

      ```yaml
      stratify: 
        - sample
        - default
      plot_size:
      	- large
      	- default
      ```

      Note: it will not work if you try to match keyword arguments of different lengths (e.g., 2 and 3) because there will not be a 1-to-1 mapping of arguments

   5. If we want both of them default size, we can just omit `plot_size`, and default value from `plot_kwargs_defaults` will be used

3. The  “same length or length 1” applies similarly to parameters of a plot. For example, if there a total of **two** plots because `stratify` has 2 options and `plot_size` has 1 or 2 options, and you want to define inidividual file names for each, there shall be 3 file names:

   ```yaml
   plot_map:
     normalize:
       _GENES_PER_CELL_HIST_:
         filename:
           - gene_per_cell_hist-a
           - gene_per_cell_hist-b
         plot_kwargs:
           stratify: 
             - sample
             - default
           plot_size: large
   ```

   1. In this case, you can let the `filename` be inferred or have to match the size of the `plot_kwargs` so that plot output files don’t override each other
   2. Same applies to `plot_method` if you happen to have multiple plotting methods

4. Sounds like a lot of parametrization? Don’t sweat because CellForest is smart and can infer the filenames and plot methods from the keyword arguments alone!

   ```yaml
   plot_kwargs_defaults:
     stratify:
       default: None
     plot_size:  # width x height in px
       default: [800, 800]
       large: [1600, 1600]
     filename_ext: png
   
   # both of the following plot maps are equivalent
   plot_map:
     normalize:
       _GENES_PER_CELL_HIST_:
         plot_kwargs:
           stratify: 
             - sample
             - default
       _UMIS_PER_CELL_HIST_: ~
           
   plot_map:
     normalize:
       _GENES_PER_CELL_HIST_:
       	plot_method: plot_genes_per_cell_hist
       	filename:
       		- genes_per_cell_hist-plot_size:800+800-stratify:sample
       		- genes_per_cell_hist-plot_size:800+800-stratify:none
         plot_kwargs:
           stratify: 
             - sample
             - default
           plot_size:
           	- default
           	- default
       _UMIS_PER_CELL_HIST_:
       	plot_method: plot_genes_per_cell_hist
       	filename: umis_per_cell_hist-plot_size:800+800-stratify:none
         plot_kwargs:
           stratify: default
           plot_size: default
   ```

   1. Did you catch it? There are multiple things that happened here! First, plot method and part of the filename is, of course, inferred from the plot name (`_GENES_PER_CELL_HIST`)
   2. Secondly, each filename has a unique suffix attached to it. There are two steps here:
      1. Keyword arguments are sorted alphabetically and stitched together with a `-`
      2. Values of these arguments are taken from `plot_kwargs_defaults` if available
      3. Everything is lowercase
   3. So now, if you have 10 different columns to stratify your plots by, you wouldn’t have to define 10 different filenames or plot methods — CellForest can infer them all

5. You can also add your own keyword arguments that the plotting function can use:

   1. HAVEN”T FINISHED HERE

6. (Optional) You can think about how this works in that `plot_map` under the hoods is parsed into three pieces:

   1. `plot_methods` that shows the mapping of plots to plotting functions

      ```python
      {
          "root": {
              "_UMIS_PER_CELL_HIST_": "plot_umis_per_cell_hist",
              "_UMIS_PER_BARCODE_RANK_CURV_": "plot_umis_per_barcode_rank_curv",
          },
          "normalize": {
              "_GENES_PER_CELL_HIST_": "plot_genes_per_cell_hist",
              "_UMIS_PER_CELL_HIST_": "some_absurd_name",
          }
      }
      ```

   2. `plot_map` (similar to `path_map`) that shows mapping of requested parameters to plot filenames

      ```python
      {
          "root": {
              "_UMIS_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "default"}': "umis_per_cell_hist-plot_size:800+800-stratify:none.png"
              },
              "_UMIS_PER_BARCODE_RANK_CURV_": {
                  '{"plot_size": "default", "stratify": "default"}': "umis_per_barcode_rank_curv-plot_size:800+800-stratify:none.png"
              },
          },
          "normalize": {
              "_GENES_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "default"}': "genes_per_cell_hist-plot_size:800+800-stratify:none.png",
                  '{"plot_size": "default", "stratify": "sample"}': "genes_per_cell_hist-plot_size:800+800-stratify:sample.png"
              },
              "_UMIS_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "lane"}': "umis_per_cell_hist-plot_size:500+500-stratify:none",
              }
      }
      ```

   3. `plot_kwargs` that shows the mapping of request parameters to parameters actually fed into plotting functions

      ```python
      {
          "root": {
              "_UMIS_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "default"}': {"stratify": "none", "plot_size": [800, 800]}
              },
              "_UMIS_PER_BARCODE_RANK_CURV_": {
                  '{"plot_size": "default", "stratify": "default"}': {"stratify": "none", "plot_size": [800, 800]}
              },
          },
          "normalize": {
              "_GENES_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "default"}': {"stratify": "none", "plot_size": [800, 800]}
                  '{"plot_size": "default", "stratify": "sample"}': {"stratify": "sample", "plot_size": [800, 800]}
              },
              "_UMIS_PER_CELL_HIST_": {
                  '{"plot_size": "default", "stratify": "default"}': {"stratify": "none", "plot_size": [800, 800]}
              }
      }
      ```

   4. Notice the overlap of keyword arguments “keys” in `plot_map` and `plot_kwargs`

### III. Adding new/customizing QC plots

1. Let’s talk a bit about how the defined parameters in `plot_map` are fed into the plotting functions:
   1. `hook_generate_plots` in DataForest checks the files in `cellforest/plot` folder and collects the plotting functions
   2. Defined plot method per plot then triggers the call of the respective plotting function
   3. Each function is decorated by either `qc_plot_py` or `qc_plot_r`. These decorators add the functionality of setting up the plot size as well as setting up the plotting functions for stratification. The rest of the keyword arguments