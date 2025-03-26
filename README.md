# laex_detect_lines_isis

This repository contains a set of functions designed for detecting and fitting Gaussian emission lines in x-ray spectra with ISIS. The code provides functions to define parameters, evaluate models, and fit Gaussian lines to the data based on a chi-square minimization approach. This function operates in keV, and requires ISISSCIPTS and Python to be properly set up in your environment.

## Functions Overview

### Public Functions

#### `create_line_model(name)`
- **Purpose**: Generates a line model using Python scripts and saves it under a specified name.
- **Key Elements**:
  - `write_plot("spec")`: Writes the plot to a file.
  - `system(cmd)`: Executes the Python script `lines.py` to generate a model.
  - `evalfile("set_line_model.sl")`: Executes a secondary script to set the line model.
- **Output**: The lines are saved in the model, and a message is displayed indicating that the lines have been saved. This model, named **linemodel**, is ready to be used within other model components.

#### `fit_line_model(evalfun)`
- **Purpose**: Fits a line model to the data by iterating through multiple candidate Gaussian models, adjusting the parameters to minimize the chi-square statistic.
- **Key Elements**:
  - `how_many_gaussian()`: Determines how many Gaussian components are required for the model.
  - `get_params()`: Retrieves the parameters for fitting the model.
  - `test_params(p, evalfun)`: Evaluates the chi-square statistic for the current set of parameters.
  - `thaw` and `freeze`: These functions control whether certain parameters are allowed to vary or are fixed during the fitting process.
- **Returns**: The function sets the parameters of the best-fitting model after iterating and finding the model with the lowest chi-square statistic. Only emission lines which improve the χ² statistic by 0.005 are kept as reliable candidates.

## Workflow

1. **`create_line_model(name)`** generates and saves an initial line model using Python scripts.
2. **`fit_line_model(evalfun)`** iterates through potential models, adjusting parameters (such as area, sigma, and center of Gaussians), and selects the best model based on the reduced chi-square statistic.

The functions make use of system commands, Python scripts, and chi-square evaluations to fit the best Gaussian line model to the data.

## Use example with ISIS:

### Detect line example:

1. We should save the files **create_line_model.sl** and **lines.py** in the same folder, for example, under the desired name (in this example, the folder is named **detect_lines_isis**) in our local models directory and indicate the path when we enter ISIS:

    ```isis
    isis> require("isisscripts.sl")
    isis> variable MODEL_DIR;
    isis> MODEL_DIR = "/our_local_models_path/detect_lines_isis/";
    isis> require(MODEL_DIR + "/create_line_model.sl");
    ```

2. As we would normally do, we load our spectra and make a plot, rebinning it to our preference:

    ```isis
    isis> load_data("spec.ds");
    isis> group(1;sn=12);
    isis> plot_data(1);
    ```

3. We can now create a model containing our emission lines, named **linemodel**, ready to use with other components as required by our spectra. Apart from this model, if we provide a name, we will save this line model in a file named `set_line_model_+name+.sl` in our folder to recover the complete model in subsequent runs:

    ```isis
    isis> create_line_model("test");
    ```

4. We define our complete model:

    ```isis
    isis> fit_fun("(tbnew(1)+constant(1)*tbnew(2))*(powerlaw(1)+bbody(1)+bbody(2)+linemodel)");
    ```

5. We fit our lines. Several line candidates are suggested, but only those which improve the fit by an absolute value of χ² = 0.005 are kept as free parameters, with an area different than 0:

    ```isis
    isis> fit_line_model(&eval_counts);
    ```

6. Once we are satisfied, we save our model as usual:

    ```isis
    isis> save_par("example_model.par");
    ```

### Recover our model example:

1. We should enter ISIS as usual, and first, we should require the saved **linemodel**:

    ```isis
    isis> require("set_line_model_test.sl");
    ```

2. We now load and recover the complete model:

    ```isis
    isis> load_par("example_model.par");
    ```

In case that we load them in the inverse order, the message "linemodel is undefined" will appear. We then must require our line model and load again.

## License

MIT license
