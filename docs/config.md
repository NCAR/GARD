# Configuring GARD

After building GARD, it is run on the command line following this syntax:

```
./gard downscale_options.txt
```

where `downscale_options.txt` is a Fortran style namelist (see example in GARD/run).

The namelist options are described in the tables below:

## Parameters

| Name                            | Type    | Required? | Default              | Description                                                                      |
|---------------------------------|---------|-----------|----------------------|----------------------------------------------------------------------------------|
| training_file                   | string  | No        | Options File Name    | Used to specify a separate namelist file for the training section shown below    |
| observation_file                | string  | No        | Options File Name    | Used to specify a separate namelist file for the observation section shown below |
| prediction_file                 | string  | No        | Options File Name    | Used to specify a separate namelist file for the prediction section shown below  |
| output_file                     | string  | No        | downscaled_output.nc | Downscaled output filename.                                                      |
| n_analogs                       | integer | No (1)    | -1                   | The number of analogs to include                                                 |
| n_log_analogs                   | integer | No (1)    | -1                   |                                                                                  |
| analog_threshold                | real    | No (1)    | -1                   |                                                                                  |
| start_date                      | string  | Yes       | n/a                  | Start date for prediction period                                                 |
| end_date                        | string  | Yes       | n/a                  | End date for prediction period                                                   |
| start_train                     | string  | Yes       | n/a                  | Start date for training period                                                   |
| end_train                       | string  | Yes       | n/a                  | End date for training period                                                     |
| start_transform                 | string  | Yes       | n/a                  | Start date for transform period                                                  |
| end_transform                   | string  | Yes       | n/a                  | End date for transform period                                                    |
| pure_regression                 | logical | No        | FALSE                | Use pure regression downscaling approach                                         |
| pure_analog                     | logical | No        | FALSE                | Use pure analog downscaling approach                                             |
| analog_regression               | logical | No        | TRUE                 | Use analog regression downscaling approach                                       |
| sample_analog                   | logical | No        | FALSE                |                                                                                  |
| logistic_from_analog_exceedance | logical | No        | FALSE                |                                                                                  |
| logistic_threshold              | logical | No        | -9999                |                                                                                  |
| weight_analogs                  | logical | No        | TRUE                 |                                                                                  |
| debug                           | logical | No        | TRUE                 |                                                                                  |

Notes:

1.  Must specify one of n_analogs, n_log_analogs, or analog_threshold

## Training Parameters

| Name                  | Type    | Required? | Default | Description                                                                              |
|-----------------------|---------|-----------|---------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a     | Name of the training parameters dataset                                                  |
| preloaded             | string  | Yes       | n/a     | filepath of the preloaded training parameters dataset                                    |
| interpolation_method  | integer | No        | 1       | nearest neighbor= 1, bilenear =2                                                         |
| time_indices          | integer | Yes       | -1      |                                                                                          |
| nvars                 | integer | Yes       | -1      | number of variables to be used in training                                               |
| data_type             | string  | Yes       | n/a     |                                                                                          |
| lat_name              | string  | Yes       | n/a     |                                                                                          |
| lon_name              | string  | Yes       | n/a     |                                                                                          |
| time_name             | string  | Yes       | n/a     |                                                                                          |
| nfiles                | integer | Yes       | -1      |                                                                                          |
| input_transformations | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a     | variables names to use in training                                                       |
| file_list             | string  | Yes       | n/a     |                                                                                          |
| selected_time         | integer | Yes       | -1      |                                                                                          |
| calendar              | string  | Yes       | n/a     |                                                                                          |
| calendar_start_year   | integer | No        | 1900    |                                                                                          |
| timezone_offset       | real    | No        | 0       | hours to offset time variable (e.g. `time_name`) to account for timezone.                |


## Prediction Parameters

| Name                  | Type    | Required? | Default | Description                                                                              |
|-----------------------|---------|-----------|---------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a     |                                                                                          |
| preloaded             | string  | Yes       | n/a     |                                                                                          |
| interpolation_method  | integer | No        | 1       | nearest neighbor= 1, bilenear =2                                                         |
| normalization_method  | integer | No        | 0       | mean/stddev from: prediction data = 0, training data = 1                                 |
| nvars                 | integer | Yes       | -1      |                                                                                          |
| data_type             | string  | Yes       | n/a     |                                                                                          |
| lat_name              | string  | Yes       | n/a     |                                                                                          |
| lon_name              | string  | Yes       | n/a     |                                                                                          |
| time_name             | string  | Yes       | n/a     |                                                                                          |
| nfiles                | integer | Yes       | -1      |                                                                                          |
| transformations       | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| input_transformations | integer | No        | 0       | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a     |                                                                                          |
| file_list             | string  | Yes       | n/a     |                                                                                          |
| selected_time         | integer | Yes       | -1      |                                                                                          |
| calendar              | string  | Yes       | n/a     |                                                                                          |
| calendar_start_year   | integer | No        | 1900    |                                                                                          |
| timezone_offset       | real    | No        | 0       | hours to offset time variable (e.g. `time_name`) to account for timezone.                |

## Observation Parameters

| Name                  | Type    | Required? | Default       | Description                                                                              |
|-----------------------|---------|-----------|---------------|------------------------------------------------------------------------------------------|
| name                  | string  | Yes       | n/a           |                                                                                          |
| preloaded             | string  | Yes       | n/a           |                                                                                          |
| nvars                 | integer | Yes       | -1            |                                                                                          |
| nfiles                | integer | Yes       | -1            |                                                                                          |
| data_type             | string  | Yes       | n/a           |                                                                                          |
| lat_name              | string  | Yes       | n/a           |                                                                                          |
| lon_name              | string  | Yes       | n/a           |                                                                                          |
| time_name             | string  | Yes       | n/a           |                                                                                          |
| input_transformations | integer | No        | 0             | no transform = 0, quantile mapping = 1, log transform = 2, cube root = 3, fifth root = 4 |
| var_names             | string  | Yes       | n/a           |                                                                                          |
| file_list             | string  | Yes       | n/a           |                                                                                          |
| calendar              | string  | Yes       | n/a           |                                                                                          |
| calendar_start_year   | integer | No        | 1900          |                                                                                          |
