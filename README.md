# nuclear_decay_analysis

Given two csv files containing the 79Rb activity data measured as a sample of 
79Sr decays to 79Rb and then to 79Kr estimates the decay constants and decay time of 79Rb and 79Sr.

The code validates the data, deleting data with errors in both numerical and physical values.
It also looks for outliers, eliminates them and shows them in a plot.
By performing a two-parameter minimised chi squared fit on the data, the code finds the 
decay constants and decay time, and calculates the errors using a contour plot, also shown.
It also calculates the reduced chi squared value to two decimal places and produces
a plot of the fit against the experimental data.

```console
>> python "3_final project_nuclear decay.py"
4 data points with errors were found and erased.
3 outliers were also found and erased from the data provided (shown in Outliers_in_scatter plot)

Reduced chi squared of the fit = 1.16

The decay constant of 79Rb is (0.000510 +- 0.000004) s^(-1)
The decay constant of 79Sr is (0.00508 +- 0.00013) s^(-1)

The half life of 79Rb is (22.6 +- 0.1) minutes
The half life of 79Sr is (2.28 +- 0.04) minutes
```
