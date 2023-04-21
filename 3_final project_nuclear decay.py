# -*- coding: utf-8 -*-
"""
@author: Pablo Simarro Gomez 
Created on Mon Nov 25 22:13:45 2019

PHYS20161 3rd project: Nuclear Decay

This code estimates the decay constants and decay time of 79Rb and 79Sr from the 79Rb 
activity data measured as a sample of 79Sr decays to 79Rb and then to 79Kr.
The code validates the data, deleting data with errors in both numerical and physical values.
It also looks for outliers, eliminates them and shows them in a plot.
By performing a two-parameter minimised chi squared fit on the data, the code finds the 
decay constants and decay time, and calculates the errors using a contour plot.
It also calculates the reduced chi squared value to two decimal places and produces
a plot of the fit against the experimental data.
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.optimize import fmin
import sys

"""parameter definitions"""
lambda_Rb_guess = 0.001
lambda_Sr_guess = 0.01
Initial_N_Sr = (1e-6*const.Avogadro)
data = np.zeros((0,3))
outliers = []
outlier_values = np.zeros((0,3))
errors = 0
new_fit = 0

"""function definitions"""
def expected_Rb_activity(t,lambda_Rb,lambda_Sr):
    """
    Returns the expected activity of 79Rb for the given values:
    t (float)
    lambda_Rb (float)
    lambda_Sr (float)
    """
    return lambda_Rb*Initial_N_Sr*((lambda_Sr/(lambda_Rb - lambda_Sr))*(np.exp(-lambda_Sr*t) - np.exp(-lambda_Rb*t)))

def chi_squared(a,b,data):
    """
    Returns chi squared for the expected_Rb_activity function, depenedent on one
    variable, t, with two parameters, lambda_Rb & lambda_Sr.
    Data type:
    data array([float, float, float])
    a (float)
    b (float)
    """
    chi_squared = 0
    for line in data:
        chi_squared += ((expected_Rb_activity(line[0], a, b) - line[1])/line[2])**2
    return chi_squared

def function(lambda_values):
    """
    Requires input array like [x,y]
    Returns chi_squared of parameters given in lambda_values
    lambda_values = [(float), (float)]
    """
    lambda_Rb = lambda_values[0]
    lambda_Sr = lambda_values[1]
    return chi_squared(lambda_Rb,lambda_Sr,data)

"""main code"""
"""
Read and combine data using np.genfromtxt
If files are not found, exits the program and warns the user
"""
try:
    data_1 = np.genfromtxt('Nuclear_data_1.csv', comments = '%', delimiter = ',')
    data_2 = np.genfromtxt('Nuclear_data_2.csv', comments = '%', delimiter = ',')
    data_combined = np.concatenate((data_1,data_2),axis=0)
except:
    print("Files not found. Nuclear_data_1 and Nuclear_data_2 must be saved in the same directory as this code's file.")
    sys.exit()
    
"""extra validation"""  
for line in data_combined:
    if np.isnan(line[0]) or np.isnan(line[1]) or np.isnan(line[2]):
        errors += 1
    elif line[2] <= 0 or line[1] < 0 or line[0] < 0:
        errors += 1
    else:
        """change values to desired units"""
        line[0] *= 3600
        line[1] *= 1e12
        line[2] *= 1e12
        data = np.vstack((data,line))
"""sort data by time"""
data = data[data[:,0].argsort()]
    
"""
Minimnised chi squared fit for two parameters
"""
fit_data = fmin(function,(lambda_Rb_guess, lambda_Sr_guess), full_output = True, disp = 0)
[lambda_Rb_min, lambda_Sr_min] = fit_data[0]
"""
Check for outliers in data. 
If any outliers are found, remove them from data and perform fit again
"""
mean = np.mean(data[:,1] - expected_Rb_activity(data[:,0],lambda_Rb_min,lambda_Sr_min))
three_std = 3*(np.std(data[:,1] - expected_Rb_activity(data[:,0],lambda_Rb_min,lambda_Sr_min)))
for i in range(len(data)):
    if abs(data[i,1] - expected_Rb_activity(data[i,0],lambda_Rb_min,lambda_Sr_min)) > (mean + three_std):
        outliers.append(i)
        new_fit = 1
        
if new_fit == 1:
    outlier_values = np.vstack((data[outliers])) #used for plot later
    data = np.delete(data, outliers, 0)
    fit_data = fmin(function,(lambda_Rb_guess, lambda_Sr_guess), full_output = True, disp = 0)
    [lambda_Rb_min, lambda_Sr_min] = fit_data[0]

chi_squared_min = fit_data[1]
reduced_chi_squared = chi_squared_min/(len(data)-2)
half_life_Rb = np.log(2)/(lambda_Rb_min*const.minute)
half_life_Sr = np.log(2)/(lambda_Sr_min*const.minute)

"""
outliers plot
"""
outliers_fig = plt.figure()
outliers_plot = outliers_fig.add_subplot(111)
outliers_plot.set_title(r'Scatter plot of fitted data (blue) and the {0} outliers found (red).'.format(len(outliers)), fontsize = 10)
outliers_plot.set_xlabel('79Rb Activity (counts)', fontsize = 12)
outliers_plot.set_ylabel('Time (seconds)', fontsize = 12)

outliers_plot.errorbar(data[:,0], data[:,1], yerr = data[:,2], fmt = 'o', color = 'blue')
outliers_plot.errorbar(outlier_values[:,0], outlier_values[:,1], yerr = outlier_values[:,2], fmt = 'o', color = 'r')
# plt.savefig('Outliers_in_scatter.png', dpi = 300)
plt.show()

"""
fit plot
"""
fitted_data_figure = plt.figure()
fitted_data_plot = fitted_data_figure.add_subplot(111)
fitted_data_plot.set_title(r'Minimised $\chi^2$ fit of the scatter data.\n Fit values: lambda 79Rb = {0:#.3g}, lambda 79Sr = {1:#.3g}, Reduced $\chi^2$ = {2:3.2f}'.format(lambda_Rb_min,lambda_Sr_min,reduced_chi_squared))
fitted_data_plot.set_xlabel('79Rb Activity (counts)', fontsize = 12)
fitted_data_plot.set_ylabel('Time (seconds)', fontsize = 12)

fitted_data_plot.plot(data[:,0], expected_Rb_activity(data[:,0], lambda_Rb_min, lambda_Sr_min), color ='r')
fitted_data_plot.errorbar(data[:,0], data[:,1], yerr = data[:,2], fmt = 'o', color = 'blue')

# plt.savefig('Min_chi_squared_fit_plot.png', dpi = 300)
plt.show()

"""
contours plot
"""
Rb_array = np.linspace(lambda_Rb_min - lambda_Rb_min*0.1, lambda_Rb_min + lambda_Rb_min*0.1 , 300)
Sr_array = np.linspace(lambda_Sr_min - lambda_Sr_min*0.1 , lambda_Sr_min + lambda_Sr_min*0.1, 300)
lambda_Rb_mesh , lambda_Sr_mesh = np.meshgrid(Rb_array,Sr_array)

contours_figure = plt.figure()
contours_plot = contours_figure.add_subplot(111)

contours_plot.set_title(r'$\chi^2$ contours against lambda 79Rb and lambda 79Sr.', fontsize = 14)
contours_plot.set_xlabel('lambda 79Rb (s^-1)', fontsize = 12)
contours_plot.set_ylabel('lambda 79Sr (s^-1)', fontsize = 12)
contours_plot.set_xlim((lambda_Rb_min - lambda_Rb_min*0.04, lambda_Rb_min + lambda_Rb_min*0.04))
contours_plot.set_ylim((lambda_Sr_min - lambda_Sr_min*0.1, lambda_Sr_min + lambda_Sr_min*0.1))
contours_plot.grid(axis='y', linestyle='-')
contours_plot.grid(axis='x', linestyle='-')

contours_plot.scatter(lambda_Rb_min, lambda_Sr_min)
chi_squared_contours = (chi_squared_min + 1.00,chi_squared_min + 2.30, chi_squared_min + 5.99)
contours = contours_plot.contour(lambda_Rb_mesh, lambda_Sr_mesh, chi_squared(lambda_Rb_mesh,lambda_Sr_mesh,data),levels = chi_squared_contours)

"""add labels"""
labels = [r'$\chi^2_{{\mathrm{{min.}}}} = $ {0:4.1f}'.format(chi_squared_min),r'$\chi^2_{{\mathrm{{min.}}}}+1.00$', r'$\chi^2_{{\mathrm{{min.}}}}+2.30$',r'$\chi^2_{{\mathrm{{min.}}}}+5.99$']
for i in range(len(labels)):
    contours_plot.collections[i].set_label(labels[i])

"""adjust scale of plot"""   
box = contours_plot.get_position()
contours_plot.set_position([box.x0, box.y0, box.width*0.7, box.height])
contours_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5),fontsize = 12)
"""get contour data for errors"""
data_contours = contours.allsegs[0][0]

# plt.savefig('Chi_squared_contours_plot.png', dpi = 300)
plt.show()

"""
error calculation using contour plot
"""
error_lambda_Rb = max(data_contours[:,0]) - lambda_Rb_min
error_lambda_Sr = max(data_contours[:,1]) - lambda_Sr_min
error_half_life_Rb = np.log(2)*(error_lambda_Rb/lambda_Rb_min)*half_life_Rb
error_half_life_Sr = np.log(2)*(error_lambda_Sr/lambda_Sr_min)*half_life_Sr

"""
print results
"""
print('{0} data points with errors were found and erased. \n{1} outliers were also found and erased from the data provided (shown in Outliers_in_scatter plot)'.format(errors,len(outliers)))
print(' ')
print('Reduced chi squared of the fit = {0:3.2f}'.format(reduced_chi_squared))
print(' ')
print('The decay constant of 79Rb is ({0:#.3g} +- {1:7.6f}) s^(-1)'.format(lambda_Rb_min,error_lambda_Rb))
print('The decay constant of 79Sr is ({0:#.3g} +- {1:6.5f}) s^(-1)'.format(lambda_Sr_min,error_lambda_Sr))
print(' ')
print('The half life of 79Rb is ({0:#.3g} +- {1:3.1f}) minutes'.format(half_life_Rb,error_half_life_Rb))
print('The half life of 79Sr is ({0:#.3g} +- {1:3.2f}) minutes'.format(half_life_Sr,error_half_life_Sr))