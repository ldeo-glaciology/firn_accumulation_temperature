
# Notes on including tetmerature in the compaction model
## introduction
I want to add temperature to the firn compaction model and publish the model and some numerical tests as a brif communication. This will allow Elizabeth to have a fully working, tested version of the model and code with which to start her third project of her thesis. She will start this when she returns from Antarctic in Jan 2023, or perhaps before she leaves in Nov 2022. 

I would like to write the code in python so that: 
- it will be accessible to anyone
- I can make a binder of the code
- I get practise with python
- we will have a python version of the original isothermal code

This document documents progress on making a (A) python version of the existing isothermal code, and (B) adding temperature.

## 7/22/22
I have a python script which setups up the model and integrates the equations with the Method Of Lines (MOL). The scipy ODE solver (solve_ivp) succesfully completes its integration, but numerical instabilities are clearly present in the solution. 

The solution is saved as an xarray. 

This is currently just a script , but I would like to turn this into a class, which class functions that:
- initalize the instance of the class
- setup the model, similar to FirnSetuoX.m in the matlab version
- integrates the model, similar to FCMY.m in the matlab version
- plots results
- does the initialization, setup and integration in one function.

First I need to get to the bottom of these numercial instabilities. 


### the first iteration of the solver
The approach is to start the first iteration of the solver and then run the equations functions in the each case (*eqns* in python and *GovEq* in matlan) and look at what is differnent. 

The initial conditions are identical between the matlab and python versions. 

dphidt is very different. Tracing this backwards through the code, w is not the same, and the root cause is that sigma is not the same. 

sigma is different because the intration functions are giving different outputs. The inputs to the integration functions are the same in each case (s_int) is identical.

python test
``` integrate.cumtrapz(z_h, np.ones_like(z_h), initial=0) = zeros```

matlab test 
``` cumtrapz(z_h, s_int*0+1) = 0:0.01:1```
DIFFERENT!

The issue was that the inputs to integrate.cumtrapz should be in the opposite order:

```integrate.cumtrapz(np.ones_like(z_h),z_h) = 0:01:1```
as expected.

I corrected this, and updated to the latest version of cumtrapz in scipy, which is called integrate.cumulative_trapezoid, in the line where sigmais calculated and everywhere else that cumulative integration is used. 

Sigma is now the same between the two version. 


The following are now the same across the two versions in the first iteration
- dphidt
- dr2dt 
- dAdt
- dTdt
- dHdt



The code now produces sensible looking results. E.g., a porosity profile which looks very similar to the matlab results.

The steady state porosity profile is the same between the two versions to within 0.1%. But I can quantify this more thoroughly.

The velocity, the Depth and the total mass, which are computed in post-processing, also appear to the be same between the two versions. 

Firn thickness: The matlab code only computes the final firn thickness. I have modifed the pythoncode to compute it at every time step. The firn thickness at the last time step matches between the two codes. 

The python code now computes firn air content as \int phi(z,t) dz

I currently have an issue with repmat and the orientation of the z_h vector. I solved this by using meshgrid for both parts of the Depth calculation instead of repmat. 


Next jobs are to
- commit recent changes (done)
- refactor the code so that it is class that can be loaded into a jupyter notebook (done)
- write a NB which runs the model and compare the results to the matlab results. 


## 7/23/22

fcm.py implements the code as a class. It saves the results as a class attribute, fcm.results. Three simple plotting routines have been added to make plots of profiles, time series, and z-t plots. 

Next jobs are to
- write a NB which runs the model and compare the results to the matlab results. (started)
- add comparison to steady state model to the NB and as a class function. 
- add class function that loops round a given list of parameters and repeatedly runs the model and saves the result in an xarray. 


I have made some progress on looping over arbitrary parameters in the NB ployying_figures.ipynb. One thing I dont understand is how xarray.concat is working when there is a variable number of time steps between different simulations.

A general comment is that the post processing seems to be taking longer than the integration in cases when beta is larger (and presumably there are lots of time steps because things evolve quickly). 


## 7/24/22
param_search_1D can now loop over prescribed parameters and saves the results in a concatenated xarray.dataset. 

## 7/25/22
fcm_demo.ipny is a start on a NB to demo the firn code

## 7/26/22
I think what I probably need to do is interpolate all the values onto a regular depth grid in the main model - i.e. present it in terms of z instead of z_h. 

