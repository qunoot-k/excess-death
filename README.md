# excess-death
In this project I will implement a demographic model for England and Wales, 
to predict the expected number of deaths per week from the beginning of 2020 
assuming death rates had stayed at the levels seen in 2017-19.
Data for annual probability of death for each one year age group, 
and the populations of each age group at the start of a period
in UK are available from the Office for National Statistics (ONS)
The model will include an adjustment to allow for seasonal variation 
in mortality rates. I will then compare this to the actual deaths per week 
that occurred over this period, to obtain an excess deaths time series, 
and then model this series using a simple Bayesian model in JAGS.
