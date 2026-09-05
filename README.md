# india-power-attribution

This is the code used for my MSc Energy Systems master's dissertation on the impact of climate change on the India power system during the 2022 India-Pakistan heatwave. The data directory is available on [OneDrive](https://unioxfordnexus-my.sharepoint.com/:f:/r/personal/kell8304_ox_ac_uk/Documents/dissertation_data?d=w0e47df9c10a04a10bc78bda05fc5f79c&csf=1&web=1&e=fXcSGI) for Nexus365 users. Ensure that it is downloaded and placed in the base directory before running any code.

1. Prepare the data by running build_outage_panel.ipynb and build_demand_shortage_panel.ipynb.

2. Fit the models and perform attribution by running fit_outage_models.R and fit_demand_shortage_models.R.

3. Visualise the results in create_figures.ipynb.

Note that climate data for attribution has already been processed and is available in data/attribution_runs.


