#################################################################################################################################################
# Basic StreamPULSE data processing pipeline                                                                                                    #
# Contact Mike Vlah (michael.vlah@duke.edu) with questions or comments.                                                                         #
#################################################################################################################################################

#####----- Install StreamPULSE pipeline tools from GitHub -----#####

remotes::install_github('streampulse/StreamPULSE', dependencies=TRUE)


#####----- Install latest version of streamMetabolizer -----#####

remotes::install_github("appling/unitted")
remotes::install_github("USGS-R/streamMetabolizer")


######----- Load packages -----#####

library(StreamPULSE)
library(streamMetabolizer)


#################################################################################################################################################
# View all variables at a site, and full available time range for that site. Note that USGS depth and discharge data may be available for sites # # that have associated USGS gage IDs, even if depth and discharge do not appear among the variables returned here. If USGS data are available,  # # they will be acquired automatically when you use prep_metabolism below. Likewise, air pressure and PAR estimates will be automatically        #
# acquired below, if necessary.                                                                                                                 #
#################################################################################################################################################

query_available_data(region = 'ID', site = 'Notus')
query_available_data(region = 'ID', site = 'Iverson')


#####----- Select site and date range for which to acquire StreamPULSE data. site_code is a combination of regionID and siteID -----#####

site_code = 'ID_Notus'
site_code_2 = 'ID_Iverson'
start_date = '2022-04-29'
end_date = '2022-10-20'


#####----- Download data from streampulse. token can be found by clicking the gear icon in the upper right of data.streampulse.org -----#####

sp_data = request_data(sitecode=site_code,
                       startdate=start_date, enddate=end_date)

sp_data_2 = request_data(sitecode=site_code_2,
                         startdate=start_date, enddate=end_date)


#####----- Choose model type for streamMetabolizer. Only "bayes" is available at this time. -----#####

model_type = 'bayes'


#####----- Which modeling framework to use. Use "streamMetabolizer" (the default); "BASE" is not available at this time. -----#####

model_name = 'streamMetabolizer'


#####----- Format data for metabolism modeling. -----#####

sp_data_prepped = prep_metabolism(d=sp_data, type=model_type,
                                  model=model_name, retrieve_air_pres=TRUE, estimate_PAR = FALSE)

sp_data_prepped_2 = prep_metabolism(d=sp_data_2, type=model_type,
                                    model=model_name, retrieve_air_pres=TRUE, estimate_PAR = FALSE)


#####----- Fit metabolism model and generate predictions (calls streamMetabolizer functions: mm_name, specs, metab, predict_metab). -----#####

model_fit = fit_metabolism(sp_data_prepped)
model_fit_2 = fit_metabolism(sp_data_prepped_2)


#################################################################################################################################################
# Plot results and diagnostics (This behaves unpredictably on some machines. If you sent your results to the data portal in the step above, you # # can view them more robustly there. If your data do not all occur within the same calendar year, the visualizations may still not work.)       #
#################################################################################################################################################

plot_output(model_fit)
plot_output(model_fit_2)


#####----- Here's where results and diagnostics live on data portal: http://data.streampulse.org:3838/streampulse_diagnostic_plots/ -----#####