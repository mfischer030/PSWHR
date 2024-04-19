# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 21:59:36 2024

@author: fism
"""

import pvlib

from pvlib.modelchain import ModelChain
from pvlib.location import Location
from pvlib.pvsystem import PVSystem
from pvlib.temperature import TEMPERATURE_MODEL_PARAMETERS

import matplotlib.pyplot as plt
import pandas as pd

location = Location(latitude=46.89193533578202, 
                    longitude=7.5428324192145295, 
                    tz="Europe/Zurich", 
                    altitude=549+5, 
                    name='Tesla Supercharger Station Rubigen')

sandia_modules = pvlib.pvsystem.retrieve_sam('SandiaMod')
cec_inverters  = pvlib.pvsystem.retrieve_sam('CECInverter')

temperature_parameters = TEMPERATURE_MODEL_PARAMETERS['sapm']['open_rack_glass_glass']

module   = sandia_modules['Canadian_Solar_CS5P_220M___2009_']
# inverter = cec_inverters['ABB__PVI_3_0_OUTD_S_US__208V_']
inverter = cec_inverters['ABB__ULTRA_750_TL_OUTD_1_US_690_x_y_z__690V_']

system = PVSystem(surface_tilt                 = 45,         # installed angle of pv panels
                  surface_azimuth              = 180,        # orientation of the pv system
                  module_parameters            = module,
                  inverter_parameters          = inverter,
                  temperature_model_parameters = temperature_parameters,
                  modules_per_string           = 150,
                  strings_per_inverter         = 20
                  )

modelchain = ModelChain(system, location) 

times = pd.date_range(start="2019-01-01", 
                       end="2019-12-31",
                      # end="2019-01-31",
                      # end="2019-01-08",
                      freq="h",
                      tz= location.tz
                      )

clear_sky = location.get_clearsky(times)

# clear_sky.plot(figsize=(16,9)) #ghi = global horizontal irradiance; dni = direct normal irradiance; dhi = diffuse horizontal irradiance
plt.show

modelchain.run_model(clear_sky)
modelchain.results.ac.plot(figsize=(16,9))
plt.show

Area_PV = system.modules_per_string * module.Area
print(Area_PV,"squaremeters")