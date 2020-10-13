import constants as c
import numpy
from pfp_utils import SeriestoMA, MAtoSeries

def absolutehumidityfromrelativehumidity(Ta,RH):
    # convert to masked arrays
    RH, WasND = SeriestoMA(RH)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    vp = RH * VPsat(Ta) / float(100)
    ah = float(1000000) * vp / ((Ta + 273.15) * c.Rv)
    # convert back to ndarray if input is not a masked array
    if WasND: ah, _ = MAtoSeries(ah)
    return ah

def co2_ppmfrommgCO2pm3(c_mgpm3,T,p):
    """
     Convert CO2 concentration units of mgCO2/m3 to umol/mol (ppm)
        Usage:
         CO2_ppm = co2_ppmfrommgCO2pm3(CO2_mgpm3, T, p)
         where
         CO2_mgpm3 (input) - CO2 concentration, mgCO2/m3
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the CO2 concentration in ppm.
    """
    # convert to masked array if required
    c_mgpm3, WasND = SeriestoMA(c_mgpm3)
    T, dummy = SeriestoMA(T)
    T = T + 273.15             # temperature in K
    p, dummy = SeriestoMA(p)
    p = p * float(1000)        # pressure in Pa
    # do the job
    c_ppm = (c_mgpm3/c.Mco2)*c.R*T/p
    # convert back to ndarray if input is not a masked array
    if WasND: c_ppm, _ = MAtoSeries(c_ppm)
    return c_ppm

def co2_mgCO2pm3fromppm(c_ppm,T,p):
    """
     Convert CO2 concentration units of umol/mol (ppm) to mgCO2/m3
        Usage:
         CO2_mgpm3 = co2_mgCO2pm3fromppm(CO2_ppm, T, p)
         where
         CO2_ppm (input) - CO2 concentration, umol/mol
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the CO2 concentration in mg/m^3.
    """
    # convert to masked array if required
    c_ppm, WasND = SeriestoMA(c_ppm)
    T, dummy = SeriestoMA(T)
    T = T + 273.15             # temperature in K
    p, dummy = SeriestoMA(p)
    p = p * float(1000)        # pressure in Pa
    # do the job
    c_mgpm3 = c_ppm*c.Mco2*p/(c.R*T)
    # convert back to ndarray if input is not a masked array
    if WasND: c_mgpm3, _ = MAtoSeries(c_mgpm3)
    return c_mgpm3

def co2_umolpm3fromppm(c_ppm, T, p):
    """
     Convert CO2 concentration units of umol/mol (ppm) to umol/m^3
        Usage:
         CO2_umolpm3 = co2_umolpm3fromppm(CO2_ppm, T, p)
         where
         CO2_ppm (input) - CO2 concentration, umol/mol
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the CO2 concentration in umol/m^3.
    """
    # convert to masked array if required
    c_ppm, WasND = SeriestoMA(c_ppm)
    T, dummy = SeriestoMA(T)
    T = T + 273.15             # temperature in K
    p, dummy = SeriestoMA(p)
    p = p * float(1000)        # pressure in Pa
    # do the job
    c_umolpm3 = c_ppm*p/(c.R*T)
    # convert back to ndarray if input is not a masked array
    if WasND: c_umolpm3, _ = MAtoSeries(c_umolpm3)
    return c_umolpm3

def densitydryair(Ta,ps,vp):
    # Calculate density of dry air from temperature, pressure and vapour pressure
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  rhod - dry air density, kg/m3
    rhod = 1000*(ps-vp)/((Ta+273.15)*c.Rd)
    return rhod

def densitymoistair(Ta,ps,vp):
    # Calculate density of moist air from temperature, pressure and vapour pressure
    #  Ta - air temperature, C
    #  ps - pressure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  rhom - moist air density, kg/m3
    rhow = vp*1000/((Ta+273.15)*c.Rv)
    rhod = (ps-vp)*1000/((Ta+273.15)*c.Rd)
    rhom = rhow + rhod
    return rhom

def densitywatervapour(Ta,vp):
    # Calculate partial density of water vapour from temperature and vapour pressure
    #  Ta - air temperature, C
    #  vp - vapour pressure, kPa
    # Returns
    #  rhow - partial density of water vapour, kg/m3
    rhow = vp*1000/((Ta+273.15)*c.Rv)
    return rhow

def VPsat(T):
    # Saturation vapour pressure.
    #  T is the air temperature, C
    #  VPsat is the saturation vapour pressure in kPa
    return 0.6106 * numpy.exp(17.27 * T / (T + 237.3))

def Fco2_gCpm2psfromumolpm2ps(Fc_umolpm2ps):
    """
    Convert Fco2 in units of umol/m^2/s to units of g/m^2/s (C, not CO2)
    Usage:
     Fco2_mgpm2ps = Fco2_gCpm2psfromumolpm2ps(Fc_umolpm2ps)
     where:
      Fco2_umolpm2ps (input) - CO2 flux in units of umol/m^2/s
    Returns the CO2 flux in units of g/m^2/s (C, not CO2)
    """
    # convert to masked array
    Fc_umolpm2ps, WasND = SeriestoMA(Fc_umolpm2ps)
    # do the job
    Fc_gCpm2ps = Fc_umolpm2ps * c.Mc/1E3
    # convert back to ndarray if input is not a masked array
    if WasND: Fc_gCpm2ps, _ = MAtoSeries(Fc_gCpm2ps)
    return Fc_gCpm2ps

def Fco2_mgCO2pm2psfromumolpm2ps(Fc_umolpm2ps):
    """
    Convert Fc in units of umol/m^2/s to units of mgCO2/m2/s
    Usage:
     Fc_mgCO2pm2ps = Fco2_mgCO2pm2psfromumolpm2ps(Fc_umolpm2ps)
     where:
      Fc_umolpm2ps (input) - CO2 flux in units of umol/m^2/s
    Returns the CO2 flux in units of mgCO2/m2/s
    """
    # convert to masked array
    Fc_umolpm2ps, WasND = SeriestoMA(Fc_umolpm2ps)
    # do the job
    Fc_mgCO2pm2ps = Fc_umolpm2ps * c.Mco2
    # convert back to ndarray if input is not a masked array
    if WasND: Fc_mgCO2pm2ps, _ = MAtoSeries(Fc_mgCO2pm2ps)
    return Fc_mgCO2pm2ps

def Fco2_umolpm2psfrommgCO2pm2ps(Fc_mgpm2ps):
    """
    Convert Fc in units of mg/m^2/s to units of umol/m^2/s
    Usage:
     Fc_umolpm2ps = Fco2_umolpm2psfrommgCO2pm2ps(Fc_mgpm2ps)
     where:
      Fc_mgpm2ps (input) - CO2 flux in units of mg/m^2/s
    Returns the CO2 flux in units of umol/m^2/s
    """
    # convert to masked array
    Fc_mgpm2ps, WasND = SeriestoMA(Fc_mgpm2ps)
    # do the job
    Fc_umolpm2ps = Fc_mgpm2ps / c.Mco2
    # convert back to ndarray if input is not a masked array
    if WasND: Fc_umolpm2ps, _ = MAtoSeries(Fc_umolpm2ps)
    return Fc_umolpm2ps

def h2o_gpm3frommmolpmol(h_mmpm,T,p):
    """
     Convert H2O concentration units of mmol/mol to g/m^3.
        Usage:
         H2O_gpm3 = h2o_gpm3frommmolpmol(H2O_mmolpmol, T, p)
         where
         H2O_mmolpmol (input) - H2O concentration, mmol/mol
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the H2O concentration in g/m^3.
    """
    # convert to masked arrays
    h_mmpm, WasND = SeriestoMA(h_mmpm)
    T, dummy = SeriestoMA(T)
    p, dummy = SeriestoMA(p)
    # do the job
    h_gpm3 = (c.Mv*h_mmpm*p*1000)/(c.R*(T+273.15))
    # convert to ndarray if input is not a masked array
    if WasND: h_gpm3, _ = MAtoSeries(h_gpm3)
    return h_gpm3

def h2o_mmolpmolfromgpm3(h_gpm3,T,p):
    """
     Convert H2O concentration units of g/m^3 to mmol/mol.
        Usage:
         H2O_mmolpmol = h2o_mmolpmolfromgpm3(H2O_gpm3, T, p)
         where
         H2O_gpm3 (input) - H2O concentration, g/m^3
         T (input) - air temperature, C
         p (input) - air pressure, kPa
        Returns the H2O concentration in mmol/mol.
    """
    # convert to masked arrays
    h_gpm3, WasND = SeriestoMA(h_gpm3)
    T, dummy = SeriestoMA(T)
    p, dummy = SeriestoMA(p)
    # do the job
    h_mmpm = (h_gpm3/c.Mv)*c.R*(T+273.15)/(p*1000)
    # convert to ndarray if input is not a masked array
    if WasND: h_mmpm, _ = MAtoSeries(h_mmpm)
    return h_mmpm

def h2o_mmolpm3fromgpm3(h_gpm3):
    """
     Convert H2O concentration units of g/m^3 to mmol/m^3.
        Usage:
         H2O_mmolpm3 = h2o_mmolpm3fromgpm3(H2O_gpm3)
         where
         H2O_gpm3 (input) - H2O concentration, g/m^3
        Returns the H2O concentration in mmol/m^3.
    """
    # convert to masked arrays
    h_gpm3, WasND = SeriestoMA(h_gpm3)
    # do the job
    h_mmolpm3 = h_gpm3/float(c.Mv)
    # convert to ndarray if input is not a masked array
    if WasND: h_mmolpm3, _ = MAtoSeries(h_mmolpm3)
    return h_mmolpm3

def Lv(Ta):
    # Calculate Lv as a function of temperature, from Stull 1988
    #  Ta - air temperature, C
    # Returns
    #  Lv - latent heat of vapourisation, J/kg
    return 2500800 - (2366.8 * Ta)

def mixingratio(ps,vp):
    # Calculate mixing ratio from vapour pressure and pressure
    #  ps - presure, kPa
    #  vp - vapour pressure, kPa
    # Returns
    #  mr - mixing ratio, kg/kg
    return 0.622*vp/(ps- vp)

def molen(T,AH,p,ustar,heatflux,fluxtype='sensible'):
    # Calculate the Monin-Obukhov length
    ustar = numpy.ma.sqrt(numpy.square(ustar))    # force the sign of ustar to be positive
    vp = vapourpressure(AH,T)       # calculate the vapour pressure
    mr = mixingratio(p,vp)          # calculate the mixing ratio
    Tv = theta(T,p)                 # calculate potential temperature
    Tvp = virtualtheta(Tv,mr)
    if fluxtype=='sensible':
        L = -Tvp*densitydryair(T, p, vp)*c.Cp*(ustar**3)/(c.g*c.k*heatflux)
    elif fluxtype=='kinematic':
        L = -Tvp*(ustar**3)/(c.g*c.k*heatflux)
    else:
        raise Exception(" meteorologicalfunctions.molen: unkown value for fluxtype (="+str(fluxtype)+") encountered")
    return L

def specifichumidityfromrelativehumidity(RH, T, p):
    # Specific humidity (kg/kg) from relative humidity, temperature and pressure
    #  RH is the relative humidity, percent
    #  T is the air temperature, degC
    #  p is the atmospheric pressure, kPa
    SH = (c.Mv / c.Md) * (0.01 * RH * VPsat(T) / p)
    return SH

def SHsat(vpsat, ps):
    return 0.622 * (vpsat / ps)

def relativehumidityfromabsolutehumidity(AH, Ta):
    # Relative humidity from absolute humidity
    #  Ta is the air temperature, degC
    #  AH is the absolute humidity, g/m^3
    #  RH is the relative humidity, percent
    # convert to masked arrays
    AH, WasND = SeriestoMA(AH)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    VP = vapourpressure(AH, Ta)
    RH = float(100)*VP/VPsat(Ta)
    # convert back to ndarray if input is not a masked array
    if WasND: RH, _ = MAtoSeries(RH)
    return RH

def relativehumidityfromdewpoint(Td,Ta):
    # Relative humidity from dew point temperature
    #  Ta is the air temperature, C
    #  Td is the dew point temperature, C
    #  RH is the relative humidity, %
    # convert to masked arrays
    Td, WasND = SeriestoMA(Td)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    RH = 100*10**(7.591386*(Td/(Td+240.7263)-Ta/(Ta+240.7263)))
    # convert back to ndarray if input is not a masked array
    if WasND: RH, _ = MAtoSeries(RH)
    return RH

def relativehumidityfromspecifichumidity(SH,Ta,ps):
    # Relative humidity from specific humidity
    #  SH is the specific humidity, kg/kg
    #  Ta is the air temperature, C
    #  ps is the pressure, kPa
    #  RH is the relative humidity, %
    # convert to masked arrays
    SH, WasND = SeriestoMA(SH)
    Ta, dummy = SeriestoMA(Ta)
    # do the job
    RH = float(100)*SH*(c.Md/c.Mv)*ps/VPsat(Ta)
    # convert back to ndarray if input is not a masked array
    if WasND: RH, _ = MAtoSeries(RH)
    return RH

def densitytimesspecificheat(rhow,Cpw,rhoa,Cpa):
    '''
    Product of air density and specific heat capacity for moist air.
    '''
    return rhow*Cpw+rhoa*Cpa

def specificheatcapacitydryair(Tv):
    '''
    Specific heat capacity of air at constant pressure.
    USEAGE:
     cpd = pfp_mf.specificheatcapacitydryair(Tv)
    INPUT:
     Tv - virtual temperature (from sonic anemometer), C
    OUTPUT:
     cpd - specific heat capacity of dry air at constant pressure, J/kg/K
    SOURCE:
     EddyPro source code
    '''
    cpd = float(1005)+((Tv+23.12)**2)/float(3364)
    return cpd

def specificheatcapacitywatervapour(Ta, AH):
    '''
    Specific heat capacity of water vapour at constant pressure.
    USEAGE:
     cpv = pfp_mf.specificheatcapacitywatervapour(Ta,AH)
    INPUT:
     Ta - air temperature, C
     AH - absolute humidity, %
    OUTPUT:
     cpv - specific heat capacity of water vapour at constant pressure, J/kg/K
    SOURCE:
     EddyPro source code
    '''
    RH = relativehumidityfromabsolutehumidity(AH,Ta)
    cpv = float(1859)+0.13*RH+(0.193+0.0056*RH)*Ta+(0.001+0.00005*RH)*Ta**2
    return cpv

def specificheatmoistair(SH):
    # Calculate Cp of moist air, from Stull 1988
    #  Cp - specific heat of dry air at constant pressure, J/kg-K
    #  SH - specific humidity
    # Returns
    #  Cpm - specific heat of moist air at constant pressure, J/kg-K
    Cpm = c.Cpd * (1 + 0.84 * SH)
    return Cpm

def specifichumidity(mr):
    # Calculate specific humidity from mixing ratio
    #  mr - mixing ration, kg/kg
    # Returns
    #  SH = specific humidity, kg/kg
    SH = mr/(1+mr)
    return SH

def specifichumidityfromRH(RH, T, p):
    # Specific humidity (kg/kg) from relative humidity, temperature and pressure
    #  RH is the relative humidity, %
    #  T is the air temperature, C
    #  p is the atmospheric pressure, kPa
    # Returns
    #  SH = specific humidity, kg/kg
    SH = (c.Mv / c.Md) * (0.01 * RH * VPsat(T) / p)
    return SH

def tafromtv(Tv,SH):
    # Calculate air temperature from virtual temperature using formula
    # from Campbell Scientific CSAT manual.
    # NOTE: this differs from the usual definition by using 0.51 not 0.61
    #  Tv - virtual temperature, C
    #  SH - specific humidity, kg/kg
    # Returns
    #  Ta - air temperature, C
    Ta = ((Tv+273.15)/(1+0.51*SH))-273.15
    return Ta

def tvfromta(Ta, mr):
    # Calculate virtual temperature from air temperature using formula
    # from Stull 1988.
    #  Ta - virtual temperature, C
    #  mr - H2O mixing ratio
    # Returns
    #  Tv - virtual temperature, C
    Tv = ((Ta + 273.15)*(1 + 0.61*mr)) - 273.15
    return Tv

def theta(T,p):
    # Calculate potential temperature from air temperature and pressure
    #  T - air temperature, C
    #  p - pressure, kPa
    # Returns
    #  theta - potential temperature, K
    return (T+273.15)*(100/p)**0.286

def vapourpressure(AH,Ta):
    # Calculate vapour pressure from absolute humidity and temperature
    #  AH - absolute humidity, g/m^3
    #  Ta - air temperature, C
    # Returns
    #  vp - vapour pressure, kPa
    vp = 0.000001*AH*(Ta+273.15)*c.R/c.Mv
    return vp

def virtualtheta(theta_in, mr):
    # Calculate virtual potential temperature
    #  theta - potential temperature, K
    #  mr - mixing ratio, kg/kg
    # Returns
    #  Tvp - virtual potential temperature, K
    Tvp = theta_in * (1 + (0.61 * mr))
    return Tvp
