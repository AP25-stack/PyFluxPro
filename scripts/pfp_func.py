# standard modules
import logging
# 3rd party
import numpy
# PFP modules
import meteorologicalfunctions as pfp_mf
import pfp_utils

logger = logging.getLogger("pfp_log")

def Convert_fraction_to_percent(ds, RH_out, RH_in):
    """
    Purpose:
     Function to convert RH in units of "frac" (0 to 1) to "percent" (1 to 100).
    Usage:
     pfp_func.Convert_fraction_to_percent(ds, RH_out, RH_in)
    Author: PRI
    Date: August 2019
    """
    var_in = pfp_utils.GetVariable(ds, RH_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "%", mode="quiet")
    var_out["Label"] = RH_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_gH2Opm3_to_percent(ds, RH_out, AH_in, Ta_in):
    """
    Purpose:
     Function to convert absolute humidity in units of g/m^3 to relative humidity in percent.
    Usage:
     pfp_func.Convert_gH2Opm3_to_percent(ds, RH_out, AH_in, Ta_in)
    Author: PRI
    Date: September 2020
    """
    nRecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nRecs,dtype=numpy.int32)
    ones = numpy.ones(nRecs,dtype=numpy.int32)
    for item in [AH_in, Ta_in]:
        if item not in ds.series.keys():
            msg = " Requested series " + item + " not found, " + RH_out + " not calculated"
            logger.error(msg)
            return 0
    AH = pfp_utils.GetVariable(ds, AH_in)
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    RH = pfp_utils.GetVariable(ds, RH_out)
    RH["Data"] = pfp_mf.relativehumidityfromabsolutehumidity(AH["Data"], Ta["Data"])
    RH["Flag"] = numpy.where(numpy.ma.getmaskarray(RH["Data"]) == True, ones, zeros)
    pfp_utils.CreateVariable(ds, RH)
    return 1

def Convert_gH2Opm3_to_mmolpm3(ds, H2O_out, AH_in):
    """
    Purpose:
     Calculate H2O molar density in mmol/m^3 from absolute humidity in g/m^3.
    Usage:
     pfp_func.Convert_gH2Opm3_to_mmolpm3(ds, MD_out, AH_in)
    Author: PRI
    Date: September 2020
    """
    #nRecs = int(ds.globalattributes["nc_nrecs"])
    #zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    #ones = numpy.ones(nRecs, dtype=numpy.int32)
    #for item in [Ah_in]:
        #if item not in ds.series.keys():
            #msg = " Requested series " + item + " not found, " + MD_out + " not calculated"
            #logger.error(msg)
            #return 0
    #Ah = pfp_utils.GetVariable(ds, Ah_in)
    #MD = pfp_utils.GetVariable(ds, MD_out)
    #MD["Data"] = pfp_mf.h2o_mmolpm3fromgpm3(Ah["Data"])
    #MD["Flag"] = numpy.where(numpy.ma.getmaskarray(MD["Data"]) == True, ones, zeros)
    #pfp_utils.CreateVariable(ds, MD)
    #return 1
    for item in [AH_in]:
        if item not in ds.series.keys():
            msg = " Requested series " + item + " not found, " + H2O_out + " not calculated"
        logger.error(msg)
        return 0
    var_in = pfp_utils.GetVariable(ds, AH_in)
    got_variance = False
    if "Vr" in var_in["Label"] and ")2" in var_in["Attr"]["units"]:
        got_variance = True
        var_in["Data"] = numpy.ma.sqrt(var_in["Data"])
        var_in["Attr"]["units"] = pfp_utils.units_variance_to_standard_deviation(var_in["Attr"]["units"])
    var_out = pfp_utils.convert_units_func(ds, var_in, "mmol/m^3", mode="quiet")
    var_out["Label"] = H2O_out
    if got_variance:
        var_out["Data"] = var_out["Data"]*var_out["Data"]
        var_out["Attr"]["units"] = pfp_utils.units_standard_deviation_to_variance(var_out["Attr"]["units"])
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_gH2Opm3_to_mmolpmol(ds, MF_out, AH_in, Ta_in, ps_in):
    """
    Purpose:
     Calculate H2O mole fraction in mml/mol from absolute humidity in g/m^3.
    Usage:
     pfp_func.Convert_gH2Opm3_to_mmolpmol(ds, MF_out, AH_in, Ta_in, ps_in)
    Author: PRI
    Date: August 2019
    """
    nRecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nRecs,dtype=numpy.int32)
    ones = numpy.ones(nRecs,dtype=numpy.int32)
    for item in [AH_in, Ta_in, ps_in]:
        if item not in ds.series.keys():
            msg = " Requested series " + item + " not found, " + MF_out + " not calculated"
            logger.error(msg)
            return 0
    AH = pfp_utils.GetVariable(ds, AH_in)
    Ta = pfp_utils.GetVariable(ds, Ta_in)
    ps = pfp_utils.GetVariable(ds, ps_in)
    MF = pfp_utils.GetVariable(ds, MF_out)
    MF["Data"] = pfp_mf.h2o_mmolpmolfromgpm3(Ah["Data"], Ta["Data"], ps["Data"])
    MF["Flag"] = numpy.where(numpy.ma.getmaskarray(MF["Data"]) == True, ones, zeros)
    pfp_utils.CreateVariable(ds, MF)
    return 1

def Convert_hPa_to_kPa(ds, ps_out, ps_in):
    """
    Purpose:
     Function to convert pressure from hPa (mb) to kPa.
    Usage:
     pfp_func.ConverthPa2kPa(ds, ps_in, ps_out)
    Author: PRI
    Date: February 2018
    """
    var_in = pfp_utils.GetVariable(ds, ps_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "kPa", mode="quiet")
    var_out["Label"] = ps_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_K_to_C(ds, T_out, T_in):
    """
    Purpose:
     Function to convert temperature from K to C.
    Usage:
     pfp_func.Convert_K_to_C(ds, T_out, T_in)
    Author: PRI
    Date: February 2018
    """
    if T_in not in list(ds.series.keys()):
        msg = " ConvertK2C: variable " + T_in + " not found, skipping ..."
        logger.warning(msg)
        return 0
    if "<" in T_out or ">" in T_out:
        logger.warning(" ***")
        msg = " *** " + T_in + ": illegal name (" + T_out + ") in function, skipping ..."
        logger.warning(msg)
        logger.warning(" ***")
        return 0
    var_in = pfp_utils.GetVariable(ds, T_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "degC", mode="quiet")
    var_out["Label"] = T_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_kgpm3_to_gpm3(ds, AH_out, AH_in):
    """
    Purpose:
     Function to convert absolute humidity from kg/m^3 to g/m^3.
    Usage:
     pfp_func.Convertkgpm32gpm3(ds, Ah_out, Ah_in)
    Author: PRI
    Date: August 2020
    """
    var_in = pfp_utils.GetVariable(ds, AH_in)
    var_out = pfp_utils.convert_units_func(ds, var_in, "g/m^3", mode="quiet")
    var_out["Label"] = AH_out
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_mmolpm3_to_gH2Opm3(ds, AH_out, H2O_in):
    """
    Purpose:
     Function to convert mmol/m^3 (molar density) to g/m^3 (mass density).
    Usage:
     pfp_func.Convert_mmolpm3_to_gpm3(ds, AH_out, H2O_in)
    Author: PRI
    Date: August 2020
    """
    for item in [H2O_in]:
        if item not in list(ds.series.keys()):
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
            return 0
    var_in = pfp_utils.GetVariable(ds, H2O_in)
    got_variance = False
    if "Vr" in var_in["Label"] and ")2" in var_in["Attr"]["units"]:
        got_variance = True
        var_in["Data"] = numpy.ma.sqrt(var_in["Data"])
        var_in["Attr"]["units"] = pfp_utils.units_variance_to_standard_deviation(var_in["Attr"]["units"])
    var_out = pfp_utils.convert_units_func(ds, var_in, "g/m^3", mode="quiet")
    var_out["Label"] = AH_out
    if got_variance:
        var_out["Data"] = var_out["Data"]*var_out["Data"]
        var_out["Attr"]["units"] = pfp_utils.units_standard_deviation_to_variance(var_out["Attr"]["units"])
    pfp_utils.CreateVariable(ds, var_out)
    return 1

def Convert_mmolpmol_to_gH2Opm3(ds, AH_out, MF_in, Ta_in, ps_in):
    """
    Purpose:
     Function to calculate absolute humidity given the water vapour mole
     fraction, air temperature and pressure.  Absolute humidity is not calculated
     if any of the input series are missing or if the specified output series
     already exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     pfp_func.Convert_mmolpmol_to_gpm3(ds,"AH_IRGA_Av","H2O_IRGA_Av","Ta_HMP_2m","ps")
    Author: PRI
    Date: September 2015
    """
    nRecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nRecs, dtype=numpy.int32)
    ones = numpy.ones(nRecs, dtype=numpy.int32)
    for item in [MF_in, Ta_in, ps_in]:
        if item not in list(ds.series.keys()):
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
        return 0
    if AH_out in ds.series.keys():
        msg = " Output series " + AH_out + " already exists, skipping ..."
        logger.error(msg)
        return 0
    MF_data,MF_flag,MF_attr = pfp_utils.GetSeriesasMA(ds,MF_in)
    Ta_data,Ta_flag,Ta_attr = pfp_utils.GetSeriesasMA(ds,Ta_in)
    ps_data,ps_flag,ps_attr = pfp_utils.GetSeriesasMA(ds,ps_in)
    AH_data = pfp_mf.h2o_gpm3frommmolpmol(MF_data,Ta_data,ps_data)
    long_name = "Absolute humidity calculated from " + MF_in + ", " + Ta_in + " and " + ps_in
    AH_attr = pfp_utils.MakeAttributeDictionary(long_name=long_name,
                                              height=MF_attr["height"],
                                              units="g/m^3")
    AH_flag = numpy.where(numpy.ma.getmaskarray(AH_data) == True, ones, zeros)
    pfp_utils.CreateSeries(ds, AH_out, AH_data, AH_flag, AH_attr)
    return 1

def Convert_percent_to_mmolpmol(ds, MF_out, RH_in, Ta_in, ps_in):
    """
    Purpose:
     Calculate H2O mole fraction from relative humidity (RH).
    """
    nRecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nRecs,dtype=numpy.int32)
    ones = numpy.ones(nRecs,dtype=numpy.int32)
    for item in [RH_in, Ta_in, ps_in]:
        if item not in list(ds.series.keys()):
            msg = " Requested series " + item + " not found, " + MF_out + " not calculated"
            logger.error(msg)
            return 0
    if MF_out in list(ds.series.keys()):
        msg = " Output series " + MF_out + " already exists, skipping ..."
        logger.error(msg)
        return 0
    RH_data, RH_flag, RH_attr = pfp_utils.GetSeriesasMA(ds, RH_in)
    Ta_data,Ta_flag,Ta_attr = pfp_utils.GetSeriesasMA(ds, Ta_in)
    AH_data = pfp_mf.absolutehumidityfromrelativehumidity(Ta_data, RH_data)
    ps_data,ps_flag,ps_attr = pfp_utils.GetSeriesasMA(ds, ps_in)
    MF_data = pfp_mf.h2o_mmolpmolfromgpm3(AH_data, Ta_data, ps_data)
    long_name = "H2O mole fraction calculated from " + RH_in + ", " + Ta_in + " and " + ps_in
    MF_attr = pfp_utils.MakeAttributeDictionary(long_name=long_name,
                                              height=RH_attr["height"],
                                              units="mmol/mol")
    MF_flag = numpy.where(numpy.ma.getmaskarray(MF_data)==True,ones,zeros)
    pfp_utils.CreateSeries(ds, MF_out, MF_data, MF_flag, MF_attr)
    return 1

def Convert_percent_to_gH2Opm3(ds, AH_out, RH_in, Ta_in):
    """
    Purpose:
     Function to calculate absolute humidity given relative humidity and
     air temperature.  Absolute humidity is not calculated if any of the
     input series are missing or if the specified output series already
     exists in the data structure.
     The calculated absolute humidity is created as a new series in the
     data structure.
    Usage:
     pfp_func.Convert_percent_to_gpm3(ds,"AH_HMP_2m","RH_HMP_2m","Ta_HMP_2m")
    Author: PRI
    Date: September 2015
    """
    nRecs = int(ds.globalattributes["nc_nrecs"])
    zeros = numpy.zeros(nRecs,dtype=numpy.int32)
    ones = numpy.ones(nRecs,dtype=numpy.int32)
    for item in [RH_in, Ta_in]:
        if item not in ds.series.keys():
            msg = " Requested series " + item + " not found, " + AH_out + " not calculated"
            logger.error(msg)
            return 0
    if AH_out in list(ds.series.keys()):
        msg = " Output series " + AH_out + " already exists, skipping ..."
        logger.error(msg)
        return 0
    RH_data,RH_flag,RH_attr = pfp_utils.GetSeriesasMA(ds, RH_in)
    Ta_data,Ta_flag,Ta_attr = pfp_utils.GetSeriesasMA(ds, Ta_in)
    AH_data = pfp_mf.absolutehumidityfromrelativehumidity(Ta_data, RH_data)
    long_name = "Absolute humidity calculated from " + RH_in + " and " + Ta_in
    AH_attr = pfp_utils.MakeAttributeDictionary(long_name=long_name,
                                              height=RH_attr["height"],
                                                units="g/m^3")
    flag = numpy.where(numpy.ma.getmaskarray(AH_data) == True, ones, zeros)
    pfp_utils.CreateSeries(ds, AH_out, AH_data, flag, AH_attr)
    return 1
