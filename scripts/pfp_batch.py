# standard modules
import datetime
import logging
import ntpath
import os
import sys
import traceback
# 3rd party modules
from configobj import ConfigObj
# PFP modules
#sys.path.append("scripts")
from scripts import pfp_cfg
from scripts import pfp_clim
from scripts import pfp_compliance
from scripts import pfp_cpd1
from scripts import pfp_cpd2
from scripts import pfp_io
from scripts import pfp_levels
from scripts import pfp_mpt
from scripts import pfp_plot
from scripts import pfp_utils

logger = logging.getLogger("pfp_log")

def do_L1_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L1")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        logger.info("Starting L1 processing with %s", cf_file_name[1])
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            ds1 = pfp_levels.l1qc(cf)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds1)
            msg = "Finished L1 processing with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during L1 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_L2_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L2")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting L2 processing with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            infilename = pfp_io.get_infilenamefromcf(cf)
            ds1 = pfp_io.nc_read_series(infilename)
            if ds1.returncodes["value"] != 0: return
            ds2 = pfp_levels.l2qc(cf, ds1)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds2)
            msg = "Finished L2 processing with " + cf_file_name[1]
            logger.info(msg)
            if "Plots" in list(cf.keys()):
                logger.info("Plotting L1 and L2 data")
                for nFig in list(cf['Plots'].keys()):
                    plt_cf = cf['Plots'][str(nFig)]
                    if 'type' in plt_cf.keys():
                        if str(plt_cf['type']).lower() == 'xy':
                            pfp_plot.plotxy(cf, nFig, plt_cf, ds1, ds2)
                        else:
                            pfp_plot.plottimeseries(cf, nFig, ds1, ds2)
                    else:
                        pfp_plot.plottimeseries(cf, nFig, ds1, ds2)
                logger.info("Finished plotting L1 and L2 data")
            logger.info("")
        except Exception:
            msg = "Error occurred during L2 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_L3_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L3")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting L3 processing with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            infilename = pfp_io.get_infilenamefromcf(cf)
            ds2 = pfp_io.nc_read_series(infilename)
            if ds2.returncodes["value"] != 0: return
            ds3 = pfp_levels.l3qc(cf, ds2)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds3)
            msg = "Finished L3 processing with " + cf_file_name[1]
            logger.info(msg)
            if "Plots" in list(cf.keys()):
                logger.info("Plotting L3 data")
                for nFig in list(cf['Plots'].keys()):
                    plt_cf = cf['Plots'][str(nFig)]
                    if 'type' in plt_cf.keys():
                        if str(plt_cf['type']).lower() == 'xy':
                            pfp_plot.plotxy(cf, nFig, plt_cf, ds2, ds3)
                        else:
                            pfp_plot.plottimeseries(cf, nFig, ds2, ds3)
                    else:
                        pfp_plot.plottimeseries(cf, nFig, ds2, ds3)
                logger.info("Finished plotting L3 data")
            logger.info("")
        except Exception:
            msg = "Error occurred during L3 processing " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_ecostress_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "ecostress")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting ECOSTRESS output with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            pfp_io.write_csv_ecostress(cf)
            msg = "Finished ECOSTRESS output with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during ECOSTRESS output with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_fluxnet_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "fluxnet")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting FluxNet output with " + cf_file_name[1]
        logger.info(msg)
        cf = pfp_io.get_controlfilecontents(cf_level[i])
        pfp_io.fn_write_csv(cf)
        msg = "Finished FluxNet output with " + cf_file_name[1]
        logger.info(msg)
        logger.info("")
    return
def do_reddyproc_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "reddyproc")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting REddyProc output with " + cf_file_name[1]
        logger.info(msg)
        cf = pfp_io.get_controlfilecontents(cf_level[i])
        pfp_io.reddyproc_write_csv(cf)
        msg = "Finished REddyProc output with " + cf_file_name[1]
        logger.info(msg)
        logger.info("")
    return
def do_concatenate_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "concatenate")
    for i in sorted(list(cf_level.keys())):
        if not os.path.isfile(cf_level[i]):
            msg = " Control file " + cf_level[i] + " not found"
            logger.error(msg)
            continue
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting concatenation with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf_cc = pfp_io.get_controlfilecontents(cf_level[i])
            info = pfp_compliance.ParseConcatenateControlFile(cf_cc)
            if not info["NetCDFConcatenate"]["OK"]:
                msg = " Error occurred parsing the control file " + cf_file_name[1]
                logger.error(msg)
                continue
            pfp_io.NetCDFConcatenate(info)
            msg = "Finished concatenation with " + cf_file_name[1]
            logger.info(msg)
            # now plot the fingerprints for the concatenated files
            opt = pfp_utils.get_keyvaluefromcf(cf_cc, ["Options"], "DoFingerprints", default="no")
            if opt.lower() == "no":
                continue
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp):
                cf_fp["Files"] = {}
            file_name = cf_cc["Files"]["Out"]["ncFileName"]
            file_path = ntpath.split(file_name)[0] + "/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            if "plot_path" in cf_cc["Files"]:
                cf_fp["Files"]["plot_path"] = cf_cc["Files"]["plot_path"]
            else:
                cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")] + "Plots/"
            if "Options" not in cf_fp:
                cf_fp["Options"] = {}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "No"
            msg = "Doing fingerprint plots using " + cf_fp["Files"]["in_filename"]
            logger.info(msg)
            pfp_plot.plot_fingerprint(cf_fp)
            msg = "Finished fingerprint plots"
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during concatenation with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_climatology_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "climatology")
    for i in list(cf_level.keys()):
        if not os.path.isfile(cf_level[i]):
            msg = " Control file " + cf_level[i] + " not found"
            logger.error(msg)
            continue
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting climatology with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            pfp_clim.climatology(cf)
            msg = "Finished climatology with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during climatology with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_cpd1_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "cpd1")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting CPD (McHugh) with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            pfp_cpd1.cpd1_main(cf)
            msg = "Finished CPD (McHugh) with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during CPD (McHugh) with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_cpd2_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "cpd2")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting CPD (Barr) with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            pfp_cpd2.cpd2_main(cf)
            msg = "Finished CPD (Barr) with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during CPD (Barr) with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_mpt_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "mpt")
    for i in list(cf_level.keys()):
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting MPT with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            pfp_mpt.mpt_main(cf)
            msg = "Finished MPT with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during MPT with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_L4_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L4")
    for i in list(cf_level.keys()):
        if not os.path.isfile(cf_level[i]):
            msg = " Control file " + cf_level[i] + " not found"
            logger.error(msg)
            continue
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting L4 processing with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf_l4 = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf_l4:
                cf_l4["Options"] = {}
            cf_l4["Options"]["call_mode"] = "batch"
            cf_l4["Options"]["show_plots"] = "No"
            infilename = pfp_io.get_infilenamefromcf(cf_l4)
            ds3 = pfp_io.nc_read_series(infilename)
            if ds3.returncodes["value"] != 0: return
            ds4 = pfp_levels.l4qc(None, cf_l4, ds3)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l4)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds4)
            msg = "Finished L4 processing with " + cf_file_name[1]
            logger.info(msg)
            # now plot the fingerprints for the L4 files
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp):
                cf_fp["Files"] = {}
            file_name = pfp_io.get_outfilenamefromcf(cf_l4)
            file_path = ntpath.split(file_name)[0] + "/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            if "plot_path" in cf_l4["Files"]:
                cf_fp["Files"]["plot_path"] = cf_l4["Files"]["plot_path"]
            else:
                cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")] + "Plots/"
            if "Options" not in cf_fp:
                cf_fp["Options"] = {}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "No"
            msg = "Doing fingerprint plots using " + cf_fp["Files"]["in_filename"]
            logger.info(msg)
            pfp_plot.plot_fingerprint(cf_fp)
            logger.info("Finished fingerprint plots")
            logger.info("")
        except Exception:
            msg = "Error occurred during L4 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_L5_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L5")
    for i in list(cf_level.keys()):
        if not os.path.isfile(cf_level[i]):
            msg = " Control file " + cf_level[i] + " not found"
            logger.error(msg)
            continue
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting L5 processing with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf_l5 = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf_l5:
                cf_l5["Options"] = {}
            cf_l5["Options"]["call_mode"] = "batch"
            cf_l5["Options"]["show_plots"] = "No"
            infilename = pfp_io.get_infilenamefromcf(cf_l5)
            ds4 = pfp_io.nc_read_series(infilename)
            if ds4.returncodes["value"] != 0: return
            ds5 = pfp_levels.l5qc(None, cf_l5, ds4)
            outfilename = pfp_io.get_outfilenamefromcf(cf_l5)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds5)
            msg = "Finished L5 processing with " + cf_file_name[1]
            logger.info(msg)
            # now plot the fingerprints for the L5 files
            cf_fp = pfp_io.get_controlfilecontents("controlfiles/standard/fingerprint.txt")
            if "Files" not in dir(cf_fp):
                cf_fp["Files"] = {}
            file_name = pfp_io.get_outfilenamefromcf(cf_l5)
            file_path = ntpath.split(file_name)[0] + "/"
            cf_fp["Files"]["file_path"] = file_path
            cf_fp["Files"]["in_filename"] = ntpath.split(file_name)[1]
            if "plot_path" in cf_l5["Files"]:
                cf_fp["Files"]["plot_path"] = cf_l5["Files"]["plot_path"]
            else:
                cf_fp["Files"]["plot_path"] = file_path[:file_path.index("Data")] + "Plots/"
            if "Options" not in cf_fp:
                cf_fp["Options"] = {}
            cf_fp["Options"]["call_mode"] = "batch"
            cf_fp["Options"]["show_plots"] = "No"
            msg = "Doing fingerprint plots using " + cf_fp["Files"]["in_filename"]
            logger.info(msg)
            pfp_plot.plot_fingerprint(cf_fp)
            msg = "Finished fingerprint plots"
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during L5 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_L6_batch(cf_level):
    #logger = pfp_log.change_logger_filename("pfp_log", "L6")
    for i in list(cf_level.keys()):
        if not os.path.isfile(cf_level[i]):
            msg = " Control file " + cf_level[i] + " not found"
            logger.error(msg)
            continue
        cf_file_name = os.path.split(cf_level[i])
        msg = "Starting L6 processing with " + cf_file_name[1]
        logger.info(msg)
        try:
            cf = pfp_io.get_controlfilecontents(cf_level[i])
            if "Options" not in cf:
                cf["Options"] = {}
            cf["Options"]["call_mode"] = "batch"
            cf["Options"]["show_plots"] = "No"
            infilename = pfp_io.get_infilenamefromcf(cf)
            ds5 = pfp_io.nc_read_series(infilename)
            if ds5.returncodes["value"] != 0: return
            ds6 = pfp_levels.l6qc(None, cf, ds5)
            outfilename = pfp_io.get_outfilenamefromcf(cf)
            nc_file = pfp_io.nc_open_write(outfilename)
            if nc_file is None: return
            pfp_io.nc_write_series(nc_file, ds6)
            msg = "Finished L6 processing with " + cf_file_name[1]
            logger.info(msg)
            logger.info("")
        except Exception:
            msg = "Error occurred during L6 with " + cf_file_name[1]
            logger.error(msg)
            error_message = traceback.format_exc()
            logger.error(error_message)
            continue
    return
def do_levels_batch(cf_batch):
    logger = logging.getLogger("pfp_log")
    #batch_log_path = pfp_log.get_batch_log_path("logfiles")
    #batch_log_file_name = os.path.join(batch_log_path, "batch.log")
    #logger = pfp_log.init_logger("pfp_log", batch_log_file_name, to_file=True, to_screen=True)
    start = datetime.datetime.now()
    msg = " Started batch processing at " + start.strftime("%Y%m%d%H%M")
    logger.info(msg)
    if "Options" in cf_batch:
        if "levels" in cf_batch["Options"]:
            levels = pfp_cfg.cfg_string_to_list(cf_batch["Options"]["levels"])
        else:
            msg = "No 'levels' entry found in [Options] section"
            logger.error(msg)
            sys.exit()
    else:
        msg = "No [Options] section in control file"
        logger.error(msg)
        sys.exit()
    processing_levels = ["l1", "l2", "l3",
                         "ecostress", "fluxnet", "reddyproc",
                         "concatenate", "climatology",
                         "cpd_mchugh", "cpd_barr", "mpt",
                         "l4", "l5", "l6"]
    for level in levels:
        if level.lower() not in processing_levels:
            msg = "Unrecognised level " + level
            logger.warning(msg)
            continue
        if level.lower() == "l1":
            # L1 processing
            do_L1_batch(cf_batch["Levels"][level])
        elif level.lower() == "l2":
            # L2 processing
            do_L2_batch(cf_batch["Levels"][level])
        elif level.lower() == "l3":
            # L3 processing
            do_L3_batch(cf_batch["Levels"][level])
        elif level.lower() == "ecostress":
            # convert netCDF files to ECOSTRESS CSV files
            do_ecostress_batch(cf_batch["Levels"][level])
        elif level.lower() == "fluxnet":
            # convert netCDF files to FluxNet CSV files
            do_fluxnet_batch(cf_batch["Levels"][level])
        elif level.lower() == "reddyproc":
            # convert netCDF files to REddyProc CSV files
            do_reddyproc_batch(cf_batch["Levels"][level])
        elif level.lower() == "concatenate":
            # concatenate netCDF files
            do_concatenate_batch(cf_batch["Levels"][level])
        elif level.lower() == "climatology":
            # climatology
            do_climatology_batch(cf_batch["Levels"][level])
        elif level.lower() == "cpd_mchugh":
            # ustar threshold from change point detection
            do_cpd1_batch(cf_batch["Levels"][level])
        elif level.lower() == "cpd_barr":
            # ustar threshold from change point detection
            do_cpd2_batch(cf_batch["Levels"][level])
        elif level.lower() == "mpt":
            # ustar threshold from change point detection
            do_mpt_batch(cf_batch["Levels"][level])
        elif level.lower() == "l4":
            # L4 processing
            do_L4_batch(cf_batch["Levels"][level])
        elif level.lower() == "l5":
            # L5 processing
            do_L5_batch(cf_batch["Levels"][level])
        elif level.lower() == "l6":
            # L6 processing
            do_L6_batch(cf_batch["Levels"][level])
    #logger = pfp_log.change_logger_filename("pfp_log", "batch")
    end = datetime.datetime.now()
    msg = " Finished batch processing at " + end.strftime("%Y%m%d%H%M")
    logger.info(msg)
    return

#if (__name__ == '__main__'):
    ## get the control file name
    #if len(sys.argv) == 1:
        ## not on the command line, so ask the user
        #cfg_file_path = raw_input("Enter the control file name: ")
        ## exit if nothing selected
        #if len(cfg_file_path) == 0:
            #sys.exit()
    #else:
        ## control file name on the command line
        #if not os.path.exists(sys.argv[1]):
            ## control file doesn't exist
            #logger.error("Control file %s does not exist", sys.argv[1])
            #sys.exit()
        #else:
            #cfg_file_path = sys.argv[1]
    ## read the control file
    #cf_batch = ConfigObj(cfg_file_path, indent_type="    ", list_values=False)
    ## call the processing
    #do_levels_batch(cf_batch)
