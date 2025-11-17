###############################################################################################################
###  CREATING OF 12-DAY MOSAICS of S1 DATA FOR S14AMAZONAS
### Script written by Dr. Neha Hunka, 10.02.2021
# The available processed backscatter files are grouped into 12-day intervals based on the observation date. The files
# which belong a particular group are merged together, according to their polarization.
### Script execute with command -
# python3 SARbackscatter_mosaic.py --in_tile "21LYG" --backscatter_folder ".../output/GRD"
# --mosaic_folder ".../output/Mosaic" --aux_folder ".../aux_data"
# --in_multipolygon_string "MultiPolygon(((-55.1650695389999 -11.753601353,-54.158203488 -11.7453511209999,-54.147570924 -12.737064393,-55.15819755 -12.746032956,-55.1650695389999 -11.753601353)))"
# --in_timeperiod "20150601-20171231"
###############################################################################################################






# import packages
import re
import logging
import os
import requests
from pathlib import Path
from tempfile import mkdtemp
import time
import subprocess
import argparse

from csv import reader
import json
from datetime import datetime, timedelta
from datetime import date, timedelta
from time import gmtime

import shutil

try:
    import ogr
except:
    from osgeo import ogr

try:
    import gdal
except:
    from osgeo import gdal

import numpy as np
import glob


import json
from pathlib import Path
from datetime import datetime








HLINE = '---------------------------------------------------------------------------------------------------------'
log = logging.getLogger(__name__)

### GDAL and SNAP settings ###
memsize = 2048
co = ['COMPRESS=DEFLATE', 'TILED=YES']

### CREATE A GRID-LIST OF DATES EVERY 6/12 DAYS (for 12-day acquisition frequency) ###
acq_frequency = 12

### DECLARE IF MOSAICS OF BACKSCATTER, COHERENCE, OR BOTH, ARE TO BE PREPARED, and IF INDIVIDUAL ORBITS ARE NEEDED ###
TYPE = ['_BAC_']                      # backscatter mosaics only
Individual_orbits_COH_ = 'FALSE'      # not output single orbit coherence images
Individual_orbits_BAC_ = 'FALSE'      # not output single-orbit backscatter images


# making a logger
log = logging.getLogger(__name__)     # create logger object
def init_logging():
    logging.Formatter.converter = gmtime             # use UTC (GMT)
    console_handler = logging.StreamHandler()        # outputs outputs log messages
    console_handler.setFormatter(logging.Formatter(fmt="{levelname} {message}", style="{"))     # set log message format
    root_log = logging.getLogger()                   # top-level logger for the app
    root_log.addHandler(console_handler)             # console handler to the root
    root_log.setLevel(logging.NOTSET)                # 
    log.info("Logging has been started.")            # logging initialized


# 
def get_orbit_list_from_files(folder):
    orbit_list = set()

    folder = Path(folder)

    for f in folder.iterdir():
        name = f.name

        # Only process Sentinel-1 GRD SAFE folders or ZIPs
        if not (name.startswith("S1A") or name.startswith("S1B")):
            continue

        parts = name.replace(".SAFE", "").split("_")

        if len(parts) < 7:
            continue  # unexpected format

        satellite_id = parts[0]      # S1A or S1B
        abs_orbit = int(parts[6])    # absolute orbit number

        # Convert to relative orbit
        if satellite_id == "S1A":
            rel_orbit = ((abs_orbit - 73) % 175) + 1
        elif satellite_id == "S1B":
            rel_orbit = ((abs_orbit - 27) % 175) + 1
        else:
            continue

        orbit_list.add(rel_orbit)

    return sorted(orbit_list)










#################################### THE START ##################################################
def run(tile, backscatter_folder, mosaic_folder, aux_folder, multipolygon_string, timeperiod):

    aux_folder = Path(aux_folder)

    # find the limits of the extent
    aoi_geom = ogr.CreateGeometryFromWkt(multipolygon_string)
    (xmin, xmax, ymin, ymax) = aoi_geom.GetEnvelope()
    ulx = xmin
    uly = ymax
    lrx = xmax
    lry = ymin

    init_logging()

    # folder where the backscatters of all tiles are stored
    BAC_folder = Path(backscatter_folder)
    # folder where the backscatter mosaics will be stored
    work_dir = Path(mosaic_folder)

    # replace spaces in the multipolygon string with %20 for sending url request
    geometry_wkt = multipolygon_string.replace(' ','%20')

    # list of dates with 12 days interval
    startdatestr = timeperiod.split('-')[0]
    enddatestr = timeperiod.split('-')[1]
    start_date =  datetime.strptime(startdatestr, '%Y%m%d')
    end_date = datetime.strptime(enddatestr, '%Y%m%d')
    days_interval = np.ndarray.tolist(
        np.arange(start_date, end_date, timedelta(days=acq_frequency), dtype='datetime64[D]'))
    log.debug("performing mosaic for date above {}".format(start_date))

    # use ascending and descending orbits for any region
    ORBIT_DIRECTIONS = ['descending', 'ascending']



    ######################################################################################################
    MOS_folder = Path(work_dir)
    log.debug("Mosaics folder is {}".format(MOS_folder))
    ######################################################################################################
    for ORBIT_DIRECTION in ORBIT_DIRECTIONS:
        log.debug("CURRENTLY PROCESSING ... {} for tile ... {}".format(ORBIT_DIRECTION, tile))
        # accept either <backscatter_folder>/<tile> or directly <backscatter_folder>
        src_candidates = [Path(BAC_folder).joinpath(tile), Path(BAC_folder)]
        SRC_folder_B = next((p for p in src_candidates if p.exists()), Path(BAC_folder).joinpath(tile))
        #SRC_folder_B = Path(BAC_folder).joinpath(tile)
        log.debug("source folder backscatter {}".format(SRC_folder_B))

        MOS_folder_tile = Path(MOS_folder).joinpath(tile)
        if not os.path.exists(MOS_folder_tile):
            os.makedirs(MOS_folder_tile)

        MOS_folder_tile_dir = Path(MOS_folder_tile).joinpath(ORBIT_DIRECTION)
        os.makedirs(MOS_folder_tile_dir, exist_ok=True)   
        ### CHECK IF ANY FILES EXIST ALREADY IN THE MOS_FOLDER AND RECORD THE BIGGEST FILE SIZE DIMENSIONS
        for rts, ds, files in os.walk(MOS_folder_tile_dir):
            if files:
                so_far = 0
                Base_file = ""
                for existing_file in files:
                    size = os.path.getsize(os.path.join(rts,existing_file))
                    if size > so_far:
                        so_far = size
                        Base_file = existing_file
                raster = gdal.Open(os.path.join(rts,Base_file))
                width = raster.RasterXSize
                height = raster.RasterYSize
                log.debug("Width and height found ... {} and {}".format(str(width), str(height)))

        Text_file = open(Path(MOS_folder_tile).joinpath(tile + '.txt'), 'w+')
        Text_file.write('FILE_NAME\n')
        # get the list of orbit track numbers of sentinel 1 satellites for given area of interest


        # create list of backscatter scenes with the current orbit direction

        # go through all the date in the date list and create mosaics from the backscatter scenes which belong to the
        # date interval.
        size_count = -1
        for date_instance in days_interval:
            log.debug(HLINE)
            log.debug(HLINE)

            # the backscatter scenes which belong to a particular date interval are first clipped to the given area of
            # interest and stored in a temporary directory. Then, these tiff files are then merged together using saga.

            # create temp directory to stored the clipped backscatter files that belong to the date interval.
            temp_folder = Path(work_dir).joinpath('temp', tile + '_' + str(date_instance))
            if not os.path.exists(temp_folder):
                os.makedirs(temp_folder)
            week_date = datetime.strptime(str(date_instance), '%Y-%m-%d')

            # name of the mosaic files that will be created for the current date.
            if '_BAC_' in TYPE:
                globals()['OUT_BAC_VH'] = str(Path(MOS_folder).joinpath(tile, ORBIT_DIRECTION, tile + '_BAC_VH_' + str(week_date).split(' ')[0].replace('-', '') + '.tiff'))
                globals()['OUT_BAC_VV'] = str(Path(MOS_folder).joinpath(tile, ORBIT_DIRECTION, tile + '_BAC_VV_' + str(week_date).split(' ')[0].replace('-', '') + '.tiff'))
            for pols in ['VV', 'VH']:
                for type in TYPE:
                    globals()[tile + '_List' + type + pols + '_' + str(date_instance).replace('-', '')] = list()
            log.debug("CURRENTLY PROCESSING  for date {}".format(str(date_instance)))
            log.debug("mosaic files are ... {} {}".format(OUT_BAC_VH, OUT_BAC_VV))







            # go through the available backscatter files to indentify those files which belong to the time interval.
            for root, dirs, files in os.walk(SRC_folder_B):
                for file in files:
                    if not os.path.exists(OUT_BAC_VH) or not os.path.exists(OUT_BAC_VV):
                        name_l = file.lower()
                        if not (name_l.endswith('.tif') or name_l.endswith('.tiff')):
                            continue
            
                        # find the observation datetime of the backscatter files

                        # --- robust extraction of date and polarization from filename (fixes 'm' before assignment) ---
                        m = re.search(r'(\d{8})(?:[tT]\d{6})?', file)
                        if m:
                            acq_date = datetime.strptime(m.group(1), '%Y%m%d')
                        else:
                            m = re.search(r'PC_(\d{8})', file)
                            if not m:
                                continue
                            acq_date = datetime.strptime(m.group(1), '%Y%m%d')

                        pm = re.search(r'(^|[_\.\-])(VV|VH)([_\.\-]|$)', file, re.IGNORECASE)
                        if not pm:
                            continue
                        pol = pm.group(2).upper()
                        # --- end of parser ---

                        # check if this file falls within the current 12 day interval
                        if (acq_date >= week_date) and (acq_date < week_date + timedelta(days=acq_frequency)):
                            # clip this file to the area of interest and store it in the temp folder.
                            temp_file = str(Path(temp_folder).joinpath(tile + '_' + file.split('.tiff')[0] + '.tiff'))
                            log.debug("temp file name ... {}".format(temp_file))
                            if not os.path.exists(temp_file) and not os.path.exists(
                                Path(MOS_folder).joinpath(
                                    tile, ORBIT_DIRECTION, tile + '_BAC_' + pol + '_' + str(week_date).split(' ')[0].replace('-', '') + '.tiff'
                                )
                            ):
                                cmd_gdal = [
                                    'gdal_translate',
                                    '-eco', '-projwin',
                                    '{}'.format(ulx), '{}'.format(uly), '{}'.format(lrx), '{}'.format(lry),
                                    '{}'.format(os.path.join(root, file)),
                                    '{}'.format(temp_file)
                                ]

                                cmd_output = subprocess.run(cmd_gdal, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                                log.debug("--exit code {} --> {}".format(cmd_output.returncode, cmd_gdal))

                            # store this temp file name in a list
                            (globals()[tile + '_List_BAC_' + pol + '_' + str(date_instance).replace('-', '')]).append(temp_file)


            # In the temp folder, files with particular polarisation are merged together.
            log.debug("----- ALL ORBITS WILL BE MOSAICED -----".format(str(date_instance)))
            for pols in ['VV', 'VH']:
                for type in TYPE:
                    if not os.path.exists(globals()[('OUT' + type + pols)]):
                        if len(globals()[tile + '_List' + type + pols + '_' + str(date_instance).replace('-', '')]) > 0:
                            count = -1

                            SAGA_list_orbs = [';'.join(globals()[tile + '_List' + type + pols + '_' + str(date_instance).replace('-', '')])]
                            TARGET_OUT_GRID = [Path(MOS_folder).joinpath(tile, ORBIT_DIRECTION, tile + type + pols + '_' + str(week_date).split(' ')[0].replace('-', ''))]
                            for SAGA_list in SAGA_list_orbs:
                                count += 1
                                size_count +=1
                                if (len(SAGA_list.split(';')) > 0) and not os.path.exists(str(TARGET_OUT_GRID[count]) + '.tiff'):

                                    cmd_gdal = ["saga_cmd",
                                                "-c=2",
                                                "grid_tools", str(3),
                                                "-GRIDS={}".format(SAGA_list),
                                                "-NAME=Mosaic",
                                                "-TYPE=9",
                                                "-RESAMPLING=0",
                                                "-OVERLAP=5",
                                                "-MATCH=0",
                                                "-TARGET_DEFINITION=0",
                                                "-TARGET_USER_SIZE=0.00017966",
                                                "-TARGET_USER_FITS=0",
                                                "-TARGET_OUT_GRID={}".format(TARGET_OUT_GRID[count])]
                                    cmd_output = subprocess.run(cmd_gdal, stdout=subprocess.PIPE,
                                                                stderr=subprocess.PIPE)
                                    log.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_gdal))
                                    if size_count == 0 and pols == 'VV' and not 'width' in globals() and not 'height' in globals():

                                        cmd_gdal = ['gdal_translate',
                                                    '-co', '{}'.format(' -co '.join(co)),
                                                    '-co', 'BIGTIFF=IF_NEEDED',
                                                    '{}'.format(str(TARGET_OUT_GRID[count]) + '.sdat'),
                                                    '{}'.format(str(TARGET_OUT_GRID[count]) + '.tiff')
                                                   ]
                                        cmd_output = subprocess.run(cmd_gdal, stdout=subprocess.PIPE,
                                                                    stderr=subprocess.PIPE)
                                        log.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_gdal))

                                        raster = gdal.Open(str(TARGET_OUT_GRID[count]) + '.tiff')
                                        globals()['width'] = raster.RasterXSize
                                        globals()['height'] = raster.RasterYSize
                                        log.debug("COUNT IS 0, width and height found ...width: {} height: {}".format(str(globals()['width']), str(globals()['height'])))
                                    else:
                                        if not 'width' in globals() and not 'height' in globals():
                                            cmd_gdal = ['gdal_translate',
                                                        '-co', '{}'.format(' -co '.join(co)),
                                                        '-co', 'BIGTIFF=IF_NEEDED',
                                                        '{}'.format(str(TARGET_OUT_GRID[count]) + '.sdat'),
                                                        '{}'.format(str(TARGET_OUT_GRID[count]) + '.tiff')
                                                        ]
                                            cmd_output = subprocess.run(cmd_gdal, stdout=subprocess.PIPE,
                                                                        stderr=subprocess.PIPE)
                                            log.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_gdal))

                                            raster = gdal.Open(str(TARGET_OUT_GRID[count]) + '.tiff')
                                            globals()['width'] = raster.RasterXSize
                                            globals()['height'] = raster.RasterYSize
                                            log.debug("COUNT IS {} width and height has now been found ...".format(str(size_count)))
                                        else:
                                            cmd_gdal = ['gdal_translate',
                                                        '-co', '{}'.format(' -co '.join(co)),
                                                        '-co', 'BIGTIFF=IF_NEEDED',
                                                        '-outsize','{}'.format(str(globals()['width'])), '{}'.format(str(globals()['height'])),
                                                        '{}'.format(str(TARGET_OUT_GRID[count]) + '.sdat'),
                                                        '{}'.format(str(TARGET_OUT_GRID[count]) + '.tiff')
                                                        ]
                                            cmd_output = subprocess.run(cmd_gdal, stdout=subprocess.PIPE,
                                                                        stderr=subprocess.PIPE)
                                            log.debug("exit code {} --> {}".format(cmd_output.returncode, cmd_gdal))

                                    #cmd_delete = 'find ' + str(Path(MOS_folder).joinpath(tile, ORBIT_DIRECTION)) + '/ -type f ! -iname "*.tiff"' + " -delete"
                                    #os.system(cmd_delete)
                                    for p in Path(MOS_folder).joinpath(tile, ORBIT_DIRECTION).glob('**/*'):
                                        if p.is_file() and p.suffix.lower() not in ('.tif', '.tiff'):
                                            try:
                                                p.unlink()
                                            except Exception:
                                                pass

                                    # Write only the output mosaic file name to the text file
                                    Text_string = [os.path.basename(str(TARGET_OUT_GRID[count]))]
                                    Text_file.write(','.join(str(x) for x in Text_string))
                                    Text_file.write('\n')


            shutil.rmtree(temp_folder)
        Text_file.close()
    log.info("SAR mosaic job has been finished.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Backscatter mosaic script.")
    parser.add_argument("--in_tile", type=str, help="tile name.", required=True, default=None)
    parser.add_argument("--backscatter_folder", type=str, help="year.", required=True, default=None)
    parser.add_argument("--mosaic_folder", type=str, help="Mosaic directory inside the output directory.", required=True, default=None)
    parser.add_argument("--aux_folder", type=str, help="Auxiliary data directory.", required=True, default=None)
    parser.add_argument("--in_multipolygon_string", type=str, help="Area of interest in multipolygon/polygon string.", required=True, default=None)
    parser.add_argument("--in_timeperiod", type=str, help="Time period.", required=True, default=None)

    args = parser.parse_args()

    run(tile=args.in_tile,
	    backscatter_folder = args.backscatter_folder,
	    mosaic_folder=args.mosaic_folder,
        aux_folder=args.aux_folder,
	    multipolygon_string=args.in_multipolygon_string,
	    timeperiod=args.in_timeperiod
	    )
