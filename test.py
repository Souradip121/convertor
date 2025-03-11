import numpy as np
import h5py as h5
import sys
import os
from osgeo import gdal
from osgeo import osr
import os
import pyproj
from watchdog.observers.polling import PollingObserver
from watchdog.events import FileSystemEventHandler

import time
from datetime import datetime
#from geos_gtif import create_geotif
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFaature
import cartopy
import cartopy.crs as ccrs
import configparser
import json
import re

logs_dir = './logs'

def write_log(rule, message, level='INFO'):
    log_file_name = get_dynamic_log_file_name(rule)
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_entry = f"{timestamp} - {level} - {rule} - {message}\n"
    # log_file_path = os.path.join(logs_dir, log_file_name)
    with open(log_file_name, 'a') as log_file:
        log_file.write(log_entry)
    log_file.close()

def get_dynamic_log_file_name(rule='COGGEN'):
    now = datetime.now()
    return os.path.join(logs_dir, f"{rule}_log_{now.strftime('%Y%m%d')}.log")

class EventHandler(FileSystemEventHandler):

    #product_code = ['DHI', 'DNI', 'GHI', 'INS', 'GHI_DLY', 'DNI_DLY', 'DHI_DLY', 'INS_DLY', 'HEM_DLY', 'AOD']

    Months = {'JAN': 1, 'FEB': 2, 'MAR': 3, 'APR': 4, 'MAY': 5, 'JUN': 6, 'JUL': 7, 'AUG': 8, 'SEP': 9, 'OCT': 10, 'NOV': 11, 'DEC': 12}

    #lpp_code = ['T2', 'T2S1', 'T2S2', 'T2S3', 'T2S1P1', 'T2S1P2', 'T2S1P3', 'T2S1P4', 'T2S1P5']
    #product_path = {
    #    'T2S1': 'insat3dr_dhi_daily_T2S1P1', 'T2S2': 'insat3dr_dhi_daily_T2S1P2', 'T2S3': 'insat3dr_dhi_daily_T2S1P3',
    #    'T2S1P1': 'insat3dr_dhi_daily_T2S1P1', 'T2S1P2': 'insat3dr_dhi_daily_T2S1P2', 'T2S1P3': 'insat3dr_dhi_daily_T2S1P3',
    #    'T2S1P4': 'insat3dr_dhi_daily_T2S1P4', 'T2S1P5': 'insat3dr_dhi_daily_T2S1P5'
    #}
    #sat_code = ['3RIMG', '3DIMG', '3SIMG']

    def __init__(self, outputdir, product_code, lpp_code, sat_code):
        self.outputdir = outputdir
        self.product_code = product_code
        self.lpp_code = lpp_code
        self.sat_code = sat_code

    def create_directory(self, filename, scode):
        pattern = r"_(\d{2})([A-Z]{3})(\d{4})_(\d{4})_"  # 29MAR2024_0000
        pattern2 = r"^([^_]+_){4}"  # 3DIMG_29MAR2024_0000_L3B_
        pattern3 = r"^([^_]+){1}"  # 3DIMG
        pattern4 = r"V(\d{2})R(\d{2})"  # V01R00

        mtch = re.search(pattern, filename)
        mtch_2 = re.search(pattern2, filename)
        mtch_3 = re.search(pattern3, filename)

        mtch_4 = re.search(pattern4, filename)

        if 'LIC' in filename:
            sector_name = re.search(r'LIC_(.*?)_V', filename).group(1)

        dirpath = None
        if mtch :
            DD, MMM, YYYY, HHMM = mtch.groups()
            dt = datetime.strptime(str.format('%s%s%s%s'%(DD,MMM,YYYY,HHMM)), '%d%b%Y%H%M')

        if mtch_3.group(0) in self.sat_code:
            dirpath = os.path.join(self.outputdir, mtch_3.group(0), str(YYYY), '%2d'%(self.Months[MMM]), '%2d'%(int(DD)))

        os.makedirs(dirpath, exist_ok=True)
        tmstr = str.format('%s%.2d%.2d%sP%s'%(YYYY,self.Months[MMM],int(DD),HHMM,HHMM))
        tfilename = f"{mtch_2.group(0)}{sector_name}_{scode}_{mtch_4.group(0)}.tif"

        if 'LIC' in filename:
            tfilename = f"{mtch_2.group(0)}{sector_name}_{scode}_{mtch_4.group(0)}.tif"

        return dirpath, tfilename

    def on_created(self, event):
        print(f"Watchdog received created event - {str(event.src_path)}")
        fl_basename = os.path.basename(str(event.src_path))
        fl, ext = os.path.splitext(fl_basename)

        if ext != '.h5':
            return

        last_size = -1
        while True:
            csize = os.path.getsize(str(event.src_path))
            if csize == last_size:
                break
            last_size = csize
            time.sleep(2)

        datasets = self.getdatasetslist(str(event.src_path))
        if len(datasets) == 0:
            return

        plevel = self.getprocessinglevel(str(event.src_path))

        for code in self.product_code:
            if code in datasets:
                dirpath, tfilename = self.create_directory(fl, code)
                tfl = os.path.join(dirpath, tfilename)

                if code == 'HEM_DLY' or code == 'LST':
                    if plevel in ['L1B', 'L2B', 'L3B']:
                        write_log('COGGEN', f"Watchdog received moved event - {event.src_path}", level='INFO')

                        create_geotiff(str(event.src_path), code, tfl)
                        write_log('COGGEN', f"Successfully written: {tfl}", level='INFO')

                elif plevel in ['L2', 'L3G']:
                    write_log('COGGEN', f"Watchdog received moved event - {event.src_path}", level='INFO')

                    create_gridded_geotiff(str(event.src_path), code, tfl)
                    write_log('COGGEN', f"Successfully written: {tfl}", level='INFO')

                else:
                    write_log('COGGEN', f"Watchdog received moved event - {event.src_path}", level='INFO')

                    write_toa_geotiff(str(event.src_path), code, tfl)
                    write_log('COGGEN', f"Successfully written: {tfl}", level='INFO')

    def on_moved(self, event):
        fl_basename = os.path.basename(str(event.src_path))
        fl, ext = os.path.splitext(fl_basename)


        datasets = self.getdatasetslist(str(event.src_path))

        for code in self.product_code:
            if code in datasets:
                dirpath, tfilename = self.create_directory(fl, code)
                tfl = os.path.join(dirpath, tfilename)
                write_log('COGGEN', f"Watchdog received moved event - {event.src_path}", level='INFO')

                create_geotiff(str(event.src_path), code, tfl)
                write_log('COGGEN', f"Successfully written: {tfl}", level='INFO')

    def getdatasetslist(self, filename):
        filename = os.path.basename(filename)
        datasets = []
        for code in self.lpp_code:
            if code in filename:
                ds = h5.File(filename, 'r')
                datasets = list(ds.keys())
                ds.close()
        return datasets

    def getprocessinglevel(self, filename):
        # filename = os.path.basename(filename)
        file_attrs = {}
        ds = h5.File(filename, 'r')
        file_attrs = ds.attrs
        ds.close()
        return file_attrs['Processing_Level'].astype(str)

def get_corner_coordinates(llc_hdf, dataset):
    write_log('COGGEN', f"Reading {llc_hdf} dataset: {dataset}", 'INFO')

    band_attrs = {}

    band = None
    try:
        dst_h5 = h5.File(llc_hdf, 'r')
        band = dst_h5.get(dataset)
        band_attrs.update(band.attrs)
        proj_info = dst_h5.get('Projection_Information')
        dst_proj = dst_h5.get(proj_info.attrs['projection'])
        dst_h5.close()
    except Exception as e:
        write_log('COGGEN', f"Exception in opening {llc_hdf}-{str(e)}", 'ERROR')

    return dst_proj, band_attrs, band

def convert_lonlat_to_xy(proj_param, step=0):
    clat = proj_param['standard_parallel'][0]
    clon = proj_param['longitude_of_projection_origin'][0]
    fe = proj_param['false_easting'][0]
    fn = proj_param['false_northing'][0]
    grid_mapping_name = proj_param['grid_mapping_name']

    if step == 0:
        lat = proj_param['upper_left_lat_lon(degrees)']
        lon = proj_param['upper_left_lat_lon(degrees)']
    elif step == 1:
        lat = proj_param['lower_right_lat_lon(degrees)']
        lon = proj_param['lower_right_lat_lon(degrees)']
        write_log('COGGEN', f"Lower right corner: Latitude {lat}, Longitude {lon}", level='INFO')

    write_log('COGGEN', f"Central Latitude: {clat}, Central Longitude: {clon}", level='INFO')
    bmap = pyproj.Proj(proj="merc", lat_ts=clat, lon_0=clon, ellps="WGS84")
    crs = pyproj.CRS.from_proj4(bmap.to_proj4())

    return bmap(lon, lat, inverse=False), crs.to_wkt()

def write_toa_geotiff(l3c_hdf_filename, dataset, l3c_geotiff_filename):
    proj_param, band_attrs, band = get_corner_coordinates(l3c_hdf_filename, dataset)


    print(f"({len(proj_param.keys())}) ({len(band_attrs.keys())})");
    if len(proj_param.keys()) == 0 or len(band_attrs.keys()) == 0:
        write_log('COGGEN', f"File {l3c_hdf_filename} either corrupted or (dataset) missing", 'ERROR')
        return

    dst_ul_xy, wkt = convert_lonlat_to_xy(proj_param)
    dst_lr_xy, _ = convert_lonlat_to_xy(proj_param, 1)

    write_log('COGGEN', f"Upper-left corner (x, y): ({dst_ul_xy[0]}, {dst_ul_xy[1]})", 'INFO')
    write_log('COGGEN', f"Lower-right corner (x, y): ({dst_lr_xy[0]}, {dst_lr_xy[1]})", 'INFO')
    write_log('COGGEN', f"Band shape: {band.shape}", 'INFO')

    xres = pixel_resolution = dst_lr_xy[0] - dst_ul_xy[0]
    yres = npixel_resolution = dst_ul_xy[1] - dst_lr_xy[1]

    nscan, npix = band.shape

    xres = xres / npix
    yres = yres / nscan

    write_log('COGGEN', f"xres: {xres}, yres: {yres}", level='INFO')
    write_log('COGGEN', f"nscan: {nscan}, npix: {npix}", level='INFO')

    dst_ds = gdal.GetDriverByName('GTiff').Create(l3c_geotiff_filename, npix, nscan, 1, gdal.GDT_Float32)
    dst_ds.SetGeoTransform([int(dst_ul_xy[0]), xres, 0.0, dst_ul_xy[1], 0.0, yres])

    bl = dst_ds.GetRasterBand(1)
    dset_attr = dset.attrs
    NoDataValue = dset_attr['_FillValue'][0]

    dset = np.squeeze(dset[()])
    dset = dset[7:2759, 47:2757]
    file.close()

    if chain in [1, 6, 5]:
        globe = ccrs.Globe(semimajor_axis=6378137.0, semiminor_axis=6356752.31414)
        chain_key = "3SIMG"
        proj = ccrs.Geostationary(central_longitude=-74.0, satellite_height=35785831, false_easting=0, false_northing=0, sweep_axis='x', globe=globe)
        if chain == 35IMG:
            proj = ccrs.Geostationary(central_longitude=-74.0, satellite_height=35785831, false_easting=0, false_northing=0, sweep_axis='x', globe=globe)

        geos = osr.SpatialReference()
        geos.SetWellKnownGeogCS('WGS84')
        geos.SetProjCS("GEOS")
        wkt = geos.ExportToWkt()

        drv = gdal.GetDriverByName('MEM')
        dst_ds = drv.Create('', dset.shape[1], dset.shape[0], 1, gdal.GDT_Float32)
        dst_ds.SetProjection(wkt)

        xres = (proj.x_limits[1] - proj.x_limits[0]) / dset.shape[1]
        yres = (proj.y_limits[1] - proj.y_limits[0]) / dset.shape[0]

        dst_ds.SetGeoTransform([proj.x_limits[0], xres, 0.0, proj.y_limits[0], 0.0, -yres])

        bl = dst_ds.GetRasterBand(1)
        bl.WriteArray(dset)
        bl.SetNoDataValue(int(NoDataValue))

    temp_dir = '/tmp'
    inter_gtif = os.path.join(temp_dir, '%s_OUT.tif' % (inpfl))

    print(f"inter_gtif) {geotiff_filename}")

    ds = gdal.Translate(inter_gtif, dst_ds, outputBounds = [proj.x_limits[0], proj.y_limits[0], proj.x_limits[1], proj.y_limits[1]], resampleAlg="near", noData=NoDataValue)

    ws = gdal.Warp(geotiff_filename, ds, dstSRS="EPSG:4326", xRes=0.03636, yRes=0.03636, resampleAlg="near", srcNoData=NoDataValue, dstNoData=NoDataValue)
    ws = None
    ds = None
    os.remove(inter_gtif)
    
    # Create COG after GeoTIFF
    # create_cog(geotiff_filename)

def create_gridded_geotiff(hdf_filename, dataset, geotiff_filename):
    flt = h5.File(hdf_filename)
    fl_attrs = {}

    fl_attrs.update(flt.attrs)
    ds = flt.get(dataset)
    ds_attrs = ds.attrs
    ds = ds[()]

    ds_attrs.update(fl.get(dataset).attrs)
    flt.close()

    ds = gdal.GetDriverByName('GTiff')

    dst_ds = drv.Create(geotiff_filename, ds.shape[1], ds.shape[0], 1, gdal.GDT_Float32)

    dst_ds.SetProjection('EPSG:4326')
    dst_ds.SetGeoTransform([fl_attrs['left_longitude'][0], fl_attrs['lon_interval'][0], 0.0, fl_attrs['upper_latitude'][0], 0.0, fl_attrs['lat_interval'][0]])

    bl = dst_ds.GetRasterBand(1)
    bl.WriteArray(ds)
    bl.SetNoDataValue(int(ds_attrs['_FillValue'][0]))

    dst_ds = None
    
    # Create COG after GeoTIFF
    create_cog(geotiff_filename)

def create_cog(input_tiff, output_cog=None):
    """
    Convert a regular GeoTIFF to a Cloud Optimized GeoTIFF (COG) with LZW compression
    and overview levels up to 6.
    
    Args:
        input_tiff: Path to the input GeoTIFF file
        output_cog: Path to the output COG file. If None, overwrite the input file
    """
    if output_cog is None:
        output_cog = input_tiff
        temp_tiff = input_tiff + '.temp.tif'
        os.rename(input_tiff, temp_tiff)
        input_tiff = temp_tiff

    try:
        write_log('COGGEN', f"Creating COG for {input_tiff}", level='INFO')
        
        # Create overview levels
        ds = gdal.Open(input_tiff, gdal.GA_ReadOnly)
        overview_levels = [2, 4, 8, 16, 32, 64]
        gdal.SetConfigOption('COMPRESS_OVERVIEW', 'LZW')
        ds.BuildOverviews("NEAREST", overview_levels)
        ds = None
        
        # Create COG with LZW compression
        translate_options = gdal.TranslateOptions(
            format='GTiff',
            creationOptions=[
                'COMPRESS=LZW',
                'TILED=YES',
                'COPY_SRC_OVERVIEWS=YES',
                'BIGTIFF=IF_NEEDED',
                'NUM_THREADS=ALL_CPUS'
            ]
        )
        
        gdal.Translate(output_cog, input_tiff, options=translate_options)
        write_log('COGGEN', f"Successfully created COG: {output_cog}", level='INFO')
        
        if input_tiff.endswith('.temp.tif'):
            os.remove(input_tiff)
            
    except Exception as e:
        write_log('COGGEN', f"Failed to create COG: {str(e)}", level='ERROR')
        if input_tiff.endswith('.temp.tif'):
            os.rename(input_tiff, output_cog)

def main():
    if (len(sys.argv) < 2):
        print(f'Usage {sys.argv[0]} <config.ini>')
        sys.exit(1)

    config = configparser.ConfigParser()
    config.read(sys.argv[1])

    watch_dirs = config['COG_GENERATION']['WATCH_DIR'].split(',')  # All but the last argument are input directories
    output_dir = config['COG_GENERATION']['OUTPUT_DIR']  # The last argument is the output directory
    sat_codes = config['COG_GENERATION']['SAT_CODES'].split(',')
    lpp_codes = config['COG_GENERATION']['LPP_CODES'].split(',')
    product_codes = config['COG_GENERATION']['PRODUCT_CODES'].split(',')
    
    # Set the global logs_dir variable
    global logs_dir
    logs_dir = config['COG_GENERATION']['LOGS']
    
    # Create logs directory if it doesn't exist
    os.makedirs(logs_dir, exist_ok=True)

    event_handler = EventHandler(output_dir, product_codes, lpp_codes, sat_codes)

    observer = PollingObserver()
    for watch_dir in watch_dirs:
        observer.schedule(event_handler, watch_dir, recursive=True)
    
    print(f"Monitoring directories: {watch_dirs}")
    write_log('COGGEN', f"Monitoring started for directories: {watch_dirs}", level='INFO')

    try:
        while True:
            time.sleep(0.5)
    except KeyboardInterrupt:
        observer.stop()
        observer.join()
        write_log('COGGEN', "Monitoring stopped by user", level='INFO')

if __name__ == "__main__":
    main()



