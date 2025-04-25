# basic packages
import os
import glob
import urllib
from numpy import newaxis
import xarray as xr
import numpy as np
import pandas as pd
import warnings

# rasterio
import rasterio
import rasterio.plot
from rasterio import features
from rasterio.windows import from_bounds
from rasterio.plot import show
from rasterio.enums import Resampling
from rasterio.mask import mask
from rasterio.fill import fillnodata
from rasterio.warp import reproject, Resampling, calculate_default_transform

# gis
from pysheds.grid import Grid
from scipy import ndimage
import geopandas as gpd
from rasterio.crs import CRS

# plotting
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import colors


def vmi_from_puhti(fd, subset, out_fd, layer='all', interpolate=None, use_center=None, max_search_distance=2, resample=False, plot=True, save_in='geotiff'):
    '''
    Downloads VMI data (all layers 'all' or one given layer e.g. 'kasvupaikka_vmi1x_1721.img') 
    from given file directory (fd) for a given subset area
    Option to interpolate over the VMI classified 'nonland' grid-cells with a given max_search_distance (1)
    Option to assign the 'nonland' grid-cells as zeros and then interpolate only those values (2)
    Option to resample the original 16m resolution rasters with a given scaling factor (e.g. 0.5 -> 32m)
    Saving in output file directory (out_fd) as 'geotiff' or 'ascii'
    '''
    
    if layer=='all':
        p = os.path.join(fd, '*.img')
    else:
        p = os.path.join(fd, layer)

    if not os.path.exists(out_fd):
        os.makedirs(out_fd)

    # file count so that mask will be only saved simulataneously with first file
    i = 0
    for file in glob.glob(p):
        print(file)
        out_fn = file.rpartition('/')[-1][:-4]
        mask_fn = 'nonland'
        if save_in == 'geotiff':
            out_fp = os.path.join(out_fd, out_fn) + '.tif'
            mask_fp = os.path.join(out_fd, mask_fn) + '.tif'
        elif save_in == 'asc':
            out_fp = os.path.join(out_fd, out_fn) + '.asc'
            mask_fp = os.path.join(out_fd, mask_fn) + '.asc'

        with rasterio.open(file) as src:
            data = src.read(1, window=from_bounds(subset[0], subset[1], subset[2], subset[3], src.transform))
            profile = src.profile
            out_meta = src.profile.copy()

            new_affine = rasterio.Affine(out_meta['transform'][0], 
                                         out_meta['transform'][1], 
                                         subset[0], 
                                         out_meta['transform'][3], 
                                         out_meta['transform'][4], 
                                         subset[3])

            # save nonland mask at first loop iteration
            if i == 0:
                data_mask = np.zeros(shape=data.shape)
                data_mask[data == 32767] = 32767
            
            if len(data.flatten()[data.flatten() == 32766]) > 0:
                print('*** Data has', len(data.flatten()[data.flatten() == 32766]), 'nan values (=32766) ***')
                #print('--> converted to 0 ***')
            if len(data.flatten()[data.flatten() == 32767]) > 0:
                print('*** Data has', len(data.flatten()[data.flatten() == 32767]), 'non land values (=32767) ***')
                #print('--> converted to 0 ***')

            # Update the metadata for geotiff
            out_meta.update({"driver": "GTiff",
                             "height": data.shape[0],
                             "width": data.shape[1],
                              "transform": new_affine,
                              "nodata": 32767,
                              "crs": CRS.from_epsg(3067)})
            
            if save_in == 'asc':
                out_meta.update({"driver": "AAIGrid"})
                
            with rasterio.Env():
                with rasterio.open(out_fp, 'w', **out_meta, force_cellsize=True) as dst:
                    src = dst.write(data, 1)
                            # if save_mask:
                if i == 0:
                    with rasterio.open(mask_fp, 'w', **out_meta, force_cellsize=True) as dst_mask:
                        src_mask = dst_mask.write(data_mask, 1)  
                    
            if interpolate == 1:
                with rasterio.open(out_fp, 'r+') as src_new:
                        data = src_new.read(1)
                        data_filled = fillnodata(data, 
                                                 mask=src_new.read_masks(1), 
                                                 max_search_distance=max_search_distance, 
                                                 smoothing_iterations=0)
                        print('*** non land interpolated ***')
                if i == 0:
                    with rasterio.open(mask_fp) as src_mask_interp:
                        mask_interp = src_mask_interp.read(1)
                        mask_interp_fill = fillnodata(mask_interp, 
                                                      mask=src_mask_interp.read_masks(1), 
                                                      max_search_distance=max_search_distance, smoothing_iterations=0)
                        
                with rasterio.Env():
                    with rasterio.open(out_fp, 'w', **out_meta, force_cellsize=True) as dst_new:
                        src_new = dst_new.write(data_filled, 1)
                    if i == 0:
                        with rasterio.open(mask_fp, 'w', **out_meta, force_cellsize=True) as dst_mi:
                            src_mi = dst_mi.write(mask_interp_fill, 1)

            elif interpolate == 2:
                with rasterio.open(out_fp, 'r+') as src_new:
                        data = src_new.read(1)
                        data_filled = interpolate_over_mask(data=data, mask=data_mask, use_center=use_center)
        
                with rasterio.Env():
                    with rasterio.open(out_fp, 'w', **out_meta, force_cellsize=False) as dst_new:
                        src_new = dst_new.write(data_filled, 1)         

            elif interpolate == 3: # does not interpolate but assigns non_land as zeros
                with rasterio.open(out_fp, 'r+') as src_new:
                        data = src_new.read(1)
                        data[data_mask == 32767] = 0.0
        
                with rasterio.Env():
                    with rasterio.open(out_fp, 'w', **out_meta, force_cellsize=False) as dst_new:
                        src_new = dst_new.write(data, 1)    
                        
            if resample:
                with rasterio.open(out_fp) as src_new:
                    # resample data to target shape
                    print(f'*** resampled with a scalefactor of {resample} ***')
                    new_width = int((subset[2]-subset[0])/(16/resample))
                    new_height = int((subset[3]-subset[1])/(16/resample))
                    data_resampled = src_new.read(1,
                                    out_shape=(src_new.count,int(new_height),int(new_width)),
                                    resampling=Resampling.bilinear
                                    )
                      
                    out_meta = src_new.profile.copy()

                    new_affine = rasterio.Affine(out_meta['transform'][0]/resample, 
                                                 out_meta['transform'][1], 
                                                 subset[0], 
                                                 out_meta['transform'][3], 
                                                 out_meta['transform'][4]/resample, 
                                                 subset[3])
                    # Update the metadata for geotiff
                    out_meta.update({"height": data_resampled.shape[0],
                                    "width": data_resampled.shape[1],
                                    "transform": new_affine})

                if i == 0:
                    with rasterio.open(mask_fp) as src_mask_new:
                        data_mask_resampled = src_mask_new.read(1,
                                        out_shape=(src_mask_new.count,int(new_height),int(new_width)),
                                        resampling=Resampling.bilinear
                                        )
                with rasterio.Env():
                    with rasterio.open(out_fp, 'w', **out_meta, force_cellsize=True) as dst:
                        src_new = dst.write(data_resampled, 1)
                    if i == 0:
                        with rasterio.open(mask_fp, 'w', **out_meta, force_cellsize=True) as dst_mask:
                            src_mask_new = dst_mask.write(data_mask_resampled, 1)

        i += 1
                                
        if plot==True:
            raster = rasterio.open(out_fp)
            show(raster)
