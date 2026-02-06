How to use NEFCUP?

CLASSIFICATION SHARES:
- area fraction of each class within the 10x10km grid
1. classified_shares.ipynb
- classified mNFI grids to 12 classes based on maaluokka, paatyyppi, ravinteikkuus
- produces one grid with 1-12 number for each 16x16m grid-cell
2. classified_shares_to_10x10.pynb
- calculates fractions of each class for 10x10km grid
- produces 12 rasters containing each class fraction within the 10x10km grid

CLASSIFICATION IKA AND PPA
- ika and ppa grids for each class
3. mnfi_to_nc.ipynb
- saves the mnfi ika and ppa files to nc
4. classify_mnfi
- classified mNFI ppa and age grids to 12 classes based on maaluokka, paatyyppi, ravinteikkuus
- produces 12 rasters for ika and ppa each containing the given class ika and ppa

REPROJECTION AND RESAMPLING
- reprojecting and resampling to common crs and grid
5. resample_and_reproj_mnfi.ipynb
- Resamples and reprojects the 'classified' mnfi .tif files. Uses reference grid
6. reproj_khkmask.ipynb
- reprojects khkmask
7. reproj_fmi.ipynb
- reproject fmi grids

MERGE
- preparing to compute
8. merge_final_layers.ipynb
- combines mnfi layers to nc file
- creates the meteo file
- creates the regional mask

COMPUTING FLUXES
- computing and scaling
9. upscaling_to_fluxes.ipynb
- does the calculation for the 12 age/ppa classes
10. scale_fluxes.ipynb
- scales fluxes based on their fractions in grid-cells
11. regional_fluxes.ipynb
- uses the regional mask to compute fluxes
12. read_final_layers.ipynb
- plots the upscaled layers