from pathlib import Path
import iris
import matplotlib.pyplot as plt
import numpy as np

import warnings

from esmvalcore.preprocessor import climate_statistics,extract_volume,concatenate, axis_statistics, extract_region, area_statistics,extract_shape,add_supplementary_variables

#variables = ["tottrdsdetoc","mort2sdetoc","zmicronoassim","ldetoc2sdetoc","doc2sdetoc",\
#             "zmesofragmldetoc","zmicroingestsdetoc","zmesoingestffsdetoc","remsdetoc",\
#             "sdetoc2ldetocagg","xadsdetoc","yadsdetoc","zadsdetoc","ldfsdetoc",\
#             "zdfsdetoc","trdexpsdetoc"]

variables = ["zingest2doc"]

diag_list= ['a67o'] #,'a683','a5xl','a5kq']

shapefile = '/esarchive/scratch/aorihuel/mask/change_mask/shapefiles_py/all_regions_1shapefile/IR_LAB_budgets'

rowName_List = ["aLAB_IR"]


ZMIN_LIST=[0, 100, 500, 100, 1000]
ZMAX_LIST=[100, 500, 1000, 1000, 2000]
layer_list = ['epipelagic', 'upper_meso','lower_meso', 'meso','bathy']

areacello = iris.load_cube('/esarchive/exp/nemo/a5gj/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Ofx/areacello/gn/v20221216/areacello_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc')

bathymask = iris.load_cube('/esarchive/scratch/aorihuel/bathy_mask/bathymask_2500_0_1.nc')
mask_data = np.ma.masked_where(bathymask.data == 0, bathymask.data)
mask_cube = bathymask.copy(mask_data)

cmorversion='v*'
yrStrt=2009
yrLast=2019
type = 'pkt' # 'small_POC_budgets'


data_region = {}
data_region_mean={}

for expid in diag_list:
    clim_region_cubes = []
    clim_area_mean = []
    clim_cubes = []
    for var in variables:
        cubes = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for yr in range(yrStrt,yrLast):
                # 1) Cargamos los ficheros para calcular las climatologias
                cube = iris.load_cube("/esarchive/scratch/aorihuel/BUDGETS_LB_IR/non_clim/"+expid+"/" + var + "/" + var +"_Omon_EC-Earth3_opera-control_r1i1p1f1_gn_"+ str(yr) + "01-" + str(yr) + "12.nc")
           
                # 2) Maskeamos los puntos shallower a 1200m depth
                cube_masked_data = cube.data*mask_cube.data
                bathy_cube = cube.copy(cube_masked_data)
                # 3) maskeamos con los shapefiles las regiones que nos interesan:
                region_cube = extract_shape(bathy_cube, shapefile, method='contains', crop=False, decomposed=True, ids=rowName_List)
                clim_region_cubes.append(region_cube)
            
                data_region[expid] = clim_region_cubes
                # 4) Hacemos la space average
                add_supplementary_variables(region_cube, [areacello])
                var_mean = area_statistics(region_cube, 'mean')
            
                clim_area_mean.append(var_mean)
            
                data_region_mean[expid] = clim_area_mean
                filo2 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/'+var+'_'+expid+'_sp_average_IR_LAB.nc')
                filo2.parent.mkdir(parents=True, exist_ok=True)
                iris.save(var_mean, filo2)
            # 5) Hacemos la integracion vertical
                for layer in layer_list:
                    if layer ==layer_list[0]:
                        zmin=ZMIN_LIST[0]
                        zmax=ZMAX_LIST[0]
                        var_layer_cube = extract_volume(var_mean, z_min=zmin, z_max=zmax)
                        var_layer_mean = axis_statistics(var_layer_cube, 'z', 'sum')
                         # 6) Guardamos
                        filo3 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/bathy_depth_int/'+layer+'/'+var+'_'+expid+'_sp_average_IR_LAB_'+str(zmin)+'_'+str(zmax)+'_'+layer+'.nc')
                        filo3.parent.mkdir(parents=True, exist_ok=True)
                        iris.save(var_layer_mean, filo3)
                    if layer ==layer_list[1]:
                        zmin=ZMIN_LIST[1]
                        zmax=ZMAX_LIST[1]
                        var_layer_cube = extract_volume(var_mean, z_min=zmin, z_max=zmax)
                        var_layer_mean = axis_statistics(var_layer_cube, 'z', 'sum')
                     # 6) Guardamos
                        filo3 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/bathy_depth_int/'+layer+'/'+var+'_'+expid+'_sp_average_IR_LAB_'+str(zmin)+'_'+str(zmax)+'_'+layer+'.nc')
                        filo3.parent.mkdir(parents=True, exist_ok=True)
                        iris.save(var_layer_mean, filo3)
                    if layer ==layer_list[2]:
                        zmin=ZMIN_LIST[2]
                        zmax=ZMAX_LIST[2]
                        var_layer_cube = extract_volume(var_mean, z_min=zmin, z_max=zmax)
                        var_layer_mean = axis_statistics(var_layer_cube, 'z', 'sum')
                         # 6) Guardamos
                        filo3 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/bathy_depth_int/'+layer+'/'+var+'_'+expid+'_sp_average_IR_LAB_'+str(zmin)+'_'+str(zmax)+'_'+layer+'.nc')
                        filo3.parent.mkdir(parents=True, exist_ok=True)
                        iris.save(var_layer_mean, filo3)
                    if layer ==layer_list[3]:
                        zmin=ZMIN_LIST[3]
                        zmax=ZMAX_LIST[3]
                        var_layer_cube = extract_volume(var_mean, z_min=zmin, z_max=zmax)
                        var_layer_mean = axis_statistics(var_layer_cube, 'z', 'sum')
                         # 6) Guardamos
                        filo3 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/bathy_depth_int/'+layer+'/'+var+'_'+expid+'_sp_average_IR_LAB_'+str(zmin)+'_'+str(zmax)+'_'+layer+'.nc')
                        filo3.parent.mkdir(parents=True, exist_ok=True)
                        iris.save(var_layer_mean, filo3)
                    if layer ==layer_list[4]:
                        zmin=ZMIN_LIST[4]
                        zmax=ZMAX_LIST[4]
                        var_layer_cube = extract_volume(var_mean, z_min=zmin, z_max=zmax)
                        var_layer_mean = axis_statistics(var_layer_cube, 'z', 'sum')
                         # 6) Guardamos
                        filo3 = Path('/esarchive/scratch/aorihuel/BUDGETS_LB_IR/'+type+'/diag_'+expid+'/yr_files/'+ str(yr) +'/sp_average/'+rowName_List[0]+'/bathy_depth_int/'+layer+'/'+var+'_'+expid+'_sp_average_IR_LAB_'+str(zmin)+'_'+str(zmax)+'_'+layer+'.nc')
                        filo3.parent.mkdir(parents=True, exist_ok=True)
                        iris.save(var_layer_mean, filo3)
                    
        
        
                
            
            

            
            

            
            