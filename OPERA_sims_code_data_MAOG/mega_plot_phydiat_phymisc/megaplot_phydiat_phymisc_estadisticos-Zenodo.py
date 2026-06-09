# Mega plot con obs phydiat & phymisc separados + estadisticos

import glob
from pathlib import Path

import iris
from esmvalcore.preprocessor import load, concatenate, climate_statistics

import matplotlib.pyplot as plt

import iris.plot as iplt
import iris.quickplot as qplt

import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfea

import iris.analysis.cartography

import numpy as np

import matplotlib as mpl
import warnings

variables =["phymiscos"] #"phydiatos"] #, "phymiscos"]

yrStrt=2009
yrLast=2019

data_sat = {}

for yr in range(yrStrt,yrLast):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cubes = []
        for var in variables:
            cube = iris.load_cube("/esarchive/obs/esa/pc_oc-cci-v5/surface_restoring/ESA_CCI_v5/ORCA1L75/fc00/"+var+"_y"+str(yr)+".nc")[0:12,:,:]
            mask = iris.load_cube("/esarchive/obs/esa/pc_oc-cci-v5/surface_restoring/masks/ESA_CCI_v5/ORCA1L75/trc_ssr_mask_y"+str(yr)+".nc")[0:12,:,:]
            cube.data = np.ma.masked_where(mask.data == 0, cube.data)
            cubes.append(cube)
            data_sat[str(yr)] = cubes

            

variables = ["phymisc"] #["phydiat"] #, "phymisc"]
#variables =["chldiat","chlmisc"]
diag_list = ['a67o','a683']


yrStrt=2009
yrLast=2019


data_expid = {}

for expid in diag_list:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        data_yr = {}
        for yr in range(yrStrt,yrLast):
            cubes = []
            for var in variables:            
                cube = iris.load_cube("/esarchive/exp/nemo/"+expid+"/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Omon/"+var+"/gn/v*/"+var+"_Omon_EC-Earth3_opera-control_r1i1p1f1_gn_"+str(yr)+"01-"+str(yr)+"12.nc")[0:12,0,:,:]
                mask = iris.load_cube("/esarchive/obs/esa/pc_oc-cci-v5/surface_restoring/masks/ESA_CCI_v5/ORCA1L75/trc_ssr_mask_y"+str(yr)+".nc")[0:12,:,:]
                cube.data = np.ma.masked_where(mask.data == 0, cube.data)
                cubes.append(cube)
                data_yr[str(yr)] = cubes
        data_expid[expid] = data_yr
        
yr_list = list(range(yrStrt,yrLast))
# Create new dictionary structure for summed data
data_sim = {}

for expid in diag_list:       
        # Initialize list to store summed cubes for each year
   
    data_sim_list = []   
        # Sum cubes year by year
    for year_idx, yr in enumerate(yr_list):
        # Start with the first variable for this year
        
        # Add the remaining variables for this year
        sim_cube = data_expid[expid][str(yr)][0]
                #summed_cube.var_name = 'sPOC'
                #summed_cube.units = 'mol m-2'
        data_sim_list.append(sim_cube)
        # Store the list of summed cubes as 'sPOC'
    data_sim[expid]= concatenate(data_sim_list)

yr_list = list(range(yrStrt,yrLast))
# Create new dictionary structure for summed data
data_sat_list = []
# Sum cubes year by year
for year_idx, yr in enumerate(yr_list):       
    # Add the remaining variables for this year
    sat_cube = data_sat[str(yr)][0]
                
    data_sat_list.append(sat_cube)
        # Store the list of summed cubes as 'sPOC'
concat_data_sat = concatenate(data_sat_list)
data_sat_cube = data_sim[expid].copy(concat_data_sat.data)
data_sat_cube.units = 'mg C m-3'

from esmvalcore.preprocessor import extract_time, extract_month, climate_statistics
years_groups = {
    '2009-2013': ['2009', '2010', '2011', '2012', '2013'],
    '2014-2018': ['2014', '2015', '2016', '2017', '2018']
}
weak_list = ['2009', '2010', '2011', '2012', '2013']
strong_list = ['2014', '2015', '2016', '2017', '2018']
extracted_months = [3,4,5,6,7,8] #months to extract 1- based
mon_list = ['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug']
data_model_weak = {}
weak_mon_data={}
for expid in diag_list:
    weak_cube = extract_time(data_sim[expid], start_year=2009, start_month=1, start_day=1, end_year=2013, end_month=12, end_day=31)
    data_model_weak[expid]=weak_cube
    weak_mon_data[expid]={}
    for mon_idx, mon in enumerate(extracted_months):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mon_cube = extract_month(data_model_weak[expid], month=mon)
            mon_mean_cube = climate_statistics(mon_cube, operator='mean', period='full')
            weak_mon_data[expid][mon_list[mon_idx]]=mon_mean_cube

data_model_strong = {}
strong_mon_data={}
for expid in diag_list:
    strong_cube = extract_time(data_sim[expid], start_year=2014, start_month=1, start_day=1, end_year=2018, end_month=12, end_day=31)
    data_model_strong[expid]=strong_cube
    strong_mon_data[expid]={}
    for mon_idx, mon in enumerate(extracted_months):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mon_cube = extract_month(data_model_strong[expid], month=mon)
            mon_mean_cube = climate_statistics(mon_cube, operator='mean', period='full')
            strong_mon_data[expid][mon_list[mon_idx]]=mon_mean_cube

weak_sat_cube = extract_time(data_sat_cube, start_year=2009, start_month=1, start_day=1, end_year=2013, end_month=12, end_day=31)
strong_sat_cube = extract_time(data_sat_cube, start_year=2014, start_month=1, start_day=1, end_year=2018, end_month=12, end_day=31)
weak_sat_mon_data = {}
strong_sat_mon_data = {}
for mon_idx, mon in enumerate(extracted_months):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        weak_mon_cube_sat = extract_month(weak_sat_cube, month=mon)
        weak_mon_mean_cube_sat = climate_statistics(weak_mon_cube_sat, operator='mean', period='full')
        weak_sat_mon_data[mon_list[mon_idx]]=weak_mon_mean_cube_sat
        strong_mon_cube_sat = extract_month(strong_sat_cube, month=mon)
        strong_mon_mean_cube_sat = climate_statistics(strong_mon_cube_sat, operator='mean', period='full')
        strong_sat_mon_data[mon_list[mon_idx]]=strong_mon_mean_cube_sat

all_obs_dict ={'2009-2013':weak_sat_mon_data,
                  '2014-2018':strong_sat_mon_data}
ratio_weak = {}
ratio_strong = {}
for expid in diag_list:
    ratio_weak[expid]={}
    ratio_strong[expid]={}
    for mon_idx, mon in enumerate(extracted_months):
        weak_ratio_data =  weak_mon_data[expid][mon_list[mon_idx]].copy(weak_mon_data[expid][mon_list[mon_idx]].data*12*1000-(weak_sat_mon_data[mon_list[mon_idx]].data))
        weak_ratio_data.units = 'mg m-2'
        ratio_weak[expid][mon_list[mon_idx]]=weak_ratio_data
        strong_ratio_data = strong_mon_data[expid][mon_list[mon_idx]].copy(strong_mon_data[expid][mon_list[mon_idx]].data*12*1000-(strong_sat_mon_data[mon_list[mon_idx]].data))
        strong_ratio_data.units = 'mg m-2'
        ratio_strong[expid][mon_list[mon_idx]]=strong_ratio_data

all_ratios_dict ={'2009-2013':ratio_weak,
                  '2014-2018':ratio_strong}

bathy_a67o = iris.load_cube ('/esarchive/exp/nemo/a67o/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Ofx/deptho/gn/v20230630/deptho_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc')
bathy_a683 = iris.load_cube ('/esarchive/exp/nemo/a683/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Ofx/deptho/gn/v*/deptho_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc')

## Calculamos ahora los estadisticos
import numpy as np
import numpy.ma as ma
from esmvalcore.preprocessor import area_statistics, extract_shape, add_supplementary_variables
shapefile = '/esarchive/scratch/aorihuel/mask/change_mask/shapefiles_py/all_regions_1shapefile/IR_LAB_budgets'
areacello = iris.load_cube('/esarchive/exp/nemo/a5gj/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Ofx/areacello/gn/v20221216/areacello_Ofx_EC-Earth3_opera-control_r1i1p1f1_gn.nc')

rowName_List = ["aLAB_IR"]


def calculate_statistics(model_data, obs_data):
    """
    Calcula estadísticos entre datos del modelo y observaciones.
    
    Parameters:
    -----------
    model_data : numpy.ma.array
        Datos del modelo (masked array)
    obs_data : numpy.ma.array
        Datos de observación (masked array)
    
    Returns:
    --------
    dict con los estadísticos calculados
    """
    # Crear máscara común donde ambos tienen datos válidos
    model_data_LAB_IR = extract_shape(model_data, shapefile, method='contains', crop=False, decomposed=True, ids=rowName_List)
    obs_data_LAB_IR = extract_shape(obs_data, shapefile, method='contains', crop=False, decomposed=True, ids=rowName_List)

    
    # Verificar que hay datos válidos
    if len(model_data_LAB_IR.data) == 0 or len(obs_data_LAB_IR.data) == 0:
        return {
            'rmsd': np.nan,
            'rmsd_rel': np.nan,
            'bias': np.nan,
            'rel_bias': np.nan
        }
        #'ubrmsd': np.nan
    
    # Diferencias
    diff = model_data_LAB_IR.copy(model_data_LAB_IR.data - obs_data_LAB_IR.data)
    
    # Mean Bias (MB)
    add_supplementary_variables(diff, [areacello])
    bias = area_statistics(diff, 'mean')
    #bias = np.mean(diff)
    
    # Root Mean Square Deviation (RMSD)
    
    #rmsd = np.sqrt(np.mean(diff**2))
    sq_diff = diff.copy(diff.data**2)
    mean_sq_diff= area_statistics(sq_diff, 'mean')
    rmsd = mean_sq_diff.copy(np.sqrt(mean_sq_diff.data))
    
    # RMSD Relative (%)
    #obs_mean = np.mean(obs_data_LAB_IR)
    add_supplementary_variables(obs_data_LAB_IR, [areacello])
    obs_mean = area_statistics(obs_data_LAB_IR, 'mean')
    if obs_mean != 0:
        rmsd_rel = rmsd.copy((rmsd.data / obs_mean.data) * 100)
    else:
        rmsd_rel = np.nan
    
    # Relative Bias (%)
    if obs_mean != 0:
        rel_bias = bias.copy((bias.data / obs_mean.data) * 100)
    else:
        rel_bias = np.nan
    
    # Unbiased Root Mean Square Deviation (ubRMSD)
    #mean_sq_diff_minus_bias= area_statistics(diff.copy((diff.data - bias.data)**2), 'mean')
    #ubrmsd =mean_sq_diff_minus_bias.copy(np.sqrt(mean_sq_diff_minus_bias.data))
    
    return {
        'rmsd': rmsd.data,
        'rmsd_rel': rmsd_rel.data,
        'bias': bias.data,
        'rel_bias': rel_bias.data
    }
    # 'ubrmsd': ubrmsd.data
# Calcular estadísticos para cada experimento, período y mes
statistics = {}

for expid in diag_list:
    statistics[expid] = {}
    
    for period_key in years_groups.keys():
        statistics[expid][period_key] = {}
        
        # Seleccionar datos del período correcto
        if period_key == '2009-2013':
            model_mon_data = weak_mon_data[expid]
            obs_mon_data = weak_sat_mon_data
        else:  # '2014-2018'
            model_mon_data = strong_mon_data[expid]
            obs_mon_data = strong_sat_mon_data
        
        for mon_idx, mon in enumerate(extracted_months):
            month_name = mon_list[mon_idx]
            
            # Convertir datos del modelo de mol m-3 a mg C m-3 (multiplicar por 1000*12)
            # y observaciones dividir por 12 para tener las mismas unidades
            #mg C m-3
            model_data = model_mon_data[month_name].copy(model_mon_data[month_name].data * 1000*12)
            obs_data =  obs_mon_data[month_name].copy(obs_mon_data[month_name].data)

            #mmol C m-3
            #model_data = model_mon_data[month_name].copy(model_mon_data[month_name].data * 1000)
            #obs_data =  obs_mon_data[month_name].copy(obs_mon_data[month_name].data / 12)
            
            # Calcular estadísticos
            stats = calculate_statistics(model_data, obs_data)
            statistics[expid][period_key][month_name] = stats

# También calcular estadísticos para las observaciones vs observaciones (será todo ceros)
# Esto es para mantener la estructura consistente en el plotting

statistics['obs_SAT'] = {}
for period_key in years_groups.keys():
    statistics['obs_SAT'][period_key] = {}
    
    if period_key == '2009-2013':
        obs_mon_data = weak_sat_mon_data
    else:
        obs_mon_data = strong_sat_mon_data
    
    for mon_idx, mon in enumerate(extracted_months):
        month_name = mon_list[mon_idx]
        obs_mon_data_LAB_IR = extract_shape(obs_mon_data[month_name], shapefile, method='contains', crop=False, decomposed=True, ids=rowName_List)
        add_supplementary_variables(obs_mon_data_LAB_IR, [areacello])
        obs_mean = area_statistics(obs_mon_data_LAB_IR, 'mean')
        ## comentado pq no es weighted obs_std_dev = area_statistics(obs_mon_data_LAB_IR, 'std_dev')
        # Obtener grid weights
        # Calcular std_dev ponderada manualmente con una funcion
        #def weighted_std_dev(cube):
            # Obtener los pesos de área
         #   weights = cube.coord('cell_area').points
    
            # Datos
          #  data = cube.data
    
            # Media ponderada
           # weighted_mean = np.average(data, weights=weights)
    
            # Varianza ponderada
            #weighted_variance = np.average((data - weighted_mean)**2, weights=weights)
    
            # Desviación estándar ponderada
            #weighted_std = np.sqrt(weighted_variance)
    
            #return weighted_std

        #obs_std_dev_weighted = weighted_std_dev(obs_mon_data_LAB_IR)
        weights = areacello.data
    
        # Datos del cubo principal
        data = obs_mon_data_LAB_IR.data
    
        # Aplanar los arrays si son multidimensionales
        data_flat = data.flatten()
        weights_flat = weights.flatten()
    
    
        # Media ponderada
        weighted_mean = np.average(data_flat, weights=weights_flat)
    
        # Varianza ponderada
        weighted_variance = np.average((data_flat - weighted_mean)**2, weights=weights_flat)
    
        # Desviación estándar ponderada
        weighted_std = np.sqrt(weighted_variance)
        
        statistics['obs_SAT'][period_key][month_name] = {
            'Mean': obs_mean.data,
            'Std Dev': weighted_std
        }

print("Estadísticos calculados!")
print(f"Experimentos: {list(statistics.keys())}")
print(f"Ejemplo de estadísticos para {diag_list[0]}, 2009-2013, Mar:")
print(statistics[diag_list[0]]['2009-2013']['Mar'])

## incluimos los estadisticos en los mapas que estan en mg C m-3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
import iris.analysis
import iris.analysis.cartography
import cartopy.crs as ccrs
import cmocean
import matplotlib.patches as mpatches

mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['xtick.major.pad'] = 15
mpl.rcParams['ytick.major.pad'] = 15

# === CONFIG ===
years_groups = {
    '2009-2013': ['2009', '2010', '2011', '2012', '2013'],
    '2014-2018': ['2014', '2015', '2016', '2017', '2018']
}

# Configuración para diferencias (en lugar de ratio)
cmap_ratio = cmocean.cm.balance  # Colormap balanceado para diferencias (centrado en 0)
vmin_ratio = -50  # Límite inferior para diferencias (mgC m-3)
vmax_ratio = 50   # Límite superior para diferencias (mgC m-3)
levels_ratio = np.linspace(vmin_ratio, vmax_ratio, 50)  # Niveles lineales en lugar de logarítmicos
ticks_ratio = [-50, -30, -15, 0, 15, 30, 50]  # Ticks apropiados para diferencias
norm_ratio = mpl.colors.Normalize(vmin=vmin_ratio, vmax=vmax_ratio)  # Normalización lineal

nx, ny = 362, 292
target_proj = ccrs.PlateCarree()

def add_statistics_box(ax, stats, type_data, fontsize=7):
    """
    Añade un recuadro con los estadísticos en la esquina inferior derecha del mapa.
    
    Parameters:
    -----------
    ax : matplotlib axis
        El eje donde añadir el recuadro
    stats : dict
        Diccionario con los estadísticos calculados
    fontsize : int
        Tamaño de la fuente
    type_data : str
        Si se trata de observacion o modelo
    """
    if type_data == 'Sim':
        # Crear texto con los estadísticos
        stats_text = (
            f"RMSD: {stats['rmsd']:.2f}\n"
            f"RMSD%: {stats['rmsd_rel']:.1f}\n"
            f"MB: {stats['bias']:.2f}\n"
            f"RB%: {stats['rel_bias']:.1f}"
        )
        #f"ubRMSD: {stats['ubrmsd']:.2f}\n"
    else:
        # Crear texto con los estadísticos
        stats_text = (
            f"Mean: {stats['Mean']:.2f}\n"
            f"Std Dev: {stats['Std Dev']:.1f}"
        )

    # Añadir recuadro de texto
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.text(0.98, 0.02, stats_text, transform=ax.transAxes,
            fontsize=fontsize, verticalalignment='bottom', horizontalalignment='right',
            bbox=props, family='monospace')

def plot_combined_ratio_figure():
    # Crear figura con 6 filas y 6 columnas
    fig, axs = plt.subplots(6, 6, figsize=(20, 11),
                            subplot_kw={'projection': target_proj}, 
                            gridspec_kw={'hspace': 0.07, 'wspace': 0.1})
    
    # Mapeo de experimentos a labels
    exp_labels = {'a67o': 'Sim_REF', 'a683': 'Sim_SAT'}
    
    # Configuración para la colorbar de observaciones (valores absolutos)
    cmap_obs = cmocean.cm.algae  # Colormap para valores positivos
    vmin_obs = 0
    vmax_obs = 100  # Ajustar según el rango de tus datos (mmol m-2) (mg C m-2)
    levels_obs = np.linspace(vmin_obs, vmax_obs, 20)
    ticks_obs =  [0, 25, 50, 75, 100]
    #ticks_obs =  [0, 10, 20, 30, 40, 50]
    norm_obs = mpl.colors.Normalize(vmin=vmin_obs, vmax=vmax_obs)
    
    # Primeras dos filas: Observaciones de satélite
    for row_idx, (period_key, years) in enumerate(years_groups.items()):
        for col_idx, month in enumerate(extracted_months):
            ax = axs[row_idx, col_idx]
            
            # Datos de observaciones de satélite
            obs_data = all_obs_dict[period_key][mon_list[col_idx]]
            
            # Proyección de datos de satélite
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                obs_proj, obs_extent = iris.analysis.cartography.project(
                    obs_data, target_proj, nx=nx, ny=ny)
            
            # Get coordinates
            x = obs_proj.coord('longitude').points
            y = obs_proj.coord('latitude').points
            
            extent = (-62.0, -20.0, 50.0, 68.0)
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            
            # Plotear observaciones de satélite
            im_obs = ax.contourf(x, y, obs_proj.data, 
                                levels=levels_obs, norm=norm_obs,
                                cmap=cmap_obs, transform=target_proj, extend='max')
            
            gl = ax.gridlines(draw_labels=False, alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False
            ax.coastlines()
            ax.add_feature(cfea.LAND, zorder=1)
            
            # Añadir estadísticos (para obs vs obs serán todos 0)
            month_name = mon_list[col_idx]
            stats = statistics['obs_SAT'][period_key][month_name]
            add_statistics_box(ax, stats, type_data = 'Obs', fontsize=6.5)
            
            # Títulos de los meses en la primera fila
            if row_idx == 0:
                ax.set_title(f"{['Mar','Apr','May','Jun','Jul','Aug'][col_idx]}")
                if col_idx == 0:
                    ax.text(0.5, 1.5, "Obs", va='center', ha='right',
                                rotation=0, transform=ax.transAxes, fontsize=15, 
                                fontweight='bold')
            
            # Labels para las observaciones en la primera columna
            if col_idx == 0:
                ax.text(-0.35, 0.5, period_key, va='center', ha='right',
                        rotation=90, transform=ax.transAxes, fontsize=12)
                gl.left_labels = True
    
    # Iterar por cada experimento (ahora empezando desde la fila 2)
    for exp_idx, (expid, exp_label) in enumerate(exp_labels.items()):
        start_row = exp_idx * 2 + 2  # a67o: filas 2-3, a683: filas 4-5
        
        for col_idx, month in enumerate(extracted_months):
            for row_idx, (period_key, years) in enumerate(years_groups.items()):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    ratio_proj, ratio_extent = iris.analysis.cartography.project(
                        all_ratios_dict[period_key][expid][mon_list[col_idx]], 
                        target_proj, nx=nx, ny=ny)
                
                # Get coordinates          
                x = ratio_proj.coord('longitude').points
                y = ratio_proj.coord('latitude').points
                ax = axs[start_row + row_idx, col_idx]
                
                extent = (-62.0, -20.0, 50.0, 68.0)
                ax.set_extent(extent, crs=ccrs.PlateCarree())
                
                im = ax.contourf(x, y, ratio_proj.data, 
                                 levels=levels_ratio, norm=norm_ratio,
                                 cmap=cmap_ratio, transform=target_proj, extend='both')
                
                # Configurar batimetría según el experimento
                if expid == 'a67o':
                    bathy_data = bathy_a67o
                    color = 'grey'
                else:  # expid == 'a683'
                    bathy_data = bathy_a683
                    color = 'dimgrey'
                
                # Proyectar batimetría al mismo grid
                bathy_proj, _ = iris.analysis.cartography.project(bathy_data, target_proj, nx=nx, ny=ny)
                
                # Dibujar contornos de batimetría (1000m, 2000m, 2500m)
                ax.contour(x, y, bathy_proj.data, levels=[1000, 2000, 2500], 
                          colors=color, linewidths=1, linestyles='-', 
                          transform=target_proj)
                
                gl = ax.gridlines(draw_labels=False, alpha=0.5)
                gl.top_labels = False
                gl.right_labels = False
                ax.coastlines()
                ax.add_feature(cfea.LAND, zorder=1)
                
                # Añadir estadísticos para este experimento y período
                month_name = mon_list[col_idx]
                stats = statistics[expid][period_key][month_name]
                add_statistics_box(ax, stats, type_data = 'Sim', fontsize=6.5)
                
                # Labels de períodos y experimentos en la primera columna
                if col_idx == 0:
                    # Label del período
                    ax.text(-0.35, 0.5, period_key, va='center', ha='right',
                            rotation=90, transform=ax.transAxes, fontsize=12)
                    
                    # Label del experimento (centrado verticalmente entre las dos filas del experimento)
                    if row_idx == 0:  # Solo en la primera fila de cada experimento
                        ax.text(0.5, 1.25, exp_label, va='center', ha='right',
                                rotation=0, transform=ax.transAxes, fontsize=15, 
                                fontweight='bold')
                    
                    gl.left_labels = True
                
                # Labels en la fila inferior de cada experimento
                if start_row + row_idx == 5:
                    gl.bottom_labels = True
    
    # Colorbar común para diferencias (filas 2-5)
    cbar = fig.colorbar(im, ax=axs[2:], orientation='vertical', shrink=0.75, pad=0.02,
                        ticks=ticks_ratio)
    cbar.set_label(r' Sim - Obs $(mg C· m^{-3})$')
    cbar.ax.set_yticklabels(["-50", "-30", "-15", "0", "15", "30", "50"])  # Etiquetas para diferencias
   
    # Colorbar para observaciones (filas 0-1)
    cbar2 = fig.colorbar(im_obs, ax=axs[0:2], orientation='vertical', shrink=0.75, pad=0.02,
                        ticks=ticks_obs)
    cbar2.set_label(f'Obs $POC_{{{var}}}$ $(mg\\ C\cdot m^{{-3}})$') 
    cbar2.ax.set_yticklabels(["0", "25", "50", "75", "100"]) 
    
    # Guardar figura
    fig.savefig('/esarchive/scratch/aorihuel/my_py/BUDGETS_LABRADOR/POC_'+var+'_diff_combined_LB_IR_a683_estadisticos2_mgCm3.png', 
                bbox_inches='tight', dpi=300)
    fig.savefig('/esarchive/scratch/aorihuel/my_py/BUDGETS_LABRADOR/POC_'+var+'_diff_combined_LB_IR_a683_estadisticos2_mgCm3.svg', 
                bbox_inches='tight')
    plt.show()

# === CREAR FIGURA COMBINADA ===
plot_combined_ratio_figure()