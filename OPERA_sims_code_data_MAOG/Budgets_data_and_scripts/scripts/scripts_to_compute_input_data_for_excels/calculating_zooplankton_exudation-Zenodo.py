import iris
import numpy as np
import os
from esmvalcore.preprocessor import extract_region, extract_levels, area_statistics, volume_statistics
from pathlib import Path
# Lista de variables necesarias
var_list = ['smssdetoc', 'smsldetoc', 'smspkt', 'smsdoc', 'pp', 'remoc', 'trdexpsdetoc', 'trdexpldetoc']

# Lista de experimentos y años
expid_list = ['a67o','a683', 'a5xl'] #, 'a5xl', 'a683', '']
#years = list(range(2009, 2019))  # 2009-2018
years = list(range(1998, 2020))
# Ruta base
base_path = "/esarchive/exp/nemo/{expid}/original_files/cmorfiles/OPERA/EC-Earth-Consortium/EC-Earth3/opera-control/r1i1p1f1/Omon/{var}/gn/v*/{var}_Omon_EC-Earth3_opera-control_r1i1p1f1_gn_{year}01-{year}12.nc"

# Directorio de salida
out_base = "/esarchive/scratch/aorihuel/BUDGETS_LB_IR/non_clim/{expid}/"
def load_variable(expid, year, var):
    path = base_path.format(expid=expid, var=var, year=year)
    try:
        cube = iris.load_cube(path)
        print(f"✓ Cargado {var} para {expid} en {year}")
        return cube
    except Exception as e:
        print(f"⚠️ Error cargando {var} para {expid} en {year}: {e}")
        return None
def compute_zingest(cubes):
    try:
        sms_total = (cubes['smssdetoc'].data + cubes['smsldetoc'].data +
                     cubes['smspkt'].data + cubes['smsdoc'].data)
        total_sinks = cubes['pp'].data - sms_total
        residual_sink = total_sinks - cubes['remoc'].data + cubes['trdexpsdetoc'].data + cubes['trdexpldetoc'].data
        
        zingest2dic = cubes['pp'].copy(data=residual_sink)
        zingest2dic.var_name = 'zingest2dic'
        zingest2dic.long_name = 'Residual sink to DIC'
        zingest2dic.units = cubes['pp'].units
        
        zingest2doc_data = residual_sink * 0.4 / 0.6
        zingest2doc = cubes['pp'].copy(data=zingest2doc_data)
        zingest2doc.var_name = 'zingest2doc'
        zingest2doc.long_name = 'Residual sink to DOC'
        zingest2doc.units = cubes['pp'].units
        
        return zingest2dic, zingest2doc
    except Exception as e:
        print(f"⚠️ Error en cálculo de zingest: {e}")
        return None, None
def save_cube_pathlib(cube, out_dir, filename):
    filo = Path(out_dir) / filename
    filo.parent.mkdir(parents=True, exist_ok=True)
    iris.save(cube, filo)
    print(f"💾 Guardado: {filo}")

for expid in expid_list:
    for year in years:
        cubes = {}
        for var in var_list:
            cube = load_variable(expid, year, var)
            if cube is None:
                break
            cubes[var] = cube
        
        if len(cubes) < len(var_list):
            print(f"⏭️ Saltando {expid} {year} por datos incompletos.")
            continue

        zingest2dic, zingest2doc = compute_zingest(cubes)
        if zingest2dic is None:
            continue
        
        out_dir = f"/esarchive/scratch/aorihuel/BUDGETS_LB_IR/non_clim/{expid}/"
        save_cube_pathlib(zingest2dic, out_dir+'zingest2dic/', f'zingest2dic_Omon_EC-Earth3_opera-control_r1i1p1f1_gn_{year}01-{year}12.nc')
        save_cube_pathlib(zingest2doc, out_dir+'zingest2doc/', f'zingest2doc_Omon_EC-Earth3_opera-control_r1i1p1f1_gn_{year}01-{year}12.nc')
