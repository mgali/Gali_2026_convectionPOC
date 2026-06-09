import os
import iris
import pandas as pd
import numpy as np
from pathlib import Path
import os
import iris
from esmvalcore.preprocessor import climate_statistics
import numpy as np

"""
def create_yearly_climatologies():
    
    ###Create individual annual climatologies for each year (1998-2018) 
    ###for oceanographic particle budget data
   
    
    # Define parameters
    years = list(range(1998, 2020))  # 1998 to 2018
    
    # Particle types and their corresponding variables
    particle_configs = {
        'pkt': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/pkt/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc']
        },
        'small_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/small_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/large_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    # Depth layers configuration
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    # Base output directory
    base_output_dir = '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies'
    
    # Process each particle type
    for particle_type, config in particle_configs.items():
        print(f"\nProcessing {particle_type}...")
        
        # Process each layer
        for i, layer in enumerate(layer_list):
            zmin = ZMIN_LIST[i]
            zmax = ZMAX_LIST[i]
            
            print(f"  Processing layer: {layer} ({zmin}-{zmax}m)")
            
            # Process each variable
            for var in config['variables']:
                print(f"    Processing variable: {var}")
                
                # Create output directory for this variable/layer combination
                output_dir = os.path.join(base_output_dir, particle_type, layer, var)
                os.makedirs(output_dir, exist_ok=True)
                
                # Process each year individually
                for year in years:
                    file_path = config['path_template'].format(
                        year=year,
                        layer=layer,
                        var=var,
                        zmin=zmin,
                        zmax=zmax
                    )
                    
                    if os.path.exists(file_path):
                        try:
                            print(f"      Processing year {year}...")
                            
                            # Load the yearly file
                            cube = iris.load_cube(file_path)
                            
                            # Calculate annual climatology (if the file contains monthly data)
                            # If it's already annual data, this will just return the same cube
                            cube_annual = climate_statistics(cube, operator='mean', period='full')
                            
                            # Create output filename
                            output_filename = f'{var}_{particle_type}_{layer}_{zmin}_{zmax}_climatology_{year}.nc'
                            full_output_path = os.path.join(output_dir, output_filename)
                            
                            # Save the yearly climatology
                            iris.save(cube_annual, full_output_path)
                            
                            print(f"        Saved: {output_filename}")
                            
                        except Exception as e:
                            print(f"        Error processing {year}: {str(e)}")
                            continue
                    else:
                        print(f"        Missing file for {year}")
    
    print("\nYearly climatology computation completed!")
"""

def create_yearly_climatologies_in_original_paths():
    """
    Create individual annual climatologies for each year (1998-2018) 
    and save them in the same directory structure as the original files
    """
    
    # Define parameters
    years = list(range(1998, 2020))  # 1998 to 2018
    
    # Particle types and their corresponding variables
    particle_configs = {
        'pkt': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/pkt/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc']
        },
        'small_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/small_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/large_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    # Depth layers configuration
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    # Process each particle type
    for particle_type, config in particle_configs.items():
        print(f"\nProcessing {particle_type}...")
        
        # Process each layer
        for i, layer in enumerate(layer_list):
            zmin = ZMIN_LIST[i]
            zmax = ZMAX_LIST[i]
            
            print(f"  Processing layer: {layer} ({zmin}-{zmax}m)")
            
            # Process each variable
            for var in config['variables']:
                print(f"    Processing variable: {var}")
                
                # Process each year individually
                for year in years:
                    file_path = config['path_template'].format(
                        year=year,
                        layer=layer,
                        var=var,
                        zmin=zmin,
                        zmax=zmax
                    )
                    
                    if os.path.exists(file_path):
                        try:
                            print(f"      Processing year {year}...")
                            
                            # Load the yearly file
                            cube = iris.load_cube(file_path)
                            
                            # Calculate annual climatology (if the file contains monthly data)
                            # If it's already annual data, this will just return the same cube
                            cube_annual = climate_statistics(cube, operator='mean', period='full')
                            
                            # Create output filename in the same directory as the original file
                            original_dir = os.path.dirname(file_path)
                            output_filename = f'{var}_{particle_type}_{layer}_{zmin}_{zmax}_climatology_{year}.nc'
                            full_output_path = os.path.join(original_dir, output_filename)
                            
                            # Save the yearly climatology
                            iris.save(cube_annual, full_output_path)
                            
                            print(f"        Saved: {full_output_path}")
                            
                        except Exception as e:
                            print(f"        Error processing {year}: {str(e)}")
                            continue
                    else:
                        print(f"        Missing file for {year}")
    
    print("\nYearly climatology computation completed!")

def check_file_availability():
    """
    Check which files are available before processing
    """
    years = list(range(1998, 2020))
    
    particle_configs = {
        'pkt': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/pkt/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc' ]
        },
        'small_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/small_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/large_POC_budgets/diag_a683/yr_files/{year}/sp_average/aLAB_IR/bathy_depth_int/{layer}/{var}_a683_sp_average_IR_LAB_{zmin}_{zmax}_{layer}.nc',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    print("File availability check:")
    print("=" * 50)
    
    for particle_type, config in particle_configs.items():
        print(f"\n{particle_type.upper()}:")
        for i, layer in enumerate(layer_list):
            zmin = ZMIN_LIST[i]
            zmax = ZMAX_LIST[i]
            print(f"  {layer} ({zmin}-{zmax}m):")
            
            for var in config['variables']:
                available_years = []
                for year in years:
                    file_path = config['path_template'].format(
                        year=year, layer=layer, var=var, zmin=zmin, zmax=zmax
                    )
                    if os.path.exists(file_path):
                        available_years.append(str(year))
                
                if available_years:
                    print(f"    {var}: {', '.join(available_years)} ({len(available_years)}/10 years)")
                else:
                    print(f"    {var}: No files found")

if __name__ == "__main__":
    # First, check file availability
    print("Checking file availability...")
    check_file_availability()
    
    # Ask user which approach they prefer
    #print("\nChoose climatology computation approach:")
    #print("1. Save yearly climatologies in organized structure (yearly_climatologies/)")
    #print("2. Save yearly climatologies in original file directories")
    
    #choice = input("Enter choice (1 or 2): ")
    
    #if choice == "1":
    #    response = input("\nDo you want to proceed with organized yearly climatology computation? (y/n): ")
    #    if response.lower() in ['y', 'yes']:
    #        create_yearly_climatologies()
    #    else:
    #        print("Climatology computation cancelled.")
    #elif choice == "2":
    #    response = input("\nDo you want to proceed with in-place yearly climatology computation? (y/n): ")
    #    if response.lower() in ['y', 'yes']:
 #           create_yearly_climatologies_in_original_paths()
#        else:
#            print("Climatology computation cancelled.")
#    else:
 #       print("Invalid choice. Exiting.")
    create_yearly_climatologies_in_original_paths()
    
def extract_cube_data(cube):
    """
    Extract data from an iris cube and return as a single value or array summary
    """
    try:
        # Get the data array
        data = cube.data
        
        # Handle masked arrays
        if hasattr(data, 'mask'):
            data = np.ma.filled(data, np.nan)
        
        # If it's a scalar, return it directly
        if data.ndim == 0:
            return float(data)
        
        # For multi-dimensional data, you might want to compute statistics
        # Here we'll return the mean, but you can modify this based on your needs
        if np.isfinite(data).any():
            return float(np.nanmean(data))
        else:
            return np.nan
            
    except Exception as e:
        print(f"    Warning: Could not extract data from cube: {str(e)}")
        return np.nan
"""
def climatology_to_single_excel():
    
    ### Convert all climatology NetCDF files to a single Excel file
    ### One sheet per layer, all particle types and variables as columns, years as rows
    
    
    # Define parameters
    years = list(range(1998, 2020))  # 1998 to 2020
    
    # Particle types and their corresponding variables
    particle_configs = {
        'pkt': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/pkt/{layer}/{var}/{var}_pkt_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc']
        },
        'small_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/small_POC_budgets/{layer}/{var}/{var}_small_POC_budgets_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/large_POC_budgets/{layer}/{var}/{var}_large_POC_budgets_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    # Depth layers configuration
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    # Output directory
    output_dir = '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/excel_exports'
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a dictionary to store data for all layers
    excel_data = {}
    
    print("Processing all particle types and creating single Excel file...")
    
    # Process each layer
    for i, layer in enumerate(layer_list):
        zmin = ZMIN_LIST[i]
        zmax = ZMAX_LIST[i]
        
        print(f"\nProcessing layer: {layer} ({zmin}-{zmax}m)")
        
        # Initialize data structure for this layer
        layer_data = {}
        
        # Initialize with years as index
        layer_data['Year'] = years
        
        # Process each particle type
        for particle_type, config in particle_configs.items():
            print(f"  Processing particle type: {particle_type}")
            
            # Process each variable for this particle type
            for var in config['variables']:
                print(f"    Processing variable: {var}")
                
                # Create column name with particle type prefix
                column_name = f"{particle_type}_{var}"
                
                # Collect data for all years for this variable
                var_data = []
                
                for year in years:
                    file_path = config['path_template'].format(
                        year=year,
                        layer=layer,
                        var=var,
                        zmin=zmin,
                        zmax=zmax
                    )
                    
                    if os.path.exists(file_path):
                        try:
                            # Load the climatology file
                            cube = iris.load_cube(file_path)
                            
                            # Extract data value
                            value = extract_cube_data(cube)
                            var_data.append(value)
                            
                        except Exception as e:
                            print(f"      Error loading {year}: {str(e)}")
                            var_data.append(np.nan)
                    else:
                        print(f"      Missing file for {year}")
                        var_data.append(np.nan)
                
                # Add variable data to layer data with particle type prefix
                layer_data[column_name] = var_data
        
        # Convert to DataFrame
        df_layer = pd.DataFrame(layer_data)
        df_layer.set_index('Year', inplace=True)
        
        # Store in excel_data
        excel_data[f"{layer}_{zmin}-{zmax}m"] = df_layer
    
    # Save to single Excel file with multiple sheets
    excel_filename = 'all_climatologies_1998-2018.xlsx'
    excel_path = os.path.join(output_dir, excel_filename)
    
    print(f"\nSaving to Excel: {excel_filename}")
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        for sheet_name, df in excel_data.items():
            df.to_excel(writer, sheet_name=sheet_name)
            print(f"  Sheet created: {sheet_name} ({df.shape[1]} variables)")
    
    print(f"\nCompleted: {excel_path}")
    print("Excel export completed!")
"""

def climatology_to_single_excel_from_original_paths():
    """
    Convert all climatology NetCDF files to a single Excel file
    This version reads from the original directory structure
    """
    
    # Define parameters
    years = list(range(1998, 2020))  # 1998 to 2018
    
    # Particle types and their corresponding variables
    particle_configs = {
        'pkt': {
            'base_path': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/pkt/diag_a5xl/yr_files',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc']
        },
        'small_POC_budgets': {
            'base_path': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/small_POC_budgets/diag_a5xl/yr_files',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'base_path': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/large_POC_budgets/diag_a5xl/yr_files',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    # Depth layers configuration
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    # Output directory
    output_dir = '/esarchive/scratch/aorihuel/my_py/BUDGETS_LABRADOR/'
    os.makedirs(output_dir, exist_ok=True)
    
    # Create a dictionary to store data for all layers
    excel_data = {}
    
    print("Processing all particle types and creating single Excel file...")
    
    # Process each layer
    for i, layer in enumerate(layer_list):
        zmin = ZMIN_LIST[i]
        zmax = ZMAX_LIST[i]
        
        print(f"\nProcessing layer: {layer} ({zmin}-{zmax}m)")
        
        # Initialize data structure for this layer
        layer_data = {}
        
        # Initialize with years as index
        layer_data['Year'] = years
        
        # Process each particle type
        for particle_type, config in particle_configs.items():
            print(f"  Processing particle type: {particle_type}")
            
            # Process each variable for this particle type
            for var in config['variables']:
                print(f"    Processing variable: {var}")
                
                # Create column name with particle type prefix
                #column_name = f"{particle_type}_{var}"
                column_name = f"{var}"
                
                # Collect data for all years for this variable
                var_data = []
                
                for year in years:
                    # Construct file path for climatology file
                    file_path = os.path.join(
                        config['base_path'],
                        str(year),
                        'sp_average/aLAB_IR/bathy_depth_int',
                        layer,
                        f'{var}_{particle_type}_{layer}_{zmin}_{zmax}_climatology_{year}.nc'
                    )
                    
                    if os.path.exists(file_path):
                        try:
                            # Load the climatology file
                            cube = iris.load_cube(file_path)
                            
                            # Extract data value
                            value = extract_cube_data(cube)
                            var_data.append(value)
                            
                        except Exception as e:
                            print(f"      Error loading {year}: {str(e)}")
                            var_data.append(np.nan)
                    else:
                        print(f"      Missing file for {year}")
                        var_data.append(np.nan)
                
                # Add variable data to layer data with particle type prefix
                layer_data[column_name] = var_data
        
        # Convert to DataFrame
        df_layer = pd.DataFrame(layer_data)
        df_layer.set_index('Year', inplace=True)
        
        # Store in excel_data
        excel_data[f"{layer}_{zmin}-{zmax}m"] = df_layer
    
    # Save to single Excel file with multiple sheets
    excel_filename = 'all_climatologies_1998-2020_a5xl.xlsx'
    excel_path = os.path.join(output_dir, excel_filename)
    
    print(f"\nSaving to Excel: {excel_filename}")
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        for sheet_name, df in excel_data.items():
            df.to_excel(writer, sheet_name=sheet_name)
            print(f"  Sheet created: {sheet_name} ({df.shape[1]} variables)")
    
    print(f"\nCompleted: {excel_path}")
    print("Excel export completed!")

def check_climatology_files():
    """
    Check which climatology files are available
    """
    years = list(range(1998, 2020))

    # Check organized structure first
    particle_configs = {
        'pkt': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/pkt/{layer}/{var}/{var}_pkt_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['tottrdpkt', 'xadpkt', 'yadpkt', 'zadpkt', 'ldfpkt', 'zdfpkt', 'zingest2doc']
        },
        'small_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/small_POC_budgets/{layer}/{var}/{var}_small_POC_budgets_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['remsdetoc', 'tottrdsdetoc', 'trdexpsdetoc', 'xadsdetoc', 'yadsdetoc', 'zadsdetoc', 'ldfsdetoc', 'zdfsdetoc']
        },
        'large_POC_budgets': {
            'path_template': '/esarchive/scratch/aorihuel/BUDGETS_LB_IR/yearly_climatologies/large_POC_budgets/{layer}/{var}/{var}_large_POC_budgets_{layer}_{zmin}_{zmax}_climatology_{year}.nc',
            'variables': ['remldetoc', 'tottrdldetoc', 'trdexpldetoc', 'xadldetoc', 'yadldetoc', 'zadldetoc', 'ldfldetoc', 'zdfldetoc']
        }
    }
    
    ZMIN_LIST = [0, 100, 500, 100, 1000]
    ZMAX_LIST = [100, 500, 1000, 1000, 2000]
    layer_list = ['epipelagic', 'upper_meso', 'lower_meso', 'meso', 'bathy']
    
    print("Climatology files availability check:")
    print("=" * 50)
    """
    organized_structure_exists = False
    
    for particle_type, config in particle_configs.items():
        print(f"\n{particle_type.upper()}:")
        for i, layer in enumerate(layer_list):
            zmin = ZMIN_LIST[i]
            zmax = ZMAX_LIST[i]
            print(f"  {layer} ({zmin}-{zmax}m):")
            
            for var in config['variables']:
                available_years = []
                for year in years:
                    file_path = config['path_template'].format(
                        year=year, layer=layer, var=var, zmin=zmin, zmax=zmax
                    )
                    if os.path.exists(file_path):
                        available_years.append(str(year))
                        organized_structure_exists = True
                
                if available_years:
                    print(f"    {var}: {', '.join(available_years)} ({len(available_years)}/10 years)")
                else:
                    print(f"    {var}: No files found")
    
    return organized_structure_exists
    """

if __name__ == "__main__":
    # Check which climatology files are available
    print("Checking climatology files availability...")
    organized_exists = check_climatology_files()
    """
    if organized_exists:
        response = input("\nFound climatology files in organized structure. Proceed with single Excel export? (y/n): ")
        if response.lower() in ['y', 'yes']:
            climatology_to_single_excel()
        else:
            print("Excel export cancelled.")
    
    else:
    
        print("\nNo climatology files found in organized structure.")
        """
    response = input("Try to export from original directory structure? (y/n): ")
    if response.lower() in ['y', 'yes']:
        climatology_to_single_excel_from_original_paths()
    else:
        print("Excel export cancelled.")