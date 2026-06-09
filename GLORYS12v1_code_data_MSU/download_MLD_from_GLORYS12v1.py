# Run this script with copernicusmarine environment
i_year = 1993
f_year = 2019
year_selection = np.arange(i_year,f_year+1) 

import copernicusmarine

for yy in year_selection:
    print(yy)

    copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
        variables=["mlotst"],
        minimum_longitude=-65,
        maximum_longitude=-15,
        minimum_latitude=40,
        maximum_latitude=70,
        start_datetime=f"{yy}-01-01T00:00:00",
        end_datetime=f"{yy}-12-31T00:00:00",
        #minimum_depth=0.49402499198913574,
        #maximum_depth=0.49402499198913574,
        output_filename = f"./daily_data/mlotst_{yy}.nc",
        output_directory = "data",
        force_download=True,  # Avoid asking is you wanna download
    )
