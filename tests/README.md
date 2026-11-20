## CI tests for MYNN-EDMF
* Cases include:
    * clr, land point (43.40°N, 73.02°W) on Aug 23, 2024 
    * non-precip StCu, ocean point (40.84°N, 66.80°W) on Mar 07, 2024
    * dycoms, DYCOMS-II RF01 nocturnal stratocumulus (Stevens et al. 2005, Mon. Wea. Rev., 133, 1443-1462), run with bl_mynn_edmf = 1 and 2
* Input:

    * Input files for clr and marine_StCu are generated using WRF
    * input_dycoms.nc is built from input_marine_StCu.nc by make_input_dycoms.py: same grid and format, with the RF01 initial profile (theta_l, q_t, winds, surface pressure and fluxes) and the RF01 longwave heating held fixed as RTHRATEN

* To run the CI tests offline:
    ```
    export NFHOME=/path/to/your/netcdf/dir
    cd tests
    make
    ./driver
    ```
