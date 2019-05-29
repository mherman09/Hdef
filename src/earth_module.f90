module earth

use trig, only: pi

! Radius of the Earth
double precision, parameter, public :: radius_earth_m = 6371.0d3
double precision, parameter, public :: radius_earth_km = 6371.0d0

! Kilometers per degree of great circle arc
double precision, parameter, public :: deg2km = 2.0d0*pi*radius_earth_km/360.0d0

! Geodetic Reference System 1980
double precision, parameter, public :: eccentricity_grs80 = sqrt(0.00669438002290d0)
double precision, parameter, public :: semimajor_grs80 = 6378137.0d0
double precision, parameter, public :: semiminor_grs80 = 6356752.3141d0
double precision, parameter, public :: radius_mean_grs80 = 6371008.7714d0
double precision, parameter, public :: radius_surface_grs80 = 6371007.181d0
double precision, parameter, public :: radius_volume_grs80 = 6371000.79d0
double precision, parameter, public :: gravity_equator_grs80 = 9.7803267715d0
double precision, parameter, public :: gravity_pole_grs80 = 9.8321863685d0

! World Geodetic System 1984
double precision, parameter, public :: semimajor_wgs84 = 6378137.0d0
double precision, parameter, public :: semiminor_wgs84 = 6356752.314245d0

public :: sphere_geo2xyz
public :: sphere_xyz2geo
public :: ellipsoid_geo2xyz
public :: ellipsoid_xyz2geo

! Milli-arc-seconds/year to degrees/million years
double precision, parameter, public :: masa2degma = 1.0d6/(1.0d3*60.0d0*60.0d0)

!----
! Plate motion models
!----
! ITRF 2005
! Altamimi, Z., Collileux, X., Legrand, J., Garayt, B., Boucher, C. (2007). ITRF2005: A new release
!     of the International Terrestrial Reference Frame based on time series of station positions
!     and Earth Orientation Parameters. Journal of Geophysical Research 112, B09401.
! Rotations relative to ITRF2005 reference frame (identical to NNR-NUVEL1A)
double precision, parameter, public :: itrf05_AM_lon = -102.789d0 ! Amur
double precision, parameter, public :: itrf05_AM_lat =   52.263d0
double precision, parameter, public :: itrf05_AM_vel =    0.269d0
double precision, parameter, public :: itrf05_AN_lon = -125.315d0 ! Antarctica
double precision, parameter, public :: itrf05_AN_lat =   59.813d0
double precision, parameter, public :: itrf05_AN_vel =    0.223d0
double precision, parameter, public :: itrf05_AR_lon =    5.061d0 ! Arabia
double precision, parameter, public :: itrf05_AR_lat =   49.642d0
double precision, parameter, public :: itrf05_AR_vel =    0.579d0
double precision, parameter, public :: itrf05_AU_lon =   37.367d0 ! Australia
double precision, parameter, public :: itrf05_AU_lat =   32.407d0
double precision, parameter, public :: itrf05_AU_vel =    0.628d0
double precision, parameter, public :: itrf05_CA_lon = -104.279d0 ! Caribbean
double precision, parameter, public :: itrf05_CA_lat =   39.318d0
double precision, parameter, public :: itrf05_CA_vel =    0.241d0
double precision, parameter, public :: itrf05_EU_lon =  -95.979d0 ! Eurasia
double precision, parameter, public :: itrf05_EU_lat =   56.330d0
double precision, parameter, public :: itrf05_EU_vel =    0.261d0
double precision, parameter, public :: itrf05_IN_lon =   21.841d0 ! India
double precision, parameter, public :: itrf05_IN_lat =   49.823d0
double precision, parameter, public :: itrf05_IN_vel =    0.614d0
double precision, parameter, public :: itrf05_NZ_lon = -101.441d0 ! Nazca
double precision, parameter, public :: itrf05_NZ_lat =   45.101d0
double precision, parameter, public :: itrf05_NZ_vel =    0.642d0
double precision, parameter, public :: itrf05_NA_lon =  -87.385d0 ! North America
double precision, parameter, public :: itrf05_NA_lat =   -4.291d0
double precision, parameter, public :: itrf05_NA_vel =    0.192d0
double precision, parameter, public :: itrf05_NU_lon =  -82.501d0 ! Nubia
double precision, parameter, public :: itrf05_NU_lat =   49.955d0
double precision, parameter, public :: itrf05_NU_vel =    0.269d0
double precision, parameter, public :: itrf05_OK_lon = -132.910d0 ! Okhotsk
double precision, parameter, public :: itrf05_OK_lat =  -32.041d0
double precision, parameter, public :: itrf05_OK_vel =    0.083d0
double precision, parameter, public :: itrf05_PA_lon =  112.873d0 ! Pacific
double precision, parameter, public :: itrf05_PA_lat =  -62.569d0
double precision, parameter, public :: itrf05_PA_vel =    0.682d0
double precision, parameter, public :: itrf05_SA_lon = -129.631d0 ! South America
double precision, parameter, public :: itrf05_SA_lat =  -16.800d0
double precision, parameter, public :: itrf05_SA_vel =    0.121d0
double precision, parameter, public :: itrf05_SM_lon =  -89.542d0 ! Somalia
double precision, parameter, public :: itrf05_SM_lat =   53.661d0
double precision, parameter, public :: itrf05_SM_vel =    0.309d0
double precision, parameter, public :: itrf05_YZ_lon = -109.737d0 ! Yangtze
double precision, parameter, public :: itrf05_YZ_lat =   59.425d0
double precision, parameter, public :: itrf05_YZ_vel =    0.310d0
double precision, parameter, public :: itrf05_ITRF05_lon = 0d0 ! ITRF 2005
double precision, parameter, public :: itrf05_ITRF05_lat = 0d0
double precision, parameter, public :: itrf05_ITRF05_vel = 0d0

! ITRF 2008
! Altamimi, Z., Metivier, L., Collilieux, X. (2012). ITRF2008 plate motion model. Journal of
!     Geophysical Research 117, B07402.
! Rotations relative to ITRF 2008 reference frame
double precision, parameter, public :: itrf08_AM_lon = -113.261d0 ! Amur
double precision, parameter, public :: itrf08_AM_lat =   62.265d0
double precision, parameter, public :: itrf08_AM_vel =    0.287d0
double precision, parameter, public :: itrf08_AN_lon = -129.261d0 ! Antarctica
double precision, parameter, public :: itrf08_AN_lat =   58.545d0
double precision, parameter, public :: itrf08_AN_vel =    0.209d0
double precision, parameter, public :: itrf08_AR_lon =   -2.572d0 ! Arabia
double precision, parameter, public :: itrf08_AR_lat =   50.984d0
double precision, parameter, public :: itrf08_AR_vel =    0.531d0
double precision, parameter, public :: itrf08_AU_lon =   37.928d0 ! Australia
double precision, parameter, public :: itrf08_AU_lat =   32.783d0
double precision, parameter, public :: itrf08_AU_vel =    0.630d0
double precision, parameter, public :: itrf08_CA_lon =  -87.421d0 ! Caribbean
double precision, parameter, public :: itrf08_CA_lat =   31.370d0
double precision, parameter, public :: itrf08_CA_vel =    0.354d0
double precision, parameter, public :: itrf08_EU_lon =  -98.835d0 ! Eurasia
double precision, parameter, public :: itrf08_EU_lat =   54.225d0
double precision, parameter, public :: itrf08_EU_vel =    0.257d0
double precision, parameter, public :: itrf08_IN_lon =   13.817d0 ! India
double precision, parameter, public :: itrf08_IN_lat =   50.517d0
double precision, parameter, public :: itrf08_IN_vel =    0.554d0
double precision, parameter, public :: itrf08_NZ_lon = -102.011d0 ! Nazca
double precision, parameter, public :: itrf08_NZ_lat =   45.701d0
double precision, parameter, public :: itrf08_NZ_vel =    0.631d0
double precision, parameter, public :: itrf08_NA_lon =  -86.974d0 ! North America
double precision, parameter, public :: itrf08_NA_lat =   -8.578d0
double precision, parameter, public :: itrf08_NA_vel =    0.186d0
double precision, parameter, public :: itrf08_NU_lon =  -80.973d0 ! Nubia
double precision, parameter, public :: itrf08_NU_lat =   50.054d0
double precision, parameter, public :: itrf08_NU_vel =    0.262d0
double precision, parameter, public :: itrf08_PA_lon =  111.639d0 ! Pacific
double precision, parameter, public :: itrf08_PA_lat =  -62.771d0
double precision, parameter, public :: itrf08_PA_vel =    0.677d0
double precision, parameter, public :: itrf08_SA_lon = -128.002d0 ! South America
double precision, parameter, public :: itrf08_SA_lat =  -21.315d0
double precision, parameter, public :: itrf08_SA_vel =    0.118d0
double precision, parameter, public :: itrf08_SM_lon =  -96.129d0 ! Somalia
double precision, parameter, public :: itrf08_SM_lat =   50.127d0
double precision, parameter, public :: itrf08_SM_vel =    0.325d0
double precision, parameter, public :: itrf08_SU_lon =  -87.309d0 ! Sunda
double precision, parameter, public :: itrf08_SU_lat =   44.243d0
double precision, parameter, public :: itrf08_SU_vel =    0.388d0
double precision, parameter, public :: itrf08_ITRF05_lon = 0d0 ! ITRF 2005
double precision, parameter, public :: itrf08_ITRF05_lat = 0d0
double precision, parameter, public :: itrf08_ITRF05_vel = 0d0
double precision, parameter, public :: itrf08_ITRF08_lon = 0d0 ! ITRF 2008
double precision, parameter, public :: itrf08_ITRF08_lat = 0d0
double precision, parameter, public :: itrf08_ITRF08_vel = 0d0

! MORVEL
! Demets, C., Gordon, R.G., Argus, D.F. (2010). Geologically current plate motions.
!     Geophysical Journal International 181, 1-80.
! Rotations relative to the Pacific (PA) plate
double precision, parameter, public :: morvel_PA_lon =    0.0d0 ! Pacific
double precision, parameter, public :: morvel_PA_lat =    0.0d0
double precision, parameter, public :: morvel_PA_vel =  0.000d0
double precision, parameter, public :: morvel_AM_lon =  -82.7d0 ! Amur
double precision, parameter, public :: morvel_AM_lat =   65.9d0
double precision, parameter, public :: morvel_AM_vel =  0.929d0
double precision, parameter, public :: morvel_AN_lon =  -78.5d0 ! Antarctica
double precision, parameter, public :: morvel_AN_lat =   65.9d0
double precision, parameter, public :: morvel_AN_vel =  0.887d0
double precision, parameter, public :: morvel_AR_lon =  -33.2d0 ! Arabia
double precision, parameter, public :: morvel_AR_lat =   60.0d0
double precision, parameter, public :: morvel_AR_vel =  1.159d0
double precision, parameter, public :: morvel_AU_lon =    6.3d0 ! Australia
double precision, parameter, public :: morvel_AU_lat =   60.1d0
double precision, parameter, public :: morvel_AU_vel =  1.079d0
double precision, parameter, public :: morvel_CA_lon =  -77.5d0 ! Caribbean
double precision, parameter, public :: morvel_CA_lat =   55.8d0
double precision, parameter, public :: morvel_CA_vel =  0.905d0
double precision, parameter, public :: morvel_CO_lon = -112.8d0 ! Cocos
double precision, parameter, public :: morvel_CO_lat =   42.2d0
double precision, parameter, public :: morvel_CO_vel =  1.676d0
double precision, parameter, public :: morvel_CP_lon =  -10.1d0 ! Capricorn (Indian Ocean)
double precision, parameter, public :: morvel_CP_lat =   62.3d0
double precision, parameter, public :: morvel_CP_vel =  1.139d0
double precision, parameter, public :: morvel_EU_lon =  -78.9d0 ! Eurasia
double precision, parameter, public :: morvel_EU_lat =   61.3d0
double precision, parameter, public :: morvel_EU_vel =  0.856d0
double precision, parameter, public :: morvel_IN_lon =  -31.2d0 ! India
double precision, parameter, public :: morvel_IN_lat =   61.4d0
double precision, parameter, public :: morvel_IN_vel =  1.141d0
double precision, parameter, public :: morvel_JF_lon =   37.8d0 ! Juan de Fuca
double precision, parameter, public :: morvel_JF_lat =   -0.6d0
double precision, parameter, public :: morvel_JF_vel =  0.625d0
double precision, parameter, public :: morvel_LW_lon =  -66.9d0 ! Lwandle (South Madagascar)
double precision, parameter, public :: morvel_LW_lat =   60.0d0
double precision, parameter, public :: morvel_LW_vel =  0.932d0
double precision, parameter, public :: morvel_MQ_lon =   -8.0d0 ! Macquarie (South of New Zealand)
double precision, parameter, public :: morvel_MQ_lat =   59.2d0
double precision, parameter, public :: morvel_MQ_vel =  1.686d0
double precision, parameter, public :: morvel_NA_lon =  -71.7d0 ! North America
double precision, parameter, public :: morvel_NA_lat =   48.9d0
double precision, parameter, public :: morvel_NA_vel =  0.750d0
double precision, parameter, public :: morvel_NU_lon =  -66.6d0 ! Nubia
double precision, parameter, public :: morvel_NU_lat =   58.7d0
double precision, parameter, public :: morvel_NU_vel =  0.935d0
double precision, parameter, public :: morvel_NZ_lon =  -87.8d0 ! Nazca
double precision, parameter, public :: morvel_NZ_lat =   55.9d0
double precision, parameter, public :: morvel_NZ_vel =  1.311d0
double precision, parameter, public :: morvel_PS_lon =  -41.9d0 ! Philippine Sea
double precision, parameter, public :: morvel_PS_lat =   -4.6d0
double precision, parameter, public :: morvel_PS_vel =  0.890d0
double precision, parameter, public :: morvel_RI_lon = -104.8d0 ! Rivera (North of Cocos Plate)
double precision, parameter, public :: morvel_RI_lat =   25.7d0
double precision, parameter, public :: morvel_RI_vel =  4.966d0
double precision, parameter, public :: morvel_SA_lon =  -77.0d0 ! South America
double precision, parameter, public :: morvel_SA_lat =   56.0d0
double precision, parameter, public :: morvel_SA_vel =  0.653d0
double precision, parameter, public :: morvel_SC_lon =  -78.0d0 ! Scotia
double precision, parameter, public :: morvel_SC_lat =   57.8d0
double precision, parameter, public :: morvel_SC_vel =  0.755d0
double precision, parameter, public :: morvel_SM_lon =  -73.5d0 ! Somalia (East Africa)
double precision, parameter, public :: morvel_SM_lat =   59.3d0
double precision, parameter, public :: morvel_SM_vel =  0.980d0
double precision, parameter, public :: morvel_SR_lon =  -75.8d0 ! Sur (East of Scotia Plate)
double precision, parameter, public :: morvel_SR_lat =   55.7d0
double precision, parameter, public :: morvel_SR_vel =  0.636d0
double precision, parameter, public :: morvel_SU_lon =  -78.0d0 ! Sunda
double precision, parameter, public :: morvel_SU_lat =   59.8d0
double precision, parameter, public :: morvel_SU_vel =  0.973d0
double precision, parameter, public :: morvel_SW_lon =  -42.4d0 ! Sandwich (East of Scotia Plate)
double precision, parameter, public :: morvel_SW_lat =   -3.8d0
double precision, parameter, public :: morvel_SW_vel =  1.444d0
double precision, parameter, public :: morvel_YZ_lon =  -82.4d0 ! Yangtze
double precision, parameter, public :: morvel_YZ_lat =   65.5d0
double precision, parameter, public :: morvel_YZ_vel =  0.968d0
double precision, parameter, public :: morvel_ITRF05_lon = -68.2d0 ! ITRF 2005 (Altamimi et al., 2007)
double precision, parameter, public :: morvel_ITRF05_lat =  63.4d0
double precision, parameter, public :: morvel_ITRF05_vel = 0.677d0

! MORVEL56
! Argus, D.F., Gordon, R.G., DeMets, C. (2011). Geologically current motion of 56 plates relative
!     to the no‐net‐rotation reference frame. Geochemistry Geophysics Geosystems 12.
! Rotations relative to the Pacific (PA) plate
double precision, parameter, public :: morvel56_PA_lon =    0.00d0 ! Pacific
double precision, parameter, public :: morvel56_PA_lat =    0.00d0
double precision, parameter, public :: morvel56_PA_vel =   0.000d0
double precision, parameter, public :: morvel56_AM_lon =  -82.68d0 ! Amur
double precision, parameter, public :: morvel56_AM_lat =   65.92d0
double precision, parameter, public :: morvel56_AM_vel =   0.929d0
double precision, parameter, public :: morvel56_AN_lon =  -78.53d0 ! Antarctica
double precision, parameter, public :: morvel56_AN_lat =   65.92d0
double precision, parameter, public :: morvel56_AN_vel =   0.887d0
double precision, parameter, public :: morvel56_AR_lon =  -33.23d0 ! Arabia
double precision, parameter, public :: morvel56_AR_lat =   60.02d0
double precision, parameter, public :: morvel56_AR_vel =   1.159d0
double precision, parameter, public :: morvel56_AU_lon =    6.33d0 ! Australia
double precision, parameter, public :: morvel56_AU_lat =   60.08d0
double precision, parameter, public :: morvel56_AU_vel =   1.079d0
double precision, parameter, public :: morvel56_CA_lon =  -77.48d0 ! Caribbean
double precision, parameter, public :: morvel56_CA_lat =   55.76d0
double precision, parameter, public :: morvel56_CA_vel =   0.905d0
double precision, parameter, public :: morvel56_CO_lon = -112.78d0 ! Cocos
double precision, parameter, public :: morvel56_CO_lat =   42.18d0
double precision, parameter, public :: morvel56_CO_vel =   1.676d0
double precision, parameter, public :: morvel56_CP_lon =  -10.12d0 ! Capricorn
double precision, parameter, public :: morvel56_CP_lat =   62.34d0
double precision, parameter, public :: morvel56_CP_vel =   1.139d0
double precision, parameter, public :: morvel56_EU_lon =  -78.87d0 ! Eurasia
double precision, parameter, public :: morvel56_EU_lat =   61.27d0
double precision, parameter, public :: morvel56_EU_vel =   0.856d0
double precision, parameter, public :: morvel56_IN_lon =  -31.21d0 ! India
double precision, parameter, public :: morvel56_IN_lat =   61.39d0
double precision, parameter, public :: morvel56_IN_vel =   1.141d0
double precision, parameter, public :: morvel56_JF_lon =   37.84d0 ! Juan de Fuca
double precision, parameter, public :: morvel56_JF_lat =   -0.62d0
double precision, parameter, public :: morvel56_JF_vel =   0.625d0
double precision, parameter, public :: morvel56_LW_lon =  -66.90d0 ! Lwandle
double precision, parameter, public :: morvel56_LW_lat =   60.03d0
double precision, parameter, public :: morvel56_LW_vel =   0.932d0
double precision, parameter, public :: morvel56_MQ_lon =   -7.98d0 ! Macquarie
double precision, parameter, public :: morvel56_MQ_lat =   59.21d0
double precision, parameter, public :: morvel56_MQ_vel =   1.686d0
double precision, parameter, public :: morvel56_NA_lon =  -71.71d0 ! North America
double precision, parameter, public :: morvel56_NA_lat =   48.89d0
double precision, parameter, public :: morvel56_NA_vel =   0.750d0
double precision, parameter, public :: morvel56_NU_lon =  -66.57d0 ! Nubia
double precision, parameter, public :: morvel56_NU_lat =   58.68d0
double precision, parameter, public :: morvel56_NU_vel =   0.935d0
double precision, parameter, public :: morvel56_NZ_lon =  -87.76d0 ! Nazca
double precision, parameter, public :: morvel56_NZ_lat =   55.86d0
double precision, parameter, public :: morvel56_NZ_vel =   1.311d0
double precision, parameter, public :: morvel56_PS_lon =  -41.87d0 ! Philippine Sea
double precision, parameter, public :: morvel56_PS_lat =   -4.63d0
double precision, parameter, public :: morvel56_PS_vel =   0.890d0
double precision, parameter, public :: morvel56_RI_lon = -104.80d0 ! Rivera
double precision, parameter, public :: morvel56_RI_lat =   25.69d0
double precision, parameter, public :: morvel56_RI_vel =   4.966d0
double precision, parameter, public :: morvel56_SA_lon =  -77.03d0 ! South America
double precision, parameter, public :: morvel56_SA_lat =   55.97d0
double precision, parameter, public :: morvel56_SA_vel =   0.653d0
double precision, parameter, public :: morvel56_SC_lon =  -78.02d0 ! Scotia
double precision, parameter, public :: morvel56_SC_lat =   57.84d0
double precision, parameter, public :: morvel56_SC_vel =   0.755d0
double precision, parameter, public :: morvel56_SM_lon =  -73.55d0 ! Somalia
double precision, parameter, public :: morvel56_SM_lat =   59.27d0
double precision, parameter, public :: morvel56_SM_vel =   0.980d0
double precision, parameter, public :: morvel56_SR_lon =  -75.77d0 ! Sur
double precision, parameter, public :: morvel56_SR_lat =   55.69d0
double precision, parameter, public :: morvel56_SR_vel =   0.636d0
double precision, parameter, public :: morvel56_SU_lon =  -77.96d0 ! Sunda
double precision, parameter, public :: morvel56_SU_lat =   59.81d0
double precision, parameter, public :: morvel56_SU_vel =   0.973d0
double precision, parameter, public :: morvel56_SW_lon =  -42.36d0 ! Sandwich
double precision, parameter, public :: morvel56_SW_lat =   -3.84d0
double precision, parameter, public :: morvel56_SW_vel =   1.444d0
double precision, parameter, public :: morvel56_YZ_lon =  -82.38d0 ! Yangtze
double precision, parameter, public :: morvel56_YZ_lat =   65.45d0
double precision, parameter, public :: morvel56_YZ_vel =   0.968d0
double precision, parameter, public :: morvel56_AS_lon =  -70.76d0 ! Aegean Sea
double precision, parameter, public :: morvel56_AS_lat =   74.36d0
double precision, parameter, public :: morvel56_AS_vel =   0.648d0
double precision, parameter, public :: morvel56_AP_lon =  -77.01d0 ! Altiplano
double precision, parameter, public :: morvel56_AP_lat =   34.56d0
double precision, parameter, public :: morvel56_AP_vel =   0.929d0
double precision, parameter, public :: morvel56_AT_lon =    9.12d0 ! Anatolia
double precision, parameter, public :: morvel56_AT_lat =   54.82d0
double precision, parameter, public :: morvel56_AT_vel =   1.667d0
double precision, parameter, public :: morvel56_BR_lon = -111.00d0 ! Balmoral Reef
double precision, parameter, public :: morvel56_BR_lat =   45.90d0
double precision, parameter, public :: morvel56_BR_vel =   0.200d0
double precision, parameter, public :: morvel56_BS_lon =  122.56d0 ! Banda Sea
double precision, parameter, public :: morvel56_BS_lat =   13.34d0
double precision, parameter, public :: morvel56_BS_vel =   2.248d0
double precision, parameter, public :: morvel56_BH_lon =   88.39d0 ! Birds Head
double precision, parameter, public :: morvel56_BH_lat =   11.59d0
double precision, parameter, public :: morvel56_BH_vel =   0.346d0
double precision, parameter, public :: morvel56_BU_lon =  -76.63d0 ! Burma
double precision, parameter, public :: morvel56_BU_lat =    7.86d0
double precision, parameter, public :: morvel56_BU_vel =   2.523d0
double precision, parameter, public :: morvel56_CL_lon =  -27.64d0 ! Caroline Sea
double precision, parameter, public :: morvel56_CL_lat =    0.99d0
double precision, parameter, public :: morvel56_CL_vel =   0.199d0
double precision, parameter, public :: morvel56_CR_lon =  174.43d0 ! Conway Reef
double precision, parameter, public :: morvel56_CR_lat =  -12.56d0
double precision, parameter, public :: morvel56_CR_vel =   3.609d0
double precision, parameter, public :: morvel56_EA_lon =   66.32d0 ! Easter
double precision, parameter, public :: morvel56_EA_lat =   28.04d0
double precision, parameter, public :: morvel56_EA_vel =  11.420d0
double precision, parameter, public :: morvel56_FT_lon = -178.82d0 ! Futuna
double precision, parameter, public :: morvel56_FT_lat =  -10.12d0
double precision, parameter, public :: morvel56_FT_vel =   4.847d0
double precision, parameter, public :: morvel56_GP_lon =   79.43d0 ! Galapagos
double precision, parameter, public :: morvel56_GP_lat =    8.94d0
double precision, parameter, public :: morvel56_GP_vel =   5.307d0
double precision, parameter, public :: morvel56_JZ_lon =   70.11d0 ! Juan Fernandez
double precision, parameter, public :: morvel56_JZ_lat =   35.77d0
double precision, parameter, public :: morvel56_JZ_vel =  22.532d0
double precision, parameter, public :: morvel56_KE_lon =   -1.83d0 ! Kermadec
double precision, parameter, public :: morvel56_KE_lat =   47.61d0
double precision, parameter, public :: morvel56_KE_vel =   2.832d0
double precision, parameter, public :: morvel56_MN_lon =  150.46d0 ! Manus
double precision, parameter, public :: morvel56_MN_lat =   -3.04d0
double precision, parameter, public :: morvel56_MN_vel =  51.300d0
double precision, parameter, public :: morvel56_MO_lon =   79.96d0 ! Maoke
double precision, parameter, public :: morvel56_MO_lat =   57.43d0
double precision, parameter, public :: morvel56_MO_vel =   0.918d0
double precision, parameter, public :: morvel56_MA_lon =  144.24d0 ! Mariana
double precision, parameter, public :: morvel56_MA_lat =   39.19d0
double precision, parameter, public :: morvel56_MA_vel =   1.319d0
double precision, parameter, public :: morvel56_MS_lon =  -56.78d0 ! Molucca Sea
double precision, parameter, public :: morvel56_MS_lat =   10.54d0
double precision, parameter, public :: morvel56_MS_vel =   3.915d0
double precision, parameter, public :: morvel56_NH_lon =  -12.00d0 ! New Hebrides
double precision, parameter, public :: morvel56_NH_lat =   13.00d0
double precision, parameter, public :: morvel56_NH_vel =   2.700d0
double precision, parameter, public :: morvel56_NI_lon = -169.62d0 ! Niuafo'ou
double precision, parameter, public :: morvel56_NI_lat =    6.95d0
double precision, parameter, public :: morvel56_NI_vel =   3.248d0
double precision, parameter, public :: morvel56_ND_lon =  -80.24d0 ! North Andes
double precision, parameter, public :: morvel56_ND_lat =   59.68d0
double precision, parameter, public :: morvel56_ND_vel =   0.716d0
double precision, parameter, public :: morvel56_NB_lon =  139.00d0 ! North Bismarck
double precision, parameter, public :: morvel56_NB_lat =   -4.00d0
double precision, parameter, public :: morvel56_NB_vel =   0.330d0
double precision, parameter, public :: morvel56_OK_lon =  -76.20d0 ! Okhotsk
double precision, parameter, public :: morvel56_OK_lat =   55.81d0
double precision, parameter, public :: morvel56_OK_vel =   0.844d0
double precision, parameter, public :: morvel56_ON_lon =  141.58d0 ! Okinawa
double precision, parameter, public :: morvel56_ON_lat =   49.30d0
double precision, parameter, public :: morvel56_ON_vel =   2.743d0
double precision, parameter, public :: morvel56_PM_lon =  -88.73d0 ! Panama
double precision, parameter, public :: morvel56_PM_lat =   55.66d0
double precision, parameter, public :: morvel56_PM_vel =   0.906d0
double precision, parameter, public :: morvel56_SL_lon =  -92.39d0 ! Shetland
double precision, parameter, public :: morvel56_SL_lat =   65.24d0
double precision, parameter, public :: morvel56_SL_vel =   0.870d0
double precision, parameter, public :: morvel56_SS_lon =  133.82d0 ! Solomon Sea
double precision, parameter, public :: morvel56_SS_lat =   19.26d0
double precision, parameter, public :: morvel56_SS_vel =   1.509d0
double precision, parameter, public :: morvel56_SB_lon =  -32.99d0 ! South Bismarck
double precision, parameter, public :: morvel56_SB_lat =   10.61d0
double precision, parameter, public :: morvel56_SB_vel =   8.440d0
double precision, parameter, public :: morvel56_TI_lon =  113.28d0 ! Timor
double precision, parameter, public :: morvel56_TI_lat =   15.62d0
double precision, parameter, public :: morvel56_TI_vel =   1.629d0
double precision, parameter, public :: morvel56_TO_lon =    2.57d0 ! Tonga
double precision, parameter, public :: morvel56_TO_lat =   28.82d0
double precision, parameter, public :: morvel56_TO_vel =   9.303d0
double precision, parameter, public :: morvel56_WL_lon =  131.23d0 ! Woodlark
double precision, parameter, public :: morvel56_WL_lat =   21.81d0
double precision, parameter, public :: morvel56_WL_vel =   1.578d0
double precision, parameter, public :: morvel56_NNR_lon = -65.30d0 ! No net rotation
double precision, parameter, public :: morvel56_NNR_lat =  63.58d0
double precision, parameter, public :: morvel56_NNR_vel =  0.651d0
double precision, parameter, public :: morvel56_ITRF08_lon = -69.696d0 ! ITRF 2008
double precision, parameter, public :: morvel56_ITRF08_lat =  64.232d0 ! computed with subroutine below
double precision, parameter, public :: morvel56_ITRF08_vel =   0.649d0


! NUVEL-1A
! DeMets, C., Gordon, R.G., Argus, D.F., Stein, S. (1994). Effect of recent revisions to the
!     geomagnetic reversal time scale on estimates of current plate motions. Geophysical Research
!     Letters 21, 2191-2194.
! Rotations relative to the Pacific (PA) plate
double precision, parameter, public :: nuvel1a_PA_lon =    0.000d0 ! Pacific
double precision, parameter, public :: nuvel1a_PA_lat =    0.000d0
double precision, parameter, public :: nuvel1a_PA_vel =   0.0000d0
double precision, parameter, public :: nuvel1a_AF_lon =  -73.174d0 ! Africa
double precision, parameter, public :: nuvel1a_AF_lat =   59.160d0
double precision, parameter, public :: nuvel1a_AF_vel =   0.9270d0
double precision, parameter, public :: nuvel1a_AN_lon =  -83.984d0 ! Antarctica
double precision, parameter, public :: nuvel1a_AN_lat =   64.315d0
double precision, parameter, public :: nuvel1a_AN_vel =   0.8695d0
double precision, parameter, public :: nuvel1a_AR_lon =  -33.193d0 ! Arabia
double precision, parameter, public :: nuvel1a_AR_lat =   59.658d0
double precision, parameter, public :: nuvel1a_AR_vel =   1.1107d0
double precision, parameter, public :: nuvel1a_AU_lon =    1.742d0 ! Australia
double precision, parameter, public :: nuvel1a_AU_lat =   60.080d0
double precision, parameter, public :: nuvel1a_AU_vel =   1.0744d0
double precision, parameter, public :: nuvel1a_CA_lon =  -80.802d0 ! Caribbean
double precision, parameter, public :: nuvel1a_CA_lat =   54.195d0
double precision, parameter, public :: nuvel1a_CA_vel =   0.8160d0
double precision, parameter, public :: nuvel1a_CO_lon =  251.371d0 ! Cocos
double precision, parameter, public :: nuvel1a_CO_lat =   36.823d0
double precision, parameter, public :: nuvel1a_CO_vel =   1.9975d0
double precision, parameter, public :: nuvel1a_EU_lon =  -85.819d0 ! Eurasia
double precision, parameter, public :: nuvel1a_EU_lat =   61.066d0
double precision, parameter, public :: nuvel1a_EU_vel =   0.8591d0
double precision, parameter, public :: nuvel1a_IN_lon =  -30.403d0 ! India
double precision, parameter, public :: nuvel1a_IN_lat =   60.494d0
double precision, parameter, public :: nuvel1a_IN_vel =   1.1034d0
double precision, parameter, public :: nuvel1a_JF_lon =     29.3d0 ! Juan de Fuca
double precision, parameter, public :: nuvel1a_JF_lat =     28.3d0
double precision, parameter, public :: nuvel1a_JF_vel =    0.520d0
double precision, parameter, public :: nuvel1a_NA_lon =  -78.167d0 ! North America
double precision, parameter, public :: nuvel1a_NA_lat =   48.709d0
double precision, parameter, public :: nuvel1a_NA_vel =   0.7486d0
double precision, parameter, public :: nuvel1a_NZ_lon =  -90.096d0 ! Nazca
double precision, parameter, public :: nuvel1a_NZ_lat =   55.578d0
double precision, parameter, public :: nuvel1a_NZ_vel =   1.3599d0
double precision, parameter, public :: nuvel1a_PS_lon =     -1.2d0 ! Philippine Sea
double precision, parameter, public :: nuvel1a_PS_lat =    -45.8d0
double precision, parameter, public :: nuvel1a_PS_vel =     0.96d0
double precision, parameter, public :: nuvel1a_RI_lon =    257.6d0 ! Rivera
double precision, parameter, public :: nuvel1a_RI_lat =     31.0d0
double precision, parameter, public :: nuvel1a_RI_vel =     2.45d0
double precision, parameter, public :: nuvel1a_SA_lon =  -85.752d0 ! South America
double precision, parameter, public :: nuvel1a_SA_lat =   54.999d0
double precision, parameter, public :: nuvel1a_SA_vel =   0.6365d0
double precision, parameter, public :: nuvel1a_SC_lon =    -81.4d0 ! Scotia
double precision, parameter, public :: nuvel1a_SC_lat =     49.1d0
double precision, parameter, public :: nuvel1a_SC_vel =     0.66d0
double precision, parameter, public :: nuvel1a_NNR_lon =   -72.6d0 ! No net rotation
double precision, parameter, public :: nuvel1a_NNR_lat =    63.0d0
double precision, parameter, public :: nuvel1a_NNR_vel =  0.6411d0

!--------------------------------------------------------------------------------------------------!
contains
!--------------------------------------------------------------------------------------------------!

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------- EARTH SHAPE SUBROUTINES ------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine sphere_xyz2geo(x,y,z,lon,lat,height,radius)
!----
! Given a point in Cartesian coordinates with origin at the center of a sphere, calculate the
! longitude, latitude, and height above the surface of the sphere.
!----

use trig, only: r2d

implicit none

! Arguments
double precision :: x, y, z, lon, lat, height, radius


! Longitude depends only on x and y
lon = atan2(y,x)*r2d

! Latitude depends on z and the horizontal distance, sqrt(x^2+y^2)
lat = atan2(z,sqrt(x*x+y*y))*r2d

! Height is the total distance from origin minus the radius of the sphere
height = sqrt(x*x+y*y+z*z) - radius

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine sphere_geo2xyz(lon,lat,height,x,y,z,radius)
!----
! Given a point in geographic coordinates and its height above the surface of a sphere, calculate
! its Cartesian coordinates (origin at the center of the sphere).
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: lon, lat, height, x, y, z, radius

! Local variables
double precision :: R


R = radius + height

x = R*cos(lat*d2r)*cos(lon*d2r)
y = R*cos(lat*d2r)*sin(lon*d2r)
z = R*sin(lat*d2r)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine ellipsoid_xyz2geo(x,y,z,lon,lat,height,semimajor,eccentricity,iterations)
!----
! Given a point in Cartesian coordinates with origin at the center of an ellipsoid, calculate the
! longitude, latitude, and perpendicular height above the surface of the ellipsoid.
!
! Equations from:
!     Osborne, P. (2008). The Mercator Projections. 212 pp.
!----

use trig, only: r2d

implicit none

! Arguments
double precision :: x, y, z, lon, lat, height, semimajor, eccentricity
integer :: iterations

! Local variables
double precision :: nu
integer :: i


! The longitude is calculated the same way as for a sphere
lon = atan2(y,x)*r2d

! The latitude is computed iteratively
! Start with the latitude on the surface of the ellipsoid
lat = atan2(z,(1.0d0-eccentricity*eccentricity)*sqrt(x*x+y*y))

! Distance between point on ellipsoid and central axis of ellipsoid, perpendicular to surface
nu = semimajor/sqrt(1.0d0-eccentricity*eccentricity*sin(lat)*sin(lat))

! Iterate over latitude
do i = 1,iterations
    lat = atan2(z+eccentricity*eccentricity*nu*sin(lat),sqrt(x*x+y*y))
    nu = semimajor/sqrt(1.0d0-eccentricity*eccentricity*sin(lat)*sin(lat))
enddo

! Elevation perpendicular to ellipsoid surface
height = 1.0d0/cos(lat)*sqrt(x*x+y*y)-nu

lat = lat*r2d

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine ellipsoid_geo2xyz(lon,lat,height,x,y,z,semimajor,eccentricity)
!----
! Given a point in geographic coordinates and its perpendicular height above the surface of an
! ellipsoid, calculate its Cartesian coordinates (origin at the center of the ellipsoid).
!
! Equations from:
!     Osborne, P. (2008). The Mercator Projections. 212 pp.
!----

use trig, only: d2r

implicit none

! Arguments
double precision :: lon, lat, height, x, y, z, semimajor, eccentricity

! Local variables
double precision :: nu


nu = semimajor/sqrt(1.0d0-eccentricity*eccentricity*sin(lat*d2r)*sin(lat*d2r))

x = (                                  nu+height)*cos(lat*d2r)*cos(lon*d2r)
y = (                                  nu+height)*cos(lat*d2r)*sin(lon*d2r)
z = ((1.0d0-eccentricity*eccentricity)*nu+height)*sin(lat*d2r)

return
end subroutine

!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!
!------------------------------ PLATE CIRCUIT SUBROUTINES -----------------------------------------!
!--------------------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------------------!

subroutine pole_xyz2geo(pole_x,pole_y,pole_z,pole_lon,pole_lat,pole_vel,earth_shape)
!----
! Given a rotation vector in Cartesian coordinates, calculate the Euler pole longitude, latitude,
! and angular velocity (in the same units as input parameters).
!----

use io, only: stderr
use trig, only: r2d
use algebra, only: normalize

implicit none

! Arguments
double precision :: pole_x, pole_y, pole_z, pole_lon, pole_lat, pole_vel
character(len=*) :: earth_shape

! Local variables
double precision :: vec(3)


pole_lon = 0.0d0
pole_lat = 0.0d0
pole_vel = 0.0d0

if (earth_shape.eq.'sphere') then
    pole_lon = atan2(pole_y,pole_x)*r2d
    pole_lat = atan2(pole_z,sqrt(pole_x*pole_x+pole_y*pole_y))*r2d
    pole_vel = sqrt(pole_x*pole_x+pole_y*pole_y+pole_z*pole_z)
elseif (earth_shape.eq.'grs80') then
    ! I am not sure this is right....
    vec(1) = pole_x
    vec(2) = pole_y
    vec(3) = pole_z
    call normalize(vec)
    vec = vec*radius_earth_m
    call ellipsoid_xyz2geo(vec(1),vec(2),vec(3),pole_lon,pole_lat,pole_vel,semimajor_grs80, &
                           eccentricity_grs80,0)
else
    write(stderr,*) 'pole_xyz2geo: no earth shape named "',trim(earth_shape),'"; returning zeros'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine pole_geo2xyz(pole_lon,pole_lat,pole_vel,pole_x,pole_y,pole_z,earth_shape)
!----
! Given a rotation vector in geographic coordinates and a rotational velocity, calculate the Euler
! pole in Cartesian coordinates (in the same units as input rotational velocity).
!----

use io, only: stderr
use trig, only: d2r
use algebra, only: normalize

implicit none

! Arguments
double precision :: pole_lon, pole_lat, pole_vel, pole_x, pole_y, pole_z
character(len=*) :: earth_shape

! Local variables
double precision :: vec(3)


if (earth_shape.eq.'sphere') then
    pole_x = pole_vel*cos(pole_lat*d2r)*cos(pole_lon*d2r)
    pole_y = pole_vel*cos(pole_lat*d2r)*sin(pole_lon*d2r)
    pole_z = pole_vel*sin(pole_lat*d2r)
elseif (earth_shape.eq.'grs80') then
    call ellipsoid_geo2xyz(pole_lon,pole_lat,0.0d0,vec(1),vec(2),vec(3),1.0d0,eccentricity_grs80)
    call normalize(vec)
    pole_x = pole_vel*vec(1)
    pole_y = pole_vel*vec(2)
    pole_z = pole_vel*vec(3)
else
    write(stderr,*) 'pole_geo2xyz: no earth shape named "',trim(earth_shape),'"; returning zeros'
endif

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine itrf2008_wrt_morvel56_PA()
!----
! Calculate the angular velocity of ITRF 2008 (Altamimi et al., 2012) with respect to the Pacific
! plate of MORVEL56 (Argus et al., 2011).
!----
implicit none

! Local variables
double precision :: NNRM56_ITRF08_x, NNRM56_ITRF08_y, NNRM56_ITRF08_z
double precision :: PAM56_NNRM56_x, PAM56_NNRM56_y, PAM56_NNRM56_z
double precision :: PAM56_ITRF2008_x, PAM56_ITRF2008_y, PAM56_ITRF2008_z
double precision :: morvel56_ITRF08_lon_calc, morvel56_ITRF08_lat_calc, morvel56_ITRF08_vel_calc

! NNR-MORVEL56 w.r.t. ITRF 2008
NNRM56_ITRF08_x =  0.083d0*masa2degma
NNRM56_ITRF08_y =  0.006d0*masa2degma
NNRM56_ITRF08_z = -0.007d0*masa2degma

! MORVEL56 PA w.r.t. NNR-MORVEL56
call pole_geo2xyz(morvel56_NNR_lon,morvel56_NNR_lat,-morvel56_NNR_vel,&
                  PAM56_NNRM56_x,PAM56_NNRM56_y,PAM56_NNRM56_z,'sphere')

! MORVEL56 PA w.r.t. ITRF 2008
PAM56_ITRF2008_x = PAM56_NNRM56_x + NNRM56_ITRF08_x
PAM56_ITRF2008_y = PAM56_NNRM56_y + NNRM56_ITRF08_y
PAM56_ITRF2008_z = PAM56_NNRM56_z + NNRM56_ITRF08_z

! ITRF 2008 w.r.t. MORVEL56 PA
call pole_xyz2geo(-PAM56_ITRF2008_x,-PAM56_ITRF2008_y,-PAM56_ITRF2008_z,&
                  morvel56_ITRF08_lon_calc,morvel56_ITRF08_lat_calc,morvel56_ITRF08_vel_calc,'sphere')

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine get_pole(fixed_plate,moving_plate,model,pole,ierr)
!----
! Calculate the pole of rotation between two plates
!----

use io, only: stderr

implicit none

! Arguments
character(len=*) :: fixed_plate, moving_plate, model
double precision :: pole(3)
integer :: ierr

! Local variables
integer :: i
double precision :: geo_pole(2,3), xyz_pole(2,3)
character(len=64) :: plate(2)


ierr = 0
plate(1) = fixed_plate
plate(2) = moving_plate

do i = 1,2
    if (model.eq.'ITRF05') then
        if (plate(i).eq.'AM') then
            geo_pole(i,1) = itrf05_AM_lon
            geo_pole(i,2) = itrf05_AM_lat
            geo_pole(i,3) = itrf05_AM_vel
        elseif (plate(i).eq.'AN') then
            geo_pole(i,1) = itrf05_AN_lon
            geo_pole(i,2) = itrf05_AN_lat
            geo_pole(i,3) = itrf05_AN_vel
        elseif (plate(i).eq.'AR') then
            geo_pole(i,1) = itrf05_AR_lon
            geo_pole(i,2) = itrf05_AR_lat
            geo_pole(i,3) = itrf05_AR_vel
        elseif (plate(i).eq.'AU') then
            geo_pole(i,1) = itrf05_AU_lon
            geo_pole(i,2) = itrf05_AU_lat
            geo_pole(i,3) = itrf05_AU_vel
        elseif (plate(i).eq.'CA') then
            geo_pole(i,1) = itrf05_CA_lon
            geo_pole(i,2) = itrf05_CA_lat
            geo_pole(i,3) = itrf05_CA_vel
        elseif (plate(i).eq.'EU') then
            geo_pole(i,1) = itrf05_EU_lon
            geo_pole(i,2) = itrf05_EU_lat
            geo_pole(i,3) = itrf05_EU_vel
        elseif (plate(i).eq.'IN') then
            geo_pole(i,1) = itrf05_IN_lon
            geo_pole(i,2) = itrf05_IN_lat
            geo_pole(i,3) = itrf05_IN_vel
        elseif (plate(i).eq.'NZ') then
            geo_pole(i,1) = itrf05_NZ_lon
            geo_pole(i,2) = itrf05_NZ_lat
            geo_pole(i,3) = itrf05_NZ_vel
        elseif (plate(i).eq.'NA') then
            geo_pole(i,1) = itrf05_NA_lon
            geo_pole(i,2) = itrf05_NA_lat
            geo_pole(i,3) = itrf05_NA_vel
        elseif (plate(i).eq.'NU') then
            geo_pole(i,1) = itrf05_NU_lon
            geo_pole(i,2) = itrf05_NU_lat
            geo_pole(i,3) = itrf05_NU_vel
        elseif (plate(i).eq.'OK') then
            geo_pole(i,1) = itrf05_OK_lon
            geo_pole(i,2) = itrf05_OK_lat
            geo_pole(i,3) = itrf05_OK_vel
        elseif (plate(i).eq.'PA') then
            geo_pole(i,1) = itrf05_PA_lon
            geo_pole(i,2) = itrf05_PA_lat
            geo_pole(i,3) = itrf05_PA_vel
        elseif (plate(i).eq.'SA') then
            geo_pole(i,1) = itrf05_SA_lon
            geo_pole(i,2) = itrf05_SA_lat
            geo_pole(i,3) = itrf05_SA_vel
        elseif (plate(i).eq.'SM') then
            geo_pole(i,1) = itrf05_SM_lon
            geo_pole(i,2) = itrf05_SM_lat
            geo_pole(i,3) = itrf05_SM_vel
        elseif (plate(i).eq.'YZ') then
            geo_pole(i,1) = itrf05_YZ_lon
            geo_pole(i,2) = itrf05_YZ_lat
            geo_pole(i,3) = itrf05_YZ_vel
        elseif (plate(i).eq.'ITRF05') then
            geo_pole(i,1) = itrf05_ITRF05_lon
            geo_pole(i,2) = itrf05_ITRF05_lat
            geo_pole(i,3) = itrf05_ITRF05_vel
        else
            write(stderr,*) 'get_pole: no plate named "',trim(plate(i)),'" in ITRF05 model'
            call list_plate_models('ITRF05')
            ierr = 1
            return
        endif

    elseif (model.eq.'ITRF08') then
        if (plate(i).eq.'AM') then
            geo_pole(i,1) = itrf08_AM_lon
            geo_pole(i,2) = itrf08_AM_lat
            geo_pole(i,3) = itrf08_AM_vel
        elseif (plate(i).eq.'AN') then
            geo_pole(i,1) = itrf08_AN_lon
            geo_pole(i,2) = itrf08_AN_lat
            geo_pole(i,3) = itrf08_AN_vel
        elseif (plate(i).eq.'AR') then
            geo_pole(i,1) = itrf08_AR_lon
            geo_pole(i,2) = itrf08_AR_lat
            geo_pole(i,3) = itrf08_AR_vel
        elseif (plate(i).eq.'AU') then
            geo_pole(i,1) = itrf08_AU_lon
            geo_pole(i,2) = itrf08_AU_lat
            geo_pole(i,3) = itrf08_AU_vel
        elseif (plate(i).eq.'CA') then
            geo_pole(i,1) = itrf08_CA_lon
            geo_pole(i,2) = itrf08_CA_lat
            geo_pole(i,3) = itrf08_CA_vel
        elseif (plate(i).eq.'EU') then
            geo_pole(i,1) = itrf08_EU_lon
            geo_pole(i,2) = itrf08_EU_lat
            geo_pole(i,3) = itrf08_EU_vel
        elseif (plate(i).eq.'IN') then
            geo_pole(i,1) = itrf08_IN_lon
            geo_pole(i,2) = itrf08_IN_lat
            geo_pole(i,3) = itrf08_IN_vel
        elseif (plate(i).eq.'NZ') then
            geo_pole(i,1) = itrf08_NZ_lon
            geo_pole(i,2) = itrf08_NZ_lat
            geo_pole(i,3) = itrf08_NZ_vel
        elseif (plate(i).eq.'NA') then
            geo_pole(i,1) = itrf08_NA_lon
            geo_pole(i,2) = itrf08_NA_lat
            geo_pole(i,3) = itrf08_NA_vel
        elseif (plate(i).eq.'NU') then
            geo_pole(i,1) = itrf08_NU_lon
            geo_pole(i,2) = itrf08_NU_lat
            geo_pole(i,3) = itrf08_NU_vel
        elseif (plate(i).eq.'PA') then
            geo_pole(i,1) = itrf08_PA_lon
            geo_pole(i,2) = itrf08_PA_lat
            geo_pole(i,3) = itrf08_PA_vel
        elseif (plate(i).eq.'SA') then
            geo_pole(i,1) = itrf08_SA_lon
            geo_pole(i,2) = itrf08_SA_lat
            geo_pole(i,3) = itrf08_SA_vel
        elseif (plate(i).eq.'SM') then
            geo_pole(i,1) = itrf08_SM_lon
            geo_pole(i,2) = itrf08_SM_lat
            geo_pole(i,3) = itrf08_SM_vel
        elseif (plate(i).eq.'SU') then
            geo_pole(i,1) = itrf08_SU_lon
            geo_pole(i,2) = itrf08_SU_lat
            geo_pole(i,3) = itrf08_SU_vel
        elseif (plate(i).eq.'ITRF05') then
            geo_pole(i,1) = itrf08_ITRF05_lon
            geo_pole(i,2) = itrf08_ITRF05_lat
            geo_pole(i,3) = itrf08_ITRF05_vel
        elseif (plate(i).eq.'ITRF08') then
            geo_pole(i,1) = itrf08_ITRF08_lon
            geo_pole(i,2) = itrf08_ITRF08_lat
            geo_pole(i,3) = itrf08_ITRF08_vel
        else
            write(stderr,*) 'get_pole: no plate named "',trim(plate(i)),'" in ITRF08 model'
            call list_plate_models('ITRF08')
            ierr = 1
            return
        endif

    elseif (model.eq.'MORVEL') then
        if (plate(i).eq.'PA') then
            geo_pole(i,1) = morvel_PA_lon
            geo_pole(i,2) = morvel_PA_lat
            geo_pole(i,3) = morvel_PA_vel
        elseif (plate(i).eq.'AM') then
            geo_pole(i,1) = morvel_AM_lon
            geo_pole(i,2) = morvel_AM_lat
            geo_pole(i,3) = morvel_AM_vel
        elseif (plate(i).eq.'AN') then
            geo_pole(i,1) = morvel_AN_lon
            geo_pole(i,2) = morvel_AN_lat
            geo_pole(i,3) = morvel_AN_vel
        elseif (plate(i).eq.'AR') then
            geo_pole(i,1) = morvel_AR_lon
            geo_pole(i,2) = morvel_AR_lat
            geo_pole(i,3) = morvel_AR_vel
        elseif (plate(i).eq.'AU') then
            geo_pole(i,1) = morvel_AU_lon
            geo_pole(i,2) = morvel_AU_lat
            geo_pole(i,3) = morvel_AU_vel
        elseif (plate(i).eq.'CA') then
            geo_pole(i,1) = morvel_CA_lon
            geo_pole(i,2) = morvel_CA_lat
            geo_pole(i,3) = morvel_CA_vel
        elseif (plate(i).eq.'CO') then
            geo_pole(i,1) = morvel_CO_lon
            geo_pole(i,2) = morvel_CO_lat
            geo_pole(i,3) = morvel_CO_vel
        elseif (plate(i).eq.'CP') then
            geo_pole(i,1) = morvel_CP_lon
            geo_pole(i,2) = morvel_CP_lat
            geo_pole(i,3) = morvel_CP_vel
        elseif (plate(i).eq.'EU') then
            geo_pole(i,1) = morvel_EU_lon
            geo_pole(i,2) = morvel_EU_lat
            geo_pole(i,3) = morvel_EU_vel
        elseif (plate(i).eq.'IN') then
            geo_pole(i,1) = morvel_IN_lon
            geo_pole(i,2) = morvel_IN_lat
            geo_pole(i,3) = morvel_IN_vel
        elseif (plate(i).eq.'JF') then
            geo_pole(i,1) = morvel_JF_lon
            geo_pole(i,2) = morvel_JF_lat
            geo_pole(i,3) = morvel_JF_vel
        elseif (plate(i).eq.'LW') then
            geo_pole(i,1) = morvel_LW_lon
            geo_pole(i,2) = morvel_LW_lat
            geo_pole(i,3) = morvel_LW_vel
        elseif (plate(i).eq.'MQ') then
            geo_pole(i,1) = morvel_MQ_lon
            geo_pole(i,2) = morvel_MQ_lat
            geo_pole(i,3) = morvel_MQ_vel
        elseif (plate(i).eq.'NA') then
            geo_pole(i,1) = morvel_NA_lon
            geo_pole(i,2) = morvel_NA_lat
            geo_pole(i,3) = morvel_NA_vel
        elseif (plate(i).eq.'NU') then
            geo_pole(i,1) = morvel_NU_lon
            geo_pole(i,2) = morvel_NU_lat
            geo_pole(i,3) = morvel_NU_vel
        elseif (plate(i).eq.'NZ') then
            geo_pole(i,1) = morvel_NZ_lon
            geo_pole(i,2) = morvel_NZ_lat
            geo_pole(i,3) = morvel_NZ_vel
        elseif (plate(i).eq.'PS') then
            geo_pole(i,1) = morvel_PS_lon
            geo_pole(i,2) = morvel_PS_lat
            geo_pole(i,3) = morvel_PS_vel
        elseif (plate(i).eq.'RI') then
            geo_pole(i,1) = morvel_RI_lon
            geo_pole(i,2) = morvel_RI_lat
            geo_pole(i,3) = morvel_RI_vel
        elseif (plate(i).eq.'SA') then
            geo_pole(i,1) = morvel_SA_lon
            geo_pole(i,2) = morvel_SA_lat
            geo_pole(i,3) = morvel_SA_vel
        elseif (plate(i).eq.'SC') then
            geo_pole(i,1) = morvel_SC_lon
            geo_pole(i,2) = morvel_SC_lat
            geo_pole(i,3) = morvel_SC_vel
        elseif (plate(i).eq.'SM') then
            geo_pole(i,1) = morvel_SM_lon
            geo_pole(i,2) = morvel_SM_lat
            geo_pole(i,3) = morvel_SM_vel
        elseif (plate(i).eq.'SR') then
            geo_pole(i,1) = morvel_SR_lon
            geo_pole(i,2) = morvel_SR_lat
            geo_pole(i,3) = morvel_SR_vel
        elseif (plate(i).eq.'SU') then
            geo_pole(i,1) = morvel_SU_lon
            geo_pole(i,2) = morvel_SU_lat
            geo_pole(i,3) = morvel_SU_vel
        elseif (plate(i).eq.'SW') then
            geo_pole(i,1) = morvel_SW_lon
            geo_pole(i,2) = morvel_SW_lat
            geo_pole(i,2) = morvel_SW_vel
        elseif (plate(i).eq.'YZ') then
            geo_pole(i,1) = morvel_YZ_lon
            geo_pole(i,2) = morvel_YZ_lat
            geo_pole(i,3) = morvel_YZ_vel
        elseif (plate(i).eq.'ITRF05') then
            geo_pole(i,1) = morvel_ITRF05_lon
            geo_pole(i,2) = morvel_ITRF05_lat
            geo_pole(i,3) = morvel_ITRF05_vel
        else
            write(stderr,*) 'get_pole: no plate named "',trim(plate(i)),'" in MORVEL model'
            call list_plate_models('MORVEL')
            ierr = 1
            return
        endif

    elseif (model.eq.'MORVEL56'.or.model.eq.'MORVEL-56') then
        if (plate(i).eq.'PA') then
            geo_pole(i,1) = morvel56_PA_lon
            geo_pole(i,2) = morvel56_PA_lat
            geo_pole(i,3) = morvel56_PA_vel
        elseif (plate(i).eq.'AM') then
            geo_pole(i,1) = morvel56_AM_lon
            geo_pole(i,2) = morvel56_AM_lat
            geo_pole(i,3) = morvel56_AM_vel
        elseif (plate(i).eq.'AN') then
            geo_pole(i,1) = morvel56_AN_lon
            geo_pole(i,2) = morvel56_AN_lat
            geo_pole(i,3) = morvel56_AN_vel
        elseif (plate(i).eq.'AR') then
            geo_pole(i,1) = morvel56_AR_lon
            geo_pole(i,2) = morvel56_AR_lat
            geo_pole(i,3) = morvel56_AR_vel
        elseif (plate(i).eq.'AU') then
            geo_pole(i,1) = morvel56_AU_lon
            geo_pole(i,2) = morvel56_AU_lat
            geo_pole(i,3) = morvel56_AU_vel
        elseif (plate(i).eq.'CA') then
            geo_pole(i,1) = morvel56_CA_lon
            geo_pole(i,2) = morvel56_CA_lat
            geo_pole(i,3) = morvel56_CA_vel
        elseif (plate(i).eq.'CO') then
            geo_pole(i,1) = morvel56_CO_lon
            geo_pole(i,2) = morvel56_CO_lat
            geo_pole(i,3) = morvel56_CO_vel
        elseif (plate(i).eq.'CP') then
            geo_pole(i,1) = morvel56_CP_lon
            geo_pole(i,2) = morvel56_CP_lat
            geo_pole(i,3) = morvel56_CP_vel
        elseif (plate(i).eq.'EU') then
            geo_pole(i,1) = morvel56_EU_lon
            geo_pole(i,2) = morvel56_EU_lat
            geo_pole(i,3) = morvel56_EU_vel
        elseif (plate(i).eq.'IN') then
            geo_pole(i,1) = morvel56_IN_lon
            geo_pole(i,2) = morvel56_IN_lat
            geo_pole(i,3) = morvel56_IN_vel
        elseif (plate(i).eq.'JF') then
            geo_pole(i,1) = morvel56_JF_lon
            geo_pole(i,2) = morvel56_JF_lat
            geo_pole(i,3) = morvel56_JF_vel
        elseif (plate(i).eq.'LW') then
            geo_pole(i,1) = morvel56_LW_lon
            geo_pole(i,2) = morvel56_LW_lat
            geo_pole(i,3) = morvel56_LW_vel
        elseif (plate(i).eq.'MQ') then
            geo_pole(i,1) = morvel56_MQ_lon
            geo_pole(i,2) = morvel56_MQ_lat
            geo_pole(i,3) = morvel56_MQ_vel
        elseif (plate(i).eq.'NA') then
            geo_pole(i,1) = morvel56_NA_lon
            geo_pole(i,2) = morvel56_NA_lat
            geo_pole(i,3) = morvel56_NA_vel
        elseif (plate(i).eq.'NU') then
            geo_pole(i,1) = morvel56_NU_lon
            geo_pole(i,2) = morvel56_NU_lat
            geo_pole(i,3) = morvel56_NU_vel
        elseif (plate(i).eq.'NZ') then
            geo_pole(i,1) = morvel56_NZ_lon
            geo_pole(i,2) = morvel56_NZ_lat
            geo_pole(i,3) = morvel56_NZ_vel
        elseif (plate(i).eq.'PS') then
            geo_pole(i,1) = morvel56_PS_lon
            geo_pole(i,2) = morvel56_PS_lat
            geo_pole(i,3) = morvel56_PS_vel
        elseif (plate(i).eq.'RI') then
            geo_pole(i,1) = morvel56_RI_lon
            geo_pole(i,2) = morvel56_RI_lat
            geo_pole(i,3) = morvel56_RI_vel
        elseif (plate(i).eq.'SA') then
            geo_pole(i,1) = morvel56_SA_lon
            geo_pole(i,2) = morvel56_SA_lat
            geo_pole(i,3) = morvel56_SA_vel
        elseif (plate(i).eq.'SC') then
            geo_pole(i,1) = morvel56_SC_lon
            geo_pole(i,2) = morvel56_SC_lat
            geo_pole(i,3) = morvel56_SC_vel
        elseif (plate(i).eq.'SM') then
            geo_pole(i,1) = morvel56_SM_lon
            geo_pole(i,2) = morvel56_SM_lat
            geo_pole(i,3) = morvel56_SM_vel
        elseif (plate(i).eq.'SR') then
            geo_pole(i,1) = morvel56_SR_lon
            geo_pole(i,2) = morvel56_SR_lat
            geo_pole(i,3) = morvel56_SR_vel
        elseif (plate(i).eq.'SU') then
            geo_pole(i,1) = morvel56_SU_lon
            geo_pole(i,2) = morvel56_SU_lat
            geo_pole(i,3) = morvel56_SU_vel
        elseif (plate(i).eq.'SW') then
            geo_pole(i,1) = morvel56_SW_lon
            geo_pole(i,2) = morvel56_SW_lat
            geo_pole(i,3) = morvel56_SW_vel
        elseif (plate(i).eq.'YZ') then
            geo_pole(i,1) = morvel56_YZ_lon
            geo_pole(i,2) = morvel56_YZ_lat
            geo_pole(i,3) = morvel56_YZ_vel
        elseif (plate(i).eq.'AS') then
            geo_pole(i,1) = morvel56_AS_lon
            geo_pole(i,2) = morvel56_AS_lat
            geo_pole(i,3) = morvel56_AS_vel
        elseif (plate(i).eq.'AP') then
            geo_pole(i,1) = morvel56_AP_lon
            geo_pole(i,2) = morvel56_AP_lat
            geo_pole(i,3) = morvel56_AP_vel
        elseif (plate(i).eq.'AT') then
            geo_pole(i,1) = morvel56_AT_lon
            geo_pole(i,2) = morvel56_AT_lat
            geo_pole(i,3) = morvel56_AT_vel
        elseif (plate(i).eq.'BR') then
            geo_pole(i,1) = morvel56_BR_lon
            geo_pole(i,2) = morvel56_BR_lat
            geo_pole(i,3) = morvel56_BR_vel
        elseif (plate(i).eq.'BS') then
            geo_pole(i,1) = morvel56_BS_lon
            geo_pole(i,2) = morvel56_BS_lat
            geo_pole(i,3) = morvel56_BS_vel
        elseif (plate(i).eq.'BH') then
            geo_pole(i,1) = morvel56_BH_lon
            geo_pole(i,2) = morvel56_BH_lat
            geo_pole(i,3) = morvel56_BH_vel
        elseif (plate(i).eq.'BU') then
            geo_pole(i,1) = morvel56_BU_lon
            geo_pole(i,2) = morvel56_BU_lat
            geo_pole(i,3) = morvel56_BU_vel
        elseif (plate(i).eq.'CL') then
            geo_pole(i,1) = morvel56_CL_lon
            geo_pole(i,2) = morvel56_CL_lat
            geo_pole(i,3) = morvel56_CL_vel
        elseif (plate(i).eq.'CR') then
            geo_pole(i,1) = morvel56_CR_lon
            geo_pole(i,2) = morvel56_CR_lat
            geo_pole(i,3) = morvel56_CR_vel
        elseif (plate(i).eq.'EA') then
            geo_pole(i,1) = morvel56_EA_lon
            geo_pole(i,2) = morvel56_EA_lat
            geo_pole(i,3) = morvel56_EA_vel
        elseif (plate(i).eq.'FT') then
            geo_pole(i,1) = morvel56_FT_lon
            geo_pole(i,2) = morvel56_FT_lat
            geo_pole(i,3) = morvel56_FT_vel
        elseif (plate(i).eq.'GP') then
            geo_pole(i,1) = morvel56_GP_lon
            geo_pole(i,2) = morvel56_GP_lat
            geo_pole(i,3) = morvel56_GP_vel
        elseif (plate(i).eq.'JZ') then
            geo_pole(i,1) = morvel56_JZ_lon
            geo_pole(i,2) = morvel56_JZ_lat
            geo_pole(i,3) = morvel56_JZ_vel
        elseif (plate(i).eq.'KE') then
            geo_pole(i,1) = morvel56_KE_lon
            geo_pole(i,2) = morvel56_KE_lat
            geo_pole(i,3) = morvel56_KE_vel
        elseif (plate(i).eq.'MN') then
            geo_pole(i,1) = morvel56_MN_lon
            geo_pole(i,2) = morvel56_MN_lat
            geo_pole(i,3) = morvel56_MN_vel
        elseif (plate(i).eq.'MO') then
            geo_pole(i,1) = morvel56_MO_lon
            geo_pole(i,2) = morvel56_MO_lat
            geo_pole(i,3) = morvel56_MO_vel
        elseif (plate(i).eq.'MA') then
            geo_pole(i,1) = morvel56_MA_lon
            geo_pole(i,2) = morvel56_MA_lat
            geo_pole(i,3) = morvel56_MA_vel
        elseif (plate(i).eq.'MS') then
            geo_pole(i,1) = morvel56_MS_lon
            geo_pole(i,2) = morvel56_MS_lat
            geo_pole(i,3) = morvel56_MS_vel
        elseif (plate(i).eq.'NH') then
            geo_pole(i,1) = morvel56_NH_lon
            geo_pole(i,2) = morvel56_NH_lat
            geo_pole(i,3) = morvel56_NH_vel
        elseif (plate(i).eq.'NI') then
            geo_pole(i,1) = morvel56_NI_lon
            geo_pole(i,2) = morvel56_NI_lat
            geo_pole(i,3) = morvel56_NI_vel
        elseif (plate(i).eq.'ND') then
            geo_pole(i,1) = morvel56_ND_lon
            geo_pole(i,2) = morvel56_ND_lat
            geo_pole(i,3) = morvel56_ND_vel
        elseif (plate(i).eq.'NB') then
            geo_pole(i,1) = morvel56_NB_lon
            geo_pole(i,2) = morvel56_NB_lat
            geo_pole(i,3) = morvel56_NB_vel
        elseif (plate(i).eq.'OK') then
            geo_pole(i,1) = morvel56_OK_lon
            geo_pole(i,2) = morvel56_OK_lat
            geo_pole(i,3) = morvel56_OK_vel
        elseif (plate(i).eq.'ON') then
            geo_pole(i,1) = morvel56_ON_lon
            geo_pole(i,2) = morvel56_ON_lat
            geo_pole(i,3) = morvel56_ON_vel
        elseif (plate(i).eq.'PM') then
            geo_pole(i,1) = morvel56_PM_lon
            geo_pole(i,2) = morvel56_PM_lat
            geo_pole(i,3) = morvel56_PM_vel
        elseif (plate(i).eq.'SL') then
            geo_pole(i,1) = morvel56_SL_lon
            geo_pole(i,2) = morvel56_SL_lat
            geo_pole(i,3) = morvel56_SL_vel
        elseif (plate(i).eq.'SS') then
            geo_pole(i,1) = morvel56_SS_lon
            geo_pole(i,2) = morvel56_SS_lat
            geo_pole(i,3) = morvel56_SS_vel
        elseif (plate(i).eq.'SB') then
            geo_pole(i,1) = morvel56_SB_lon
            geo_pole(i,2) = morvel56_SB_lat
            geo_pole(i,3) = morvel56_SB_vel
        elseif (plate(i).eq.'TI') then
            geo_pole(i,1) = morvel56_TI_lon
            geo_pole(i,2) = morvel56_TI_lat
            geo_pole(i,3) = morvel56_TI_vel
        elseif (plate(i).eq.'TO') then
            geo_pole(i,1) = morvel56_TO_lon
            geo_pole(i,2) = morvel56_TO_lat
            geo_pole(i,3) = morvel56_TO_vel
        elseif (plate(i).eq.'WL') then
            geo_pole(i,1) = morvel56_WL_lon
            geo_pole(i,2) = morvel56_WL_lat
            geo_pole(i,3) = morvel56_WL_vel
        elseif (plate(i).eq.'NNR') then
            geo_pole(i,1) = morvel56_NNR_lon
            geo_pole(i,2) = morvel56_NNR_lat
            geo_pole(i,3) = morvel56_NNR_vel
        elseif (plate(i).eq.'ITRF08') then
            geo_pole(i,1) = morvel56_ITRF08_lon
            geo_pole(i,2) = morvel56_ITRF08_lat
            geo_pole(i,3) = morvel56_ITRF08_vel
        else
            write(stderr,*) 'get_pole: no plate named "',trim(plate(i)),'" in MORVEL56 model'
            call list_plate_models('MORVEL56')
            ierr = 1
            return
        endif

    elseif (model.eq.'NUVEL1A'.or.model.eq.'NUVEL-1A') then
        if (plate(i).eq.'PA') then
            geo_pole(i,1) = nuvel1a_PA_lon
            geo_pole(i,2) = nuvel1a_PA_lat
            geo_pole(i,3) = nuvel1a_PA_vel
        elseif (plate(i).eq.'AF') then
            geo_pole(i,1) = nuvel1a_AF_lon
            geo_pole(i,2) = nuvel1a_AF_lat
            geo_pole(i,3) = nuvel1a_AF_vel
        elseif (plate(i).eq.'AN') then
            geo_pole(i,1) = nuvel1a_AN_lon
            geo_pole(i,2) = nuvel1a_AN_lat
            geo_pole(i,3) = nuvel1a_AN_vel
        elseif (plate(i).eq.'AT') then
            geo_pole(i,1) = nuvel1a_AR_lon
            geo_pole(i,2) = nuvel1a_AR_lat
            geo_pole(i,3) = nuvel1a_AR_vel
        elseif (plate(i).eq.'AU') then
            geo_pole(i,1) = nuvel1a_AU_lon
            geo_pole(i,2) = nuvel1a_AU_lat
            geo_pole(i,3) = nuvel1a_AU_vel
        elseif (plate(i).eq.'CA') then
            geo_pole(i,1) = nuvel1a_CA_lon
            geo_pole(i,2) = nuvel1a_CA_lat
            geo_pole(i,3) = nuvel1a_CA_vel
        elseif (plate(i).eq.'CO') then
            geo_pole(i,1) = nuvel1a_CO_lon
            geo_pole(i,2) = nuvel1a_CO_lat
            geo_pole(i,3) = nuvel1a_CO_vel
        elseif (plate(i).eq.'EU') then
            geo_pole(i,1) = nuvel1a_EU_lon
            geo_pole(i,2) = nuvel1a_EU_lat
            geo_pole(i,3) = nuvel1a_EU_vel
        elseif (plate(i).eq.'IN') then
            geo_pole(i,1) = nuvel1a_IN_lon
            geo_pole(i,2) = nuvel1a_IN_lat
            geo_pole(i,3) = nuvel1a_IN_vel
        elseif (plate(i).eq.'JF') then
            geo_pole(i,1) = nuvel1a_JF_lon
            geo_pole(i,2) = nuvel1a_JF_lat
            geo_pole(i,3) = nuvel1a_JF_vel
        elseif (plate(i).eq.'NA') then
            geo_pole(i,1) = nuvel1a_NA_lon
            geo_pole(i,2) = nuvel1a_NA_lat
            geo_pole(i,3) = nuvel1a_NA_vel
        elseif (plate(i).eq.'NZ') then
            geo_pole(i,1) = nuvel1a_NZ_lon
            geo_pole(i,2) = nuvel1a_NZ_lat
            geo_pole(i,3) = nuvel1a_NZ_vel
        elseif (plate(i).eq.'PS') then
            geo_pole(i,1) = nuvel1a_PS_lon
            geo_pole(i,2) = nuvel1a_PS_lat
            geo_pole(i,3) = nuvel1a_PS_vel
        elseif (plate(i).eq.'RI') then
            geo_pole(i,1) = nuvel1a_RI_lon
            geo_pole(i,2) = nuvel1a_RI_lat
            geo_pole(i,3) = nuvel1a_RI_vel
        elseif (plate(i).eq.'SA') then
            geo_pole(i,1) = nuvel1a_SA_lon
            geo_pole(i,2) = nuvel1a_SA_lat
            geo_pole(i,3) = nuvel1a_SA_vel
        elseif (plate(i).eq.'SC') then
            geo_pole(i,1) = nuvel1a_SC_lon
            geo_pole(i,2) = nuvel1a_SC_lat
            geo_pole(i,3) = nuvel1a_SC_vel
        elseif (plate(i).eq.'NNR') then
            geo_pole(i,1) = nuvel1a_NNR_lon
            geo_pole(i,2) = nuvel1a_NNR_lat
            geo_pole(i,3) = nuvel1a_NNR_vel
        else
            write(stderr,*) 'get_pole: no plate named "',trim(plate(i)),'" in NUVEL1A model'
            call list_plate_models('NUVEL1A')
            ierr = 1
            return
        endif
    else
        write(stderr,*) 'get_pole: no model named "',trim(model),'"'
        ierr = 1
        return
    endif
enddo

! Angular velocity of plate2 with respect to plate1: pole2-pole1 = pole
call pole_geo2xyz(geo_pole(1,1),geo_pole(1,2),geo_pole(1,3),xyz_pole(1,1),xyz_pole(1,2),&
                  xyz_pole(1,3),'sphere')
call pole_geo2xyz(geo_pole(2,1),geo_pole(2,2),geo_pole(2,3),xyz_pole(2,1),xyz_pole(2,2),&
                  xyz_pole(2,3),'sphere')
pole(1) = xyz_pole(2,1) - xyz_pole(1,1)
pole(2) = xyz_pole(2,2) - xyz_pole(1,2)
pole(3) = xyz_pole(2,3) - xyz_pole(1,3)

! Convert back to geographic coordinates and rotational velocity
call pole_xyz2geo(pole(1),pole(2),pole(3),geo_pole(1,1),geo_pole(1,2),geo_pole(1,3),'sphere')
pole(1) = geo_pole(1,1)
pole(2) = geo_pole(1,2)
pole(3) = geo_pole(1,3)

return
end subroutine

!--------------------------------------------------------------------------------------------------!

subroutine list_plate_models(model)

use io, only: stdout
implicit none
character(len=*) :: model

if (model.eq.'all') then
    model = 'ALL'
endif

if (model.eq.'ITRF05'.or.model.eq.'ALL') then
    write(stdout,*) 'ITRF 2005 plates (Altamimi et al., 2007)'
    write(stdout,*) '    AM: Amur'
    write(stdout,*) '    AN: Antarctica'
    write(stdout,*) '    AR: Arabia'
    write(stdout,*) '    AU: Australia'
    write(stdout,*) '    CA: Caribbean'
    write(stdout,*) '    EU: Eurasia'
    write(stdout,*) '    IN: India'
    write(stdout,*) '    ITRF05: International Reference Frame 2005'
    write(stdout,*) '    NA: North America'
    write(stdout,*) '    NU: Nubia'
    write(stdout,*) '    NZ: Nazca'
    write(stdout,*) '    OK: Okhotsk'
    write(stdout,*) '    PA: Pacific'
    write(stdout,*) '    SA: South America'
    write(stdout,*) '    SM: Somalia'
    write(stdout,*) '    YZ: Yangtze'
    write(stdout,*)
endif
if (model.eq.'ITRF08'.or.model.eq.'ALL') then
    write(stdout,*) 'ITRF 2008 plates (Altamimi et al., 2012)'
    write(stdout,*) '    AM: Amur'
    write(stdout,*) '    AN: Antarctica'
    write(stdout,*) '    AR: Arabia'
    write(stdout,*) '    AU: Australia'
    write(stdout,*) '    CA: Caribbean'
    write(stdout,*) '    EU: Eurasia'
    write(stdout,*) '    IN: India'
    write(stdout,*) '    ITRF05: International Reference Frame 2005'
    write(stdout,*) '    ITRF08: International Reference Frame 2008'
    write(stdout,*) '    NA: North America'
    write(stdout,*) '    NU: Nubia'
    write(stdout,*) '    NZ: Nazca'
    write(stdout,*) '    PA: Pacific'
    write(stdout,*) '    SA: South America'
    write(stdout,*) '    SM: Somalia'
    write(stdout,*) '    SU: Sunda'
    write(stdout,*)
endif
if (model.eq.'MORVEL'.or.model.eq.'ALL') then
    write(stdout,*) 'MORVEL plates (DeMets et al., 2010)'
    write(stdout,*) '    AM: Amur'
    write(stdout,*) '    AN: Antarctica'
    write(stdout,*) '    AR: Arabia'
    write(stdout,*) '    AU: Australia'
    write(stdout,*) '    CA: Caribbean'
    write(stdout,*) '    CO: Cocos'
    write(stdout,*) '    CP: Capricorn'
    write(stdout,*) '    EU: Eurasia'
    write(stdout,*) '    IN: India'
    write(stdout,*) '    ITRF05: International Reference Frame 2005'
    write(stdout,*) '    JF: Juan de Fuca'
    write(stdout,*) '    LW: Lwandle'
    write(stdout,*) '    MQ: Macquarie'
    write(stdout,*) '    NA: North America'
    write(stdout,*) '    NU: Nubia'
    write(stdout,*) '    NZ: Nazca'
    write(stdout,*) '    PA: Pacific'
    write(stdout,*) '    PS: Philippine Sea'
    write(stdout,*) '    RI: Rivera'
    write(stdout,*) '    SA: South America'
    write(stdout,*) '    SC: Scotia'
    write(stdout,*) '    SM: Somalia'
    write(stdout,*) '    SR: Sur'
    write(stdout,*) '    SU: Sunda'
    write(stdout,*) '    SW: Sandwich'
    write(stdout,*) '    YZ: Yangtze'
    write(stdout,*)
endif
if (model.eq.'MORVEL56'.or.model.eq.'MORVEL-56'.or.model.eq.'ALL') then
    write(stdout,*) 'MORVEL56 plates (Argus et al., 2011)'
    write(stdout,*) '    MORVEL plates (DeMets et al., 2010)'
    write(stdout,*) '    AM: Amur'
    write(stdout,*) '    AN: Antarctica'
    write(stdout,*) '    AR: Arabia'
    write(stdout,*) '    AU: Australia'
    write(stdout,*) '    CA: Caribbean'
    write(stdout,*) '    CO: Cocos'
    write(stdout,*) '    CP: Capricorn'
    write(stdout,*) '    EU: Eurasia'
    write(stdout,*) '    IN: India'
    write(stdout,*) '    JF: Juan de Fuca'
    write(stdout,*) '    LW: Lwandle'
    write(stdout,*) '    MQ: Macquarie'
    write(stdout,*) '    NA: North America'
    write(stdout,*) '    NU: Nubia'
    write(stdout,*) '    NZ: Nazca'
    write(stdout,*) '    PA: Pacific'
    write(stdout,*) '    PS: Philippine Sea'
    write(stdout,*) '    RI: Rivera'
    write(stdout,*) '    SA: South America'
    write(stdout,*) '    SC: Scotia'
    write(stdout,*) '    SM: Somalia'
    write(stdout,*) '    SR: Sur'
    write(stdout,*) '    SU: Sunda'
    write(stdout,*) '    SW: Sandwich'
    write(stdout,*) '    YZ: Yangtze'
    write(stdout,*)
    write(stdout,*) '    Bird (2002) plates'
    write(stdout,*) '    AS: Aegean Sea'
    write(stdout,*) '    AP: Altiplano'
    write(stdout,*) '    AT: Anatolia'
    write(stdout,*) '    BH: Birds Head'
    write(stdout,*) '    BR: Balmoral Reef'
    write(stdout,*) '    BS: Banda Sea'
    write(stdout,*) '    BU: Burma'
    write(stdout,*) '    CL: Caroline'
    write(stdout,*) '    CR: Conway Reef'
    write(stdout,*) '    EA: Easter'
    write(stdout,*) '    FT: Futuna'
    write(stdout,*) '    GP: Galapagos'
    write(stdout,*) '    JZ: Juan Fernandez'
    write(stdout,*) '    KE: Kermadec'
    write(stdout,*) '    MA: Mariana'
    write(stdout,*) '    MN: Manus'
    write(stdout,*) '    MO: Maoke'
    write(stdout,*) '    MS: Molucca Sea'
    write(stdout,*) '    NB: North Bismarck'
    write(stdout,*) '    ND: North Andes'
    write(stdout,*) '    NH: New Hebrides'
    write(stdout,*) '    NI: Niuafoou'
    write(stdout,*) '    OK: Okhotsk'
    write(stdout,*) '    ON: Okinawa'
    write(stdout,*) '    PM: Panama'
    write(stdout,*) '    SB: South Bismarck'
    write(stdout,*) '    SL: Shetland'
    write(stdout,*) '    SS: Solomon Sea'
    write(stdout,*) '    TI: Timor'
    write(stdout,*) '    TO: Tonga'
    write(stdout,*) '    WL: Woodlark'
    write(stdout,*)
    write(stdout,*) '    NNR: No net rotation'
    write(stdout,*) '    ITRF08: International Reference Frame 2008'
    write(stdout,*)
endif
if (model.eq.'NUVEL1A'.or.model.eq.'NUVEL-1A'.or.model.eq.'ALL') then
    write(stdout,*) 'NUVEL-1A plates (DeMets et al., 1994)'
    write(stdout,*) '    AF: Africa'
    write(stdout,*) '    AN: Antarctica'
    write(stdout,*) '    AR: Arabia'
    write(stdout,*) '    AU: Australia'
    write(stdout,*) '    CA: Caribbean'
    write(stdout,*) '    CO: Cocos'
    write(stdout,*) '    EU: Eurasia'
    write(stdout,*) '    IN: India'
    write(stdout,*) '    NNR: No net rotation'
    write(stdout,*) '    JF: Juan de Fuca'
    write(stdout,*) '    NA: North America'
    write(stdout,*) '    NZ: Nazca'
    write(stdout,*) '    PA: Pacific'
    write(stdout,*) '    PS: Philippine Sea'
    write(stdout,*) '    RI: Rivera'
    write(stdout,*) '    SA: South America'
    write(stdout,*) '    SC: Scotia'
    write(stdout,*)
endif


return
end subroutine

end module
