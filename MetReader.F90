      module MetReader

      integer, parameter,private :: sp        = 4 ! single precision
      integer, parameter,private :: dp        = 8 ! double precision

      real(kind=sp), parameter :: RAD_EARTH_MET = 6371.229_sp ! Radius of Earth in km
      real(kind=sp), parameter :: DEG2RAD_MET   = 1.7453292519943295e-2_sp

      integer,public :: MR_iwind       !     MR_IWIND specifies the type of wind input to the model:
                             !   MR_IWIND=1 if a 1-D wind sounding is use, 
                             !           =2 if a 3-D grid is read from a ASCII file.
                             !           =3 if a single, multistep 3-D NetCDF file is used
                             !           =4 if multiple 3-D NetCDF files are used
                             !           =5 if multiple file with multiple steps are used
      integer,public :: MR_iwindformat !      MR_iwindformat specifies the format of the met data
                                         !  0 Custom format based on template
                                         !  1 ASCII profile
                                         !  2 Radiosonde data
                                         !  3 NARR3D NAM221 32 km North America files
                                         !  4 NARR3D NAM221 32 km North America files
                                         !  5 NAM216 AK 45km
                                         !  6 NAM Regional 90 km grid 104
                                         !  7 CONUS 212 40km
                                         !  8 CONUS 218 (12km)
                                         !  9 Unassigned
                                         ! 10 NAM 242 11.25 km AK
                                         ! 11 NAM 196 2.5 km HI
                                         ! 12 NAM 198 5.953 km AK
                                         ! 13 NAM 91 2.976 km AK
                                         ! 20 GFS 0.5
                                         ! 21 Old format GFS 0.5-degree
                                         ! 22 GFS 0.25
                                         ! 23 NCEP / DOE reanalysis 2.5 degree files
                                         ! 24 NASA-MERRA reanalysis 1.25 degree files
                                         ! 25 NCEP/NCAR reanalysis 2.5 degree files
                                         ! 27 NOAA-CIRES reanalysis 2.5 degree files
                                         ! 28 ECMWF Interim Reanalysis (ERA-Interim)
                                         ! 31 Catania forecast
                                         ! 32 Air Force Weather Agency subcenter = 0
                                         ! 33 CCSM3.0 Community Atmosphere Model (CAM)
                                         ! 40 NASA-GEOS Cp
                                         ! 41 NASA-GEOS Np
                                         ! 50 WRF - output

      integer,public :: MR_iGridCode   !   MR_iGridCode specifies the NCEP grid described in:
                                    !   http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
      integer,public :: MR_idataFormat !   Specifies the data model used
                                    !    =1 ASCII
                                    !    =2 netcdf
                                    !    =3 grib
                                    !    =4 hdf

        ! These variables describe the full list of windfiles read
      integer,public                                        :: MR_iwindfiles           ! number of files provided
      character(len=130), allocatable,dimension(:)  ,public :: MR_windfiles            ! name of file
      real(kind=dp)     , allocatable,dimension(:)  ,public :: MR_windfile_starthour   ! start hour of the file
      integer           , allocatable,dimension(:)  ,public :: MR_windfiles_nt_fullmet ! number of steps in files
      real(kind=dp)     , allocatable,dimension(:,:),public :: MR_windfile_stephour    ! offset hours of step
      character(len=80),public                              :: MR_iwf_template         ! name of the template file
      logical           , allocatable,dimension(:)  ,public :: MR_windfiles_Have_GRIB_index
      character(len=130), allocatable,dimension(:)  ,public :: MR_windfiles_GRIB_index ! name of grib index file

        ! These variables are a list of the same data above, but specific to the simulation
        ! duration
      integer,public                                       :: MR_MetSteps_Total
      integer,public                                       :: MR_iMetStep_Now
      character(len=130), allocatable,dimension(:),public  :: MR_MetStep_File
      integer           , allocatable,dimension(:),public  :: MR_MetStep_findex
      integer           , allocatable,dimension(:),public  :: MR_MetStep_tindex
      real(kind=dp)     , allocatable,dimension(:),public  :: MR_MetStep_Hour_since_baseyear
      real(kind=dp)     , allocatable,dimension(:),public  :: MR_MetStep_Interval
      integer           , allocatable,dimension(:),public  :: MR_MetStep_year
      integer           , allocatable,dimension(:),public  :: MR_MetStep_month
      integer           , allocatable,dimension(:),public  :: MR_MetStep_day
      integer           , allocatable,dimension(:),public  :: MR_MetStep_DOY
      real(kind=dp)     , allocatable,dimension(:),public  :: MR_MetStep_Hour_Of_Day

      real(kind=dp),public                                 :: MR_ForecastInterval
      logical                                              :: MR_runAsForecast = .false.
      real(kind=dp)                                        :: MR_FC_Offset     = 0.0_dp
      logical                                              :: MR_Reannalysis   = .false.
      integer           ,dimension(:), allocatable,public  :: MR_iwind5_year

      !    Native grid of Met file using Pressure as vertical coordinate
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_met_int   ! Used for catagorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_met       ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d2_metP
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_last  ! These are needed for compH interpolation
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_geoH_metP_next

      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_last  !These might need to be stored to avoid a 
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vx_metP_next  !second reading
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_last
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_vy_metP_next
      logical, public :: MR_Save_Velocities = .false.

      real(kind=sp), public :: MR_Max_geoH_metP_predicted   ! Highest expected height in Met file, based on pressure
      real(kind=sp) :: Max_geoH_metP_last
      real(kind=sp) :: Max_geoH_metP_next
      real(kind=sp) :: Max_geoH_metP
      real(kind=sp) :: Suppl_H
      !Parameter iHeightHandler specifies what to do if the maximum height
      !of the simulation region exceeds the maximum height in the wind files.
      !If iHeightHandler = 1, stop the program if the plume height exceeds mesoscale height
      !                    2, wind direction at levels above the highest node 
      !                       equal that of the highest node.  Temperatures in the
      !                       upper nodes don't change between 11 and 20 km; above
      !                       20 km they increase by 2 C/km, as in the Standard
      !                       atmosphere.  A warning is written to the log file.
      integer       :: MR_iHeightHandler = 2

        ! The following are used for sonde data
      integer :: MR_nSnd_Locs      = 1  ! Number of Sonde locations
      integer :: MR_Snd_nt_fullmet = 1  ! Number of times at the Sonde locations
      integer :: MR_Snd_nvars      = 5  ! Number of Sonde variables (P,H,U,V,T,+user-specified)
      logical :: Snd_Have_PT = .false.
      real(kind=sp),dimension(:,:,:,:),allocatable, public :: MR_SndVars_metP   ! (MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,300)
      integer      ,dimension(:,:)    ,allocatable, public :: MR_Snd_np_fullmet ! Number of pressure values for each location/time
      integer      ,dimension(:)      ,allocatable, public :: MR_SndVarsID      ! Lists which vars are in which columns of  MR_SndVars_metP
      real(kind=sp),dimension(:,:,:)  ,allocatable, public :: MR_Snd2Comp_tri_map_wgt ! weights of nearby sondes for every comp point
      integer      ,dimension(:,:,:)  ,allocatable, public :: MR_Snd2Comp_tri_map_idx ! sonde index of weights

      !    Native grid of Met file using Height as vertical coordinate
      !    (resampled onto z-gridpoints of computational grid)
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_metH
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_u_ER_metP ! For the cases where Met is proj and comp
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_v_ER_metP !  different we need to rotate so these
                                                                         !  store Earth-Relative velocities on MetP

      !    Computations grid
      integer      ,dimension(:,:)  ,allocatable, public :: MR_dum2d_comp_int  ! Used for catagorical variables
      real(kind=sp),dimension(:,:)  ,allocatable, public :: MR_dum2d_comp      ! Used for surface variables
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH
      real(kind=sp),dimension(:,:,:),allocatable, public :: MR_dum3d_compH_2   ! Used only when a vector field
                                                                               ! rotation is needed

      integer       :: nt_fullmet    ! length of t of file
      integer       :: nx_fullmet    ! length of x or lon of full met grid
      integer       :: ny_fullmet    ! length of y or lat of full met grid
      integer       :: np_fullmet    ! length of pressure of full met grid
      integer       :: np_fullmet_Vz ! length of pressure of full met grid for Vz
      integer       :: np_fullmet_RH ! length of pressure of full met grid for RH
      integer       :: np_fullmet_P0 ! length of pressure of full met grid for Precip 
                                     !   (=1 for all cases except iwindformat=24)
      integer       :: np_fullmet_pad = 1 ! We might need to pad the top of the pressure grid.
      integer       :: neta_fullmet  ! Only used by WRF

      real(kind=sp),dimension(:), allocatable :: x_fullmet_sp    ! x-coordinates of full met grid
      real(kind=sp),dimension(:), allocatable :: y_fullmet_sp    ! y-coordinates of full met grid
      real(kind=sp),dimension(:), allocatable :: p_fullmet_sp    ! z-coordinates of full met grid for U,V,T,H
      real(kind=sp),dimension(:), allocatable :: p_fullmet_Vz_sp ! z-coordinates of full met grid for Vz
      real(kind=sp),dimension(:), allocatable :: p_fullmet_RH_sp ! z-coordinates of full met grid for RH

      logical       :: IsLatLon_MetGrid
      logical       :: IsGlobal_MetGrid     = .false.  ! Not all Lon/Lat grids are periodic
      logical       :: IsLatLon_CompGrid
      logical       :: IsPeriodic_CompGrid  = .false.
      logical       :: UseFullMetGrid       = .false.  ! This is the special case where the comp grid
                                                       ! equals the Met grid

      integer       :: nx_submet ! length of x or lon of sub-grid
      integer       :: ny_submet ! length of y or lat of sub-grid
      real(kind=sp),dimension(:), allocatable :: x_submet_sp ! x-coordinates of met sub-grid
      real(kind=sp),dimension(:), allocatable :: y_submet_sp ! y-coordinates of met sub-grid
      real(kind=sp),dimension(:), allocatable :: z_approx ! zpproximate altidue from STD Atmos and press (in km)

      logical :: x_inverted     = .false.
      logical :: y_inverted     = .false. ! Some LatLon grids start at the North Pole and increment down
      logical :: z_inverted     = .false. ! Some grids give top pressure first
      logical :: y_pad_North    = .false. ! Some computational grids will require values above the top met point
      logical :: y_pad_South    = .false. !   
      logical :: isGridRelative = .true.  ! Most windfiles, whether Lat/Lon or projected, give
                                          ! velocities relative to the grid of the file.  Some (NARR)
                                          ! give velocities relative to earth coordinates and need to
                                          ! be rotated
      real(kind=dp),dimension(:,:) ,allocatable :: theta_Met  ! Earth to grid
      real(kind=dp),dimension(:,:) ,allocatable :: theta_Comp ! Earth to grid

      ! Met copies of projection variables, used for proj call on Met Grid
      integer      :: Met_iprojflag
      real(kind=8) :: Met_Re
      real(kind=8) :: Met_k0
      real(kind=8) :: Met_phi0    != 90.0_ip        ! latitude of projection point
      real(kind=8) :: Met_lam0 != -135.0_ip   ! longitude of projection point
      real(kind=8) :: Met_lam1,Met_phi1
      real(kind=8) :: Met_lam2,Met_phi2

      integer      :: Comp_iprojflag
      real(kind=8) :: Comp_Re
      real(kind=8) :: Comp_k0
      real(kind=8) :: Comp_phi0    != 90.0_ip        ! latitude of projection point
      real(kind=8) :: Comp_lam0 != -135.0_ip   ! longitude of projection point
      real(kind=8) :: Comp_lam1,Comp_phi1
      real(kind=8) :: Comp_lam2,Comp_phi2

      integer      :: Map_Case

      ! Some geometry terms
      real(kind=sp),dimension(:,:)     ,allocatable,private :: rdphi_MetP_sp
      real(kind=sp),dimension(:,:,:)   ,allocatable,private :: rdlambda_MetP_sp
      real(kind=sp),dimension(:)       ,allocatable         :: MR_dx_met, MR_dx_submet
      real(kind=sp),dimension(:)       ,allocatable         :: MR_dy_met, MR_dy_submet

      logical       :: IsRegular_MetGrid
      real(kind=sp) :: dx_met_const, dy_met_const

        ! There are some Ash3d variables we might need, so make local copies
      integer       :: MR_BaseYear            = 1900       ! This should be reset in calling program
      logical       :: MR_useLeap             = .true.  ! This too
      integer       :: MR_Comp_StartYear
      real(kind=dp) :: MR_Comp_StartHour
      real(kind=dp) :: MR_Comp_Time_in_hours
      integer       :: nx_comp
      integer       :: ny_comp
      integer       :: nz_comp
      real(kind=sp),dimension(:), allocatable :: x_comp_sp ! x-coordinates of computational grid
      real(kind=sp),dimension(:), allocatable :: y_comp_sp ! y-coordinates of computational grid
      real(kind=sp),dimension(:), allocatable :: z_comp_sp ! z-coordinates of computational grid
      real(kind=sp) :: dx_comp,dy_comp
      real(kind=sp) :: MaxZ_comp_sp

      real(kind=sp),dimension(:,:), allocatable :: CompPoint_X_on_Met_sp   ! x-coord (on Met grid) of comp point
      real(kind=sp),dimension(:,:), allocatable :: CompPoint_Y_on_Met_sp   ! y-coord (on Met grid) of comp point
      integer,dimension(:,:,:), allocatable     :: CompPoint_on_subMet_idx ! index on met sub-grid of comp point
      real(kind=sp),dimension(:,:,:), allocatable :: bilin_map_wgt
      logical :: Have_Vz = .false.
      real(kind=sp),dimension(0:50)  :: fill_value_sp = -9999.0_sp
      character(len=30),dimension(9) :: Met_dim_names
      logical          ,dimension(9) :: Met_dim_IsAvailable
      real(kind=sp)    ,dimension(9) :: Met_dim_fac  = 1.0_sp
      ! Here is the list of variables that can be read.  Each iwindformat will
      ! have just a sub-set availible with specific names.  For now, allocate
      ! space for 50 variable names
        ! Mechanical / State variables
        !   1 = Geopotential Height
        !   2 = Vx
        !   3 = Vy
        !   4 = Vz
        !   5 = Temperature
        ! Surface
        !  10 = Planetary Boundary Layer Height
        !  11 = U @ 10m
        !  12 = V @ 10m
        !  13 = Friction velocity
        !  14 = Displacement Height
        !  15 = Snow cover
        !  16 = Soil moisture
        !  17 = Surface Roughness
        !  18 = Wind gust speed
        !  19 = surface temperature
        ! Atmospheric Structure
        !  20 = pressure at lower cloud base
        !  21 = pressure at lower cloud top
        !  22 = temperature at lower cloud top
        !  23 = Total Cloud cover
        !  24 = Cloud cover (low)
        !  25 = Cloud cover (convective)
        ! Moisture
        !  30 = Rel. Hum
        !  31 = QV (specific humidity)
        !  32 = QL (liquid)
        !  33 = QI (ice)
        ! Precipitation
        !  40 = Categorical rain
        !  41 = Categorical snow
        !  42 = Categorical frozen rain
        !  43 = Categorical ice
        !  44 = Precipitation rate large-scale (liquid)
        !  45 = Precipitation rate convective  (liquid)
        !  46 = Precipitation rate large-scale (ice)
        !  47 = Precipitation rate convective  (ice)
      logical          ,dimension(50)   :: Met_var_IsAvailable
      character(len=71),dimension(50)   :: Met_var_names
      character(len=5) ,dimension(50)   :: Met_var_names_WMO
      integer          ,dimension(50)   :: Met_var_ndim
      integer          ,dimension(50)   :: Met_var_zdimID
      integer          ,dimension(50,4) :: Met_var_GRIB2_DPcPnSt

      !logical          ,dimension(50) :: Met_var_IsFloat  ! true if kind=4 otherwise false
      real(kind=sp)    ,dimension(50) :: Met_var_conversion_factor

        ! Variables needed by netcdf reader
      real(kind=sp) :: iwf25_scale_facs(50)
      real(kind=sp) :: iwf25_offsets(50)
      real(kind=sp) :: x_in_iwf25_sp(192)
      real(kind=sp) :: y_in_iwf25_sp(94)
        ! Here is the mapping for bilinear weighting coeffiecients (amap) and
        ! indices (imap) from the 1.875-deg 2d grid to the 2.5-deg 
      real(kind=sp),dimension(:,:,:)   ,allocatable :: amap_iwf25
      integer      ,dimension(:,:,:)   ,allocatable :: imap_iwf25
      integer(kind=sp),dimension(:,:,:),allocatable :: tmpsurf2d_short

      integer :: istart, iend
      integer :: jstart, jend
      !integer :: ilhalf,irhalf
      integer :: ilhalf_fm_l, ilhalf_fm_r
      integer :: irhalf_fm_l, irhalf_fm_r
      integer :: ilhalf_nx, irhalf_nx
      logical :: wrapgrid

      real(kind=sp)   ,dimension(:)       ,allocatable :: temp1d_sp
      real(kind=sp)   ,dimension(:,:,:)   ,allocatable :: temp2d_sp
      real(kind=sp)   ,dimension(:,:,:,:) ,allocatable :: temp3d_sp
      integer(kind=sp),dimension(:,:,:)   ,allocatable :: temp2d_int
      integer(kind=sp),dimension(:,:,:)   ,allocatable :: temp2d_short
      integer(kind=sp),dimension(:,:,:,:) ,allocatable :: temp3d_short

      real(kind=4),dimension(:,:),allocatable :: Met_Proj_lat
      real(kind=4),dimension(:,:),allocatable :: Met_Proj_lon

      contains

!##############################################################################
!
!     MR_Reset_Memory
!
!     This subroutine reinitializes MetReader by deallocating all MR variables.
!     This is useful if a program needs to use multiple types of wind files.
!
!##############################################################################

      subroutine MR_Reset_Memory

      implicit none

       write(*,*)"-------------------------------------------------------"
       write(*,*)"-------- Resetting all MetReader Memory ---------------"
       write(*,*)"-------------------------------------------------------"

       IF(allocated(MR_windfiles                  ))deallocate(MR_windfiles)
       IF(allocated(MR_windfile_starthour         ))deallocate(MR_windfile_starthour)
       IF(allocated(MR_windfiles_nt_fullmet       ))deallocate(MR_windfiles_nt_fullmet)
       IF(allocated(MR_windfile_stephour          ))deallocate(MR_windfile_stephour)
       IF(allocated(MR_windfiles_Have_GRIB_index  ))deallocate(MR_windfiles_Have_GRIB_index)
       IF(allocated(MR_windfiles_GRIB_index       ))deallocate(MR_windfiles_GRIB_index)
       IF(allocated(MR_MetStep_File               ))deallocate(MR_MetStep_File)
       IF(allocated(MR_MetStep_findex             ))deallocate(MR_MetStep_findex)
       IF(allocated(MR_MetStep_tindex             ))deallocate(MR_MetStep_tindex)
       IF(allocated(MR_MetStep_Hour_since_baseyear))deallocate(MR_MetStep_Hour_since_baseyear)
       IF(allocated(MR_MetStep_Interval           ))deallocate(MR_MetStep_Interval)
       IF(allocated(MR_MetStep_year               ))deallocate(MR_MetStep_year)
       IF(allocated(MR_MetStep_month              ))deallocate(MR_MetStep_month)
       IF(allocated(MR_MetStep_day                ))deallocate(MR_MetStep_day)
       IF(allocated(MR_MetStep_DOY                ))deallocate(MR_MetStep_DOY)
       IF(allocated(MR_MetStep_Hour_Of_Day        ))deallocate(MR_MetStep_Hour_Of_Day)
       IF(allocated(MR_iwind5_year                ))deallocate(MR_iwind5_year)
       IF(allocated(MR_dum2d_met_int              ))deallocate(MR_dum2d_met_int)
       IF(allocated(MR_dum2d_met                  ))deallocate(MR_dum2d_met)
       IF(allocated(MR_dum3d_metP                 ))deallocate(MR_dum3d_metP)
       IF(allocated(MR_dum3d2_metP                ))deallocate(MR_dum3d2_metP)
       IF(allocated(MR_geoH_metP_last             ))deallocate(MR_geoH_metP_last)
       IF(allocated(MR_geoH_metP_next             ))deallocate(MR_geoH_metP_next)
       IF(allocated(MR_vx_metP_last               ))deallocate(MR_vx_metP_last)
       IF(allocated(MR_vx_metP_next               ))deallocate(MR_vx_metP_next)
       IF(allocated(MR_vy_metP_last               ))deallocate(MR_vy_metP_last)
       IF(allocated(MR_vy_metP_next               ))deallocate(MR_vy_metP_next)
       IF(allocated(MR_dum3d_metH                 ))deallocate(MR_dum3d_metH)
       IF(allocated(MR_u_ER_metP                  ))deallocate(MR_u_ER_metP)
       IF(allocated(MR_v_ER_metP                  ))deallocate(MR_v_ER_metP)
       IF(allocated(MR_dum2d_comp_int             ))deallocate(MR_dum2d_comp_int)
       IF(allocated(MR_dum2d_comp                 ))deallocate(MR_dum2d_comp)
       IF(allocated(MR_dum3d_compH                ))deallocate(MR_dum3d_compH)
       IF(allocated(MR_dum3d_compH_2              ))deallocate(MR_dum3d_compH_2)
       IF(allocated(x_fullmet_sp                  ))deallocate(x_fullmet_sp)
       IF(allocated(y_fullmet_sp                  ))deallocate(y_fullmet_sp)
       IF(allocated(p_fullmet_sp                  ))deallocate(p_fullmet_sp)
       IF(allocated(p_fullmet_Vz_sp               ))deallocate(p_fullmet_Vz_sp)
       IF(allocated(p_fullmet_RH_sp               ))deallocate(p_fullmet_RH_sp)
       IF(allocated(x_submet_sp                   ))deallocate(x_submet_sp)
       IF(allocated(y_submet_sp                   ))deallocate(y_submet_sp)
       IF(allocated(z_approx                      ))deallocate(z_approx)
       IF(allocated(theta_Met                     ))deallocate(theta_Met)
       IF(allocated(theta_Comp                    ))deallocate(theta_Comp)
       IF(allocated(rdphi_MetP_sp                 ))deallocate(rdphi_MetP_sp)
       IF(allocated(rdlambda_MetP_sp              ))deallocate(rdlambda_MetP_sp)
       IF(allocated(MR_dx_met                     ))deallocate(MR_dx_met)
       IF(allocated(MR_dx_submet                  ))deallocate(MR_dx_submet)
       IF(allocated(MR_dy_met                     ))deallocate(MR_dy_met)
       IF(allocated(MR_dy_submet                  ))deallocate(MR_dy_submet)
       IF(allocated(x_comp_sp                     ))deallocate(x_comp_sp)
       IF(allocated(y_comp_sp                     ))deallocate(y_comp_sp)
       IF(allocated(z_comp_sp                     ))deallocate(z_comp_sp)
       IF(allocated(CompPoint_X_on_Met_sp         ))deallocate(CompPoint_X_on_Met_sp)
       IF(allocated(CompPoint_Y_on_Met_sp         ))deallocate(CompPoint_Y_on_Met_sp)
       IF(allocated(CompPoint_on_subMet_idx       ))deallocate(CompPoint_on_subMet_idx)
       IF(allocated(bilin_map_wgt                 ))deallocate(bilin_map_wgt)
       IF(allocated(amap_iwf25                    ))deallocate(amap_iwf25)
       IF(allocated(imap_iwf25                    ))deallocate(imap_iwf25)
       IF(allocated(tmpsurf2d_short               ))deallocate(tmpsurf2d_short)
       IF(allocated(temp1d_sp                     ))deallocate(temp1d_sp)
       IF(allocated(temp2d_sp                     ))deallocate(temp2d_sp)
       IF(allocated(temp3d_sp                     ))deallocate(temp3d_sp)
       IF(allocated(temp2d_int                    ))deallocate(temp2d_int)
       IF(allocated(temp2d_short                  ))deallocate(temp2d_short)
       IF(allocated(temp3d_short                  ))deallocate(temp3d_short)
       IF(allocated(Met_Proj_lat                  ))deallocate(Met_Proj_lat)
       IF(allocated(Met_Proj_lon                  ))deallocate(Met_Proj_lon)

      end subroutine MR_Reset_Memory

!##############################################################################
!
!     MR_Allocate_FullMetFileList
!
!     This subroutine allocates the list of windfiles and does some
!     error-checking based on iwind and iwindfiles.
!
!     From the calling program, this is called once the information about
!     the NWP files is available (e.g. the iwind, iwindformat, grid ID, data
!     format, and number of windfiles.
!
!     Variables allocated:
!        MR_windfiles(MR_iwindfiles)
!        MR_windfiles_nt_fullmet(MR_iwindfiles)
!       and possibly:
!        MR_windfiles_Have_GRIB_index(MR_iwindfiles)
!        MR_windfiles_GRIB_index(MR_iwindfiles)
!
!     The next step is for the calling program to populate MR_windfiles
!
!##############################################################################

      subroutine MR_Allocate_FullMetFileList(iw,iwf,igrid,idf,iwfiles)

      implicit none

      integer,intent(in)      :: iw
      integer,intent(in)      :: iwf
      integer,intent(in)      :: igrid
      integer,intent(in)      :: idf
      integer,intent(in)      :: iwfiles

      MR_iwind              = iw
      MR_iwindformat        = iwf
      IF(igrid.eq.0)THEN
        MR_iGridCode        = 1000+iwf ! Some grids are not in the NCEP list.  For those, reset igrid to
                                       ! the iwindformat code + 1000 
      ELSE
        MR_iGridCode          = igrid
      ENDIF
      MR_idataFormat        = idf
      MR_iwindfiles         = iwfiles

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------      MR_Allocate_FullMetFileList                           ----------"
      write(*,*)"--------------------------------------------------------------------------------"

      if ((MR_iwind.ne.1).and.(MR_iwind.ne.2).and. &
          (MR_iwind.ne.3).and.(MR_iwind.ne.4).and. &
          (MR_iwind.ne.5)) then
        write(6,*)'MR_iwind must be between 1 and 5. Program stopped'
        write(6,*)' MR_iwind = ',MR_iwind
        write(6,*)' MR_IWIND OPTIONS:'
        write(6,*)' MR_iwind = 1 read from a 1-D wind sounding'
        write(6,*)'            2 read from 3D gridded ASCII files'
        write(6,*)'            3 read from a single, multistep file'
        write(6,*)'            4 read from multiple files'
        write(6,*)'            5 read variables from separate files'
        stop 1
      endif

      if     (MR_iwindformat.eq.0) THEN  ! Custom format based on template
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Custom format based on template"
      elseif (MR_iwindformat.eq.1) THEN  ! ASCII profile
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ASCII profile"
      elseif (MR_iwindformat.eq.2) THEN  ! Radiosonde data
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Radiosonde data"
      elseif (MR_iwindformat.eq.3) THEN  ! NARR3D NAM221 32 km North America files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NARR3D NAM221 32 km North America files"
      elseif (MR_iwindformat.eq.4) THEN  ! NARR3D NAM221 32 km North America files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NARR3D NAM221 32 km North America files"
      elseif (MR_iwindformat.eq.5) THEN  ! NAM216 AK 45km
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM216 AK 45km"
      elseif (MR_iwindformat.eq.6) THEN  ! NAM Regional 90 km grid 104
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM Regional 90 km grid 104"
      elseif (MR_iwindformat.eq.7) THEN  ! CONUS 212 40km
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS 212 40km"
      elseif (MR_iwindformat.eq.8) THEN  ! CONUS 218 (12km)
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CONUS 218 (12km)"
      elseif (MR_iwindformat.eq.9) THEN  ! Unassigned
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Unassigned"
      elseif (MR_iwindformat.eq.10)THEN  ! NAM 242 11.25 km AK
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 242 11.25 km AK"
      elseif (MR_iwindformat.eq.11)THEN  ! NAM 196 2.5 km HI
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  " NAM 196 2.5 km HI"
      elseif (MR_iwindformat.eq.12)THEN  ! NAM 198 5.953 km AK
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 198 5.953 km AK"
      elseif (MR_iwindformat.eq.13)THEN  ! NAM 91 2.976 km AK
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NAM 91 2.976 km AK"
      elseif (MR_iwindformat.eq.20)THEN  ! GFS 0.5
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.5"
      elseif (MR_iwindformat.eq.21)THEN  ! Old format GFS 0.5-degree
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.5-degree"
      elseif (MR_iwindformat.eq.22)THEN  ! GFS 0.25
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.25"
      elseif (MR_iwindformat.eq.23)THEN  ! NCEP / DOE reanalysis 2.5 degree files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NCEP / DOE reanalysis 2.5 degree files"
      elseif (MR_iwindformat.eq.24)THEN  ! NASA-MERRA reanalysis 1.25 degree files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-MERRA reanalysis 1.25 degree files"
      elseif (MR_iwindformat.eq.25)THEN  ! NCEP/NCAR reanalysis 2.5 degree files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NCEP/NCAR reanalysis 2.5 degree files"
      elseif (MR_iwindformat.eq.26)THEN  ! GFS 0.5 with 47 pressure levels
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "GFS 0.5 with 47 pressure levels"
      elseif (MR_iwindformat.eq.27)THEN  ! NOAA-CIRES reanalysis 2.5 degree files
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NOAA-CIRES reanalysis 2.5 degree files"
      elseif (MR_iwindformat.eq.28)THEN  ! ECMWF Interim Reanalysis (ERA-Interim)
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "ECMWF Interim Reanalysis (ERA-Interim)"
      elseif (MR_iwindformat.eq.31)THEN  ! Catania forecast
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Catania forecast"
      elseif (MR_iwindformat.eq.32)THEN  ! Air Force Weather Agency subcenter = 0
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "Air Force Weather Agency subcenter = 0"
      elseif (MR_iwindformat.eq.33)THEN  ! CCSM3.0 Community Atmosphere Model (CAM)
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "CCSM3.0 Community Atmosphere Model (CAM)"
      elseif (MR_iwindformat.eq.40)THEN  ! NASA-GEOS Cp
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-GEOS Cp"
      elseif (MR_iwindformat.eq.41)THEN  ! NASA-GEOS Np
        MR_Reannalysis = .false.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "NASA-GEOS Np"
      elseif (MR_iwindformat.eq.50)THEN   ! WRF - output
        MR_Reannalysis = .true.
        write(*,*)"  NWP format to be used = ",MR_iwindformat,&
                  "WRF"
      else
        write(6,*)'MR_iwindformat not supported'
        write(6,*)'MR_iwindformat=',MR_iwindformat,'. Program stopped in MetReader.f90'
        stop 1
      endif
      write(*,*)"                grid ID = ",MR_iGridCode
      select case (MR_idataFormat)
      case(1)
        write(*,*)"            data format = ",MR_idataFormat,"ASCII"
      case(2)
        write(*,*)"            data format = ",MR_idataFormat,"NETCDF"
#ifndef USENETCDF
        write(*,*)"MR ERROR: Met files are netcdf, but MetReader was not"
        write(*,*)"          compiled with netcdf support."
        write(*,*)"          Please recompile MetReader with netcdf support"
        stop 1
#endif
      case(3)
        write(*,*)"            data format = ",MR_idataFormat,"GRIB"
#ifndef USEGRIB2
        write(*,*)"MR ERROR: Met files are grib2, but MetReader was not"
        write(*,*)"          compiled with grib2 support."
        write(*,*)"          Please recompile MetReader with grib2 support"
        stop 1
#endif
      case(4)
        write(*,*)"            data format = ",MR_idataFormat,"HDF"
#ifndef USEHDF
        write(*,*)"MR ERROR: Met files are hdf, but MetReader was not"
        write(*,*)"          compiled with hdf support."
        write(*,*)"          Please recompile MetReader with hdf support"
        stop 1
#endif
      case DEFAULT
        write(*,*)"MR ERROR:  MR_idataFormat not 1:4."
        write(*,*)"           Exiting in MR_Allocate_FullMetFileList"
        stop 1
      end select

      write(*,*)"     Allocating space for ",MR_iwindfiles,"files."

      if (MR_iwindformat.eq.25.or.MR_iwindformat.eq.27)then
        write(*,*)"For NCEP1 or NOAA, only the directory should be listed."
        ! Reset MR_iwindfiles to 2: only one "file" will be read,
        ! but we need this to be 2 to accommodate runs that might span two years
        MR_iwindfiles = 2
      endif

      allocate (MR_windfiles(MR_iwindfiles))
      allocate (MR_windfiles_nt_fullmet(MR_iwindfiles))
      IF(MR_idataFormat.eq.3)THEN
        allocate (MR_windfiles_Have_GRIB_index(MR_iwindfiles))
          ! This will be reset to true if the index files are found
        MR_windfiles_Have_GRIB_index = .false.
        allocate (MR_windfiles_GRIB_index(MR_iwindfiles))
      ENDIF

      If(MR_iwind.eq.1)THEN
        ! For the 1d profile or radiosonde case, igrid is used for the number of sonde
        ! locations.  If it is not provided, set it to one
        IF(igrid.eq.0)THEN
          MR_iGridCode = 1
          MR_nSnd_Locs = MR_iGridCode
        ELSE
          MR_iGridCode = igrid
          MR_nSnd_Locs = MR_iGridCode
        ENDIF
          ! Now make sure that the number of windfiles is a multiple of the number of locations
          ! since this will be the number of timesteps
        IF(mod(MR_iwindfiles,MR_nSnd_Locs).eq.0)THEN
          MR_Snd_nt_fullmet = MR_iwindfiles / MR_nSnd_Locs
        ELSE
          write(*,*)"MR ERROR:  The grid code for 1d ascii input is interpreted to be"
          write(*,*)"           the number of sonde locations.  Each group of sondes can"
          write(*,*)"           be repeated, correspoding to multiple timesteps.  The"
          write(*,*)"           number of windfiles must be a multiple of the number of"
          write(*,*)"           locations"
          write(*,*)"                   MR_iwind = ",MR_iwind
          write(*,*)"             MR_iwindformat = ",MR_iwindformat
          write(*,*)"               MR_iGridCode = ",MR_iGridCode
          write(*,*)"             MR_idataFormat = ",MR_idataFormat
          write(*,*)"               MR_nSnd_Locs = ",MR_iGridCode
          write(*,*)"          MR_Snd_nt_fullmet = ",real(MR_iwindfiles)/real(MR_nSnd_Locs)
          stop 1
        ENDIF
        IF(MR_nSnd_Locs.gt.1)THEN
          write(*,*)"MR ERROR: Currently, only one sonde location is supported."
          stop 1
        ENDIF
        IF(MR_Snd_nt_fullmet.gt.1)THEN
          write(*,*)"MR ERROR: Currently, only one sonde time is supported."
          stop 1
        ENDIF
      ENDIF
      write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Allocate_FullMetFileList


!##############################################################################
!
!     MR_Read_Met_DimVars
!
!     This subroutine expects that MR_windfiles has been filled by the calling
!     program, and checks for their existance.
!
!     Depending on the NWP data format, a subroutine is called to read the spatial
!     structure of the NWP files, followed by one for the temporal span
!      e.g. MR_Read_Met_DimVars_netcdf
!           MR_Read_Met_Times_netcdf
!
!     From the calling program, this is called once the names of the NWP files
!     is specified.  If a custom netcdf template file is to be used (iwf=0), then
!     the variable MR_iwf_template must also be filled so that it can be read
!     from MR_Read_Met_DimVars_[].  
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!     The next step is for the calling program to specify the projection parameters
!     of the computational grid (i.e. the grid MetReader should be returing values to)
!
!##############################################################################

      subroutine MR_Read_Met_DimVars(iy)

      implicit none

      integer, optional,intent(in) :: iy  ! Note: this is only needed for iwf=25 or 27
                                          !       since we need to know how many
                                          !       metsteps to allocate

      integer      :: i
      logical      :: IsThere

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------      MR_Read_Met_DimVars                                   ----------"
      write(*,*)"--------------------------------------------------------------------------------"


      IF(.not.allocated(MR_windfiles))THEN
        write(*,*)"MR ERROR:  The list of windfile names, MR_windfiles, has not been"
        write(*,*)"           allocated.  The calling program must allocate this"
        write(*,*)"           array via the subroutine MR_Allocate_FullMetFileList."
        stop 1
      ENDIF

      IF(MR_iwindformat.eq.25.or.MR_iwindformat.eq.27)THEN
        ! For NCEP 2.5 degree reanalysis and the NOAA product, only the directory
        ! was read into slot 1 of MR_windfiles(:).  We need to copy to slot 2
        ! to make sure we don't throw an error
        MR_windfiles(2)   = MR_windfiles(1)
        if(present(iy)) then
          ! This is needed at this point for allocating the number
          ! of steps per file (this depends on the year), but this
          ! is reset in MR_Set_Met_Times when the actual start time
          ! is given
          MR_Comp_StartYear = iy
        else
          write(*,*)"MR WARNING: If the iwf=25 (NCEP Reannalysis) or iwf=27 (NOAA"
          write(*,*)"            Reannalysis) are used, then MR_Read_Met_DimVars"
          write(*,*)"            should be called with a start year.  This is needed"
          write(*,*)"            to allocate the correct number of time steps per file."
          write(*,*)"            Setting MR_Comp_StartYear to 2017 for a non-leap year."
          write(*,*)"            However, MR_Comp_StartYear will be checked later to"
          write(*,*)"            verify that the start year is not a leap year."
          write(*,*)"            If there is an inconsistancy, the program will stop."
          write(*,*)"            If MR_Comp_StartYear is changed to a leap year outside"
          write(*,*)"            of MetReader, then the results will be incorrect."
          MR_Comp_StartYear = 2017
        endif

      ENDIF

      ! Check the existance of the wind files
      write(*,*)"  Verifying existance of windfiles:"
      do i=1,MR_iwindfiles
        INQUIRE( FILE=adjustl(trim(MR_windfiles(i))), EXIST=IsThere )
        write(*,*)"     ",adjustl(trim(MR_windfiles(i))),IsThere
        IF(.not.IsThere)THEN
          write(*,*)"MR ERROR: Could not find windfile ",i
          write(*,*)"          adjustl(trim(MR_windfiles(i)))"
          stop 1
        ENDIF
      enddo

      ! Now set up the full spatial and temporal grids
      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_Met_DimVars_ASCII_1d
      case(2)
        !call Set_Read_Met_DimVars_ASCII_3d
      case (3:5)
         ! call routine to populate variable lists and dimension axies based on
         ! MR_iwindformat; also do a bit of error checking by reading netcdf files
        IF(MR_idataFormat.eq.2)THEN
#ifdef USENETCDF
          call MR_Read_Met_DimVars_netcdf
          call MR_Read_Met_Times_netcdf
#endif
        ELSEIF(MR_idataFormat.eq.3)THEN
#ifdef USEGRIB2
          call MR_Read_Met_DimVars_GRIB
          call MR_Read_Met_Times_GRIB
#endif
        ELSEIF(MR_idataFormat.eq.4)THEN
#ifdef USEHDF
          !call MR_Read_Met_DimVars_HDF
#endif
          write(*,*)"HDF reader not yet implemented"
          stop 1
        ELSE
          write(*,*)"ERROR: Unknown MR_idataFormat:",MR_idataFormat
          stop 1
        ENDIF
      case DEFAULT
        write(*,*)"MR ERROR:  MR_iwind not 1:5."
        write(*,*)"           Exiting in MR_Read_Met_DimVars"
        stop 1
      end select

      write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_Met_DimVars

!##############################################################################


!##############################################################################
!
!     MR_Set_CompProjection
!
!     This subroutine only sets the projection parameters for the computational
!     grid for MetReader.  The library needs to know these parameters in order
!     to know how to map the met grid onto the comp grid and whether or not to
!     rotate wind vectors.
!
!     This subroutine can be called from the calling program at any time, but
!     must be called before MR_Initialize_Met_Grids
!
!##############################################################################

      subroutine MR_Set_CompProjection(LL_flag,ipf,lam0,phi0,phi1,phi2,ko,Re)

      implicit none

      logical     ,intent(in) :: LL_flag
      integer     ,intent(in) :: ipf
      real(kind=8),intent(in) :: lam0,phi0
      real(kind=8),intent(in) :: phi1
      real(kind=8),intent(in) :: phi2
      real(kind=8),intent(in) :: ko
      real(kind=8),intent(in) :: Re

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------      MR_Set_CompProjection                                 ----------"
      write(*,*)"--------------------------------------------------------------------------------"

      isLatLon_CompGrid = LL_flag

      Comp_lam0 = 0.0
      Comp_phi0 = 0.0
      Comp_phi1 = 0.0
      Comp_phi2 = 0.0
      Comp_k0   = 0.0
      Comp_Re   = 0.0

      IF(.not.isLatLon_CompGrid)THEN
        Comp_iprojflag = ipf 
        IF(Comp_iprojflag.eq.1)THEN
          ! Polar stereographic
          Comp_lam0 = lam0
          Comp_phi0    = phi0
          Comp_k0      = ko
          Comp_Re      = Re
        ELSEIF(Comp_iprojflag.eq.2)THEN
          ! Albers Equal Area
          Comp_lam0 = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
        ELSEIF(Comp_iprojflag.eq.3)THEN
          ! UTM
          write(6,*)"WARNING: UTM not yet verified"
          stop 1
        ELSEIF(Comp_iprojflag.eq.4)THEN
          ! Lambert conformal conic (NARR, NAM218, NAM221)
          Comp_lam0 = lam0
          Comp_phi0    = phi0
          Comp_phi1    = phi1
          Comp_phi2    = phi2
          Comp_Re      = Re
        ELSEIF(Comp_iprojflag.eq.5)THEN
          ! Mercator (NAM196)
          Comp_lam0 = lam0
          Comp_phi0    = phi0
          Comp_Re      = Re
        ENDIF
      ENDIF

      write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Set_CompProjection

!##############################################################################
!
!     MR_Initialize_Met_Grids
!
!     This subroutine sets up local copies of the computational grid and
!     evaluates the full NWP grid for the subgrid needed, through MR_Set_MetComp_Grids
!     for NWP cases and MR_Set_MetComp_Grids_1dascii for 1-d or sonde cases.
!
!     From the calling program, this must be called after all of:
!       MR_Allocate_FullMetFileList
!       MR_Read_Met_DimVars
!       MR_Set_CompProjection
!
!     In some cases, it is useful to have velocity values saved within MetReader.
!     For example, velocities might be read to calculate diffusivities, stored
!     locally, then reused last for advection.  If the calling program sets the
!     variable MR_Save_Velocities=.true. , then space is allocated for the local
!     copies of velocity.
!
!     Takes as input :: specs of computational grid defined in the calling program
!     Sets: n[x,y,z]_comp and [x,y,z]_comp_sp :: local variables holding computational grid info
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           rdphi_MetP_sp      :: length scale along y (in meters)
!           rdlambda_MetP_sp   :: length scale along x (in meters)
!
!##############################################################################

      subroutine MR_Initialize_Met_Grids(nx,ny,nz, &
                                      dumx_sp,dumy_sp,dumz_sp,periodic)

      implicit none

      integer      ,intent(in) :: nx,ny,nz
      real(kind=sp),intent(in) :: dumx_sp(nx)
      real(kind=sp),intent(in) :: dumy_sp(ny)
      real(kind=sp),intent(in) :: dumz_sp(nz)
      logical      ,intent(in) :: periodic

      integer       :: i,j,k

      !                allocates *_Met_P for subset of Met grid on pressure levels
      !                          *_Met_H for subset of Met grid on height levels
      !                          *_comp_H for data regridded to computational grid

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------      MR_Initialize_Met_Grids                               ----------"
      write(*,*)"--------------------------------------------------------------------------------"

      nx_comp = nx
      ny_comp = ny
      nz_comp = nz
      allocate(x_comp_sp(nx_comp)); x_comp_sp = dumx_sp
      allocate(y_comp_sp(ny_comp)); y_comp_sp = dumy_sp
      allocate(z_comp_sp(nz_comp)); z_comp_sp = dumz_sp
      dx_comp = x_comp_sp(2) - x_comp_sp(1)
      dy_comp = abs(y_comp_sp(2) - y_comp_sp(1))
      MaxZ_comp_sp = maxval(z_comp_sp)

      IsPeriodic_CompGrid = periodic

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Set_MetComp_Grids_ASCII_1d
      case(2)
        !call Set_MetComp_Grids_3dascii
      case (3:5)
        ! Now that we have the full grids defined in MR_Read_Met_DimVars_netcdf,
        ! calculate the subgrid needed for the simulation
        call MR_Set_MetComp_Grids

        IF(MR_Save_Velocities)THEN
          write(*,*)"Velocities will be saved on the metP grid"
          allocate(MR_vx_metP_last(nx_submet,ny_submet,np_fullmet))
          allocate(MR_vx_metP_next(nx_submet,ny_submet,np_fullmet))
          allocate(MR_vy_metP_last(nx_submet,ny_submet,np_fullmet))
          allocate(MR_vy_metP_next(nx_submet,ny_submet,np_fullmet))
        ELSE
          write(*,*)"Velocities not saved"
        ENDIF

      case DEFAULT
        write(*,*)"MR ERROR:  MR_iwind not 1:5."
        write(*,*)"           Exiting in MR_Initialize_Met_Grids"
        stop 1
      end select

      ! For LatLon Met grids, calculate some additional geometry terms
      !  These are currently only needed for calculating DelMetP_Dx and
      !  DelMetP_Dy
      IF(IsLatLon_MetGrid)THEN
        allocate(rdphi_MetP_sp(ny_submet,np_fullmet))
        allocate(rdlambda_MetP_sp(nx_submet,ny_submet,np_fullmet))
        do k=1,np_fullmet
          IF(IsRegular_MetGrid)THEN
            ! length scale along y (in meters)
            rdphi_MetP_sp(:,k) = dy_met_const*DEG2RAD_MET * (RAD_EARTH_MET+z_approx(k))*1000.0_sp
            do j=1,ny_submet
              ! length scale along x (in meters)
              rdlambda_MetP_sp(:,j,k) =(RAD_EARTH_MET+z_approx(k))*1000.0_sp * &
                                      cos(DEG2RAD_MET*(y_submet_sp(j)-0.5_sp*dy_met_const)) * &
                                      dx_met_const*DEG2RAD_MET
            enddo
          ELSE
            do i=1,nx_submet
              do j=1,ny_submet
                ! length scale along y (in meters)
                rdphi_MetP_sp(:,k) = MR_dy_submet(j)*DEG2RAD_MET * (RAD_EARTH_MET+z_approx(k))*1000.0_sp
                ! length scale along x (in meters)
                rdlambda_MetP_sp(i,j,k) =(RAD_EARTH_MET+z_approx(k))*1000.0_sp * &
                                        cos(DEG2RAD_MET*(y_submet_sp(j)-0.5_sp*MR_dy_submet(j))) * &
                                        MR_dx_submet(i)*DEG2RAD_MET
              enddo
            enddo
          ENDIF
        enddo
      ENDIF

      write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Initialize_Met_Grids

!##############################################################################
!
!     MR_Set_Met_Times
!
!     This subroutine opens each file of the MR_windfile list, determines the hour
!     of each time step and logs the files/time-steps needed to bracket the
!     simulation duration.
!
!     This subroutine is called once from Initialize_Met_Grids.
!
!     Sets: MetStep_File(i=1:MR_MetSteps_Total)   :: contains file name for met step i
!           MetStep_findex(i=1:MR_MetSteps_Total) :: contains the index of the file in
!                                            the MR_windfiles(:) list
!           MetStep_tindex(i=1:MR_MetSteps_Total) :: contains the index within the file
!                                            for met step i
!
!##############################################################################

      subroutine MR_Set_Met_Times(eStartHour,Duration)

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      real(kind=8),intent(in) :: eStartHour
      real(kind=8),intent(in) :: Duration

      integer :: i
      integer :: iw
      real(kind=8)       :: HS_HourOfDay
      integer            :: HS_YearOfEvent
      integer            :: HS_MonthOfEvent
      integer            :: HS_DayOfEvent
      integer            :: HS_DayOfYear
      integer            :: iwstep
      integer            :: istep
      real(kind=8)       :: stephour
      logical :: Found_First_Step
      logical :: Found_Last_Step
      integer :: nMetSteps_Comp
      real(kind=8) :: met_t1,met_dt1
      logical      :: prestep, poststep

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------      MR_Set_Met_Times                                      ----------"
      write(*,*)"--------------------------------------------------------------------------------"

      MR_Comp_StartHour        = eStartHour
      MR_Comp_Time_in_hours    = Duration
      MR_Comp_StartYear        = HS_YearOfEvent(MR_Comp_StartHour,MR_BaseYear,MR_useLeap)

      ! Now we want to loop through all the steps of each windfile and find the
      ! file/step that immediately preceeds the time needed (MR_Comp_StartHour)
      !  First, define MR_Comp_StartHour if we are running a forecast run

      ! Now is a good time to check to make sure the MR_Comp_StartHour is within
      ! the range of data
      !   Checking if a prestep is needed
      IF(MR_iwind.eq.1.and.MR_iwindformat.eq.1)THEN
        ! For the ASCII profile cases (not the radiosonde), hours are given as offset
        ! from the start time.  So reset all hours to relative to eStartHour
        write(*,*)"Note:  The hours value in 1d ascii profiles are interpreted as offset"
        write(*,*)"       hours.  Shifting the file time to reflect the requested start"
        write(*,*)"       time plus offset."
        MR_windfile_starthour = MR_windfile_starthour + MR_Comp_StartHour
      ENDIF
      met_t1  = MR_windfile_starthour(1)+MR_windfile_stephour(1,1)
      met_dt1 = MR_ForecastInterval
      IF(MR_Comp_StartHour.lt.met_t1)THEN
        ! Start time requested is before that available, check if we can extrapolate
        IF(MR_Comp_StartHour.ge.met_t1-met_dt1)THEN
          write(*,*)"WARNING: Start time is before the first time step"
          write(*,*)"         However, it is within one interval so we will"
          write(*,*)"         apply the value at MR_windfile_starthour(1)."
          write(*,*)"         Using a time step interval of ",MR_ForecastInterval
          IF(MR_ForecastInterval.gt.720.0_dp)THEN
            IF(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)THEN
              write(*,*)"         Note: This interval is set high since only one"
              write(*,*)"               sonde file was provided."
            ENDIF
          ENDIF
          prestep = .true.
        ELSE
          write(*,*)"MR ERROR: Start time is before the first available data and"
          write(*,*)"       cannot be extrapolated."
          stop 1
        ENDIF
      ELSE
        write(*,*)"Start hour of simulation is in range of NWP data"
        prestep  = .false.
      ENDIF
      !   Checking if a poststep is needed
      met_t1  = MR_windfile_starthour(MR_iwindfiles)+MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
      met_dt1 = MR_ForecastInterval
      IF(MR_Comp_StartHour+MR_Comp_Time_in_hours.ge.met_t1)THEN
        ! Start time requested is after that available, check if we can extrapolate
        IF(MR_Comp_StartHour+MR_Comp_Time_in_hours.le.met_t1+met_dt1)THEN
          write(*,*)"WARNING: End time is at or after the last time step."
          write(*,*)"         However, it is within one interval so we will"
          write(*,*)"         apply the value at MR_windfile_starthour(MR_iwindfiles)"
          IF(MR_ForecastInterval.gt.720.0_dp)THEN
            IF(MR_iwind.eq.1.and.MR_iwindformat.eq.1.and.nt_fullmet.eq.1)THEN
              write(*,*)"         Note: This interval is set high since only one"
              write(*,*)"               sonde file was provided."
            ENDIF
          ENDIF
          poststep = .true.
        ELSE
          write(*,*)"MR ERROR: End time is after the last available data and"
          write(*,*)"       cannot be extrapolated."
          stop 1
        ENDIF
      ELSE
        write(*,*)"End hour of simulation is in range of NWP data"
        poststep  = .false.
      ENDIF

      ! Loop through all the files and steps and count how many are needed to bracket the time
      ! needed (MR_Comp_StartHour -> MR_Comp_StartHour+MR_Comp_Time_in_hours)
      ! Note: If preset or poststep is needed, istep will be incremented accordingly
      ! Once we know the number of steps needed, we will allocate space, then fill the variables
      ! with just the step info needed
      write(6,*)'    File num  | Step in file  |  stephour   |   SimStartHour   | nMetStep | Note'
      If(prestep)THEN
        Found_First_Step = .true.
        istep = 1
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - MR_ForecastInterval
        write(6,150)1,1,stephour,MR_Comp_StartHour,istep,"Prestep before MR_Comp_StartHour "
      ELSE
        Found_First_Step = .false.
        istep = 0
      ENDIF
      Found_Last_Step  = .false.
      DO iw = 1,MR_iwindfiles
        DO iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          IF(stephour.lt.MR_Comp_StartHour)THEN
            Found_First_Step = .false.
            istep = 1
            write(6,150)iw,iwstep,stephour,MR_Comp_StartHour,istep,"Before or at MR_Comp_StartHour."
          ELSE
            !  Otherwise, increment index if we are still in the needed time bracket
            IF(.not.Found_Last_Step)THEN
              Found_First_Step = .true.
              istep = istep + 1
              write(6,150)iw,iwstep,stephour,MR_Comp_StartHour,istep,"After MR_Comp_StartHour "
            ENDIF
          ENDIF

          ! Check if the current step hour is exclusively after the end of the simulation
          IF(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)THEN
            IF(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          ENDIF
          IF(Found_Last_Step)  &
            write(6,150)iw,iwstep,stephour,MR_Comp_StartHour,nMetSteps_Comp,"At or after END OF SIM  "
        ENDDO
      ENDDO
      ! If we went through all the steps and didn't find the last step, and if poststep=T, then
      ! increment istep and set nMetSteps_Comp
      IF(.not.Found_Last_Step)THEN
        IF(poststep)THEN
          nMetSteps_Comp = istep+1
          stephour = MR_windfile_starthour(MR_iwindfiles) + &
                       MR_windfile_stephour(MR_iwindfiles,MR_windfiles_nt_fullmet(MR_iwindfiles)) + &
                       MR_ForecastInterval
          write(6,150)1,1,stephour,MR_Comp_StartHour,nMetSteps_Comp,"Poststep after END OF SIM  "

        ELSE
          write(*,*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
          write(*,*)"           the simulation."
          stop 1
        ENDIF
      ENDIF
 150  FORMAT(8x,i3,9x,i4,4x,f15.2,2x,f13.2,10x,i3,8x,a25)

      ! We now have the number of steps needed for the computation
      ! Allocate the lists
      MR_MetSteps_Total = nMetSteps_Comp
      write(*,*)"MR: Allocating space for ",MR_MetSteps_Total,"steps"
      IF(prestep.or.poststep)THEN
        write(*,*)"         Including:"
        IF(prestep) write(*,*)"             1 prestep"
        IF(poststep)write(*,*)"             1 poststep"
      ENDIF
      allocate(MR_MetStep_File(MR_MetSteps_Total))
      allocate(MR_MetStep_findex(MR_MetSteps_Total))
      allocate(MR_MetStep_tindex(MR_MetSteps_Total))
      allocate(MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total))
      allocate(MR_MetStep_Interval(MR_MetSteps_Total))
      allocate(MR_MetStep_year(MR_MetSteps_Total))
      allocate(MR_MetStep_month(MR_MetSteps_Total))
      allocate(MR_MetStep_day(MR_MetSteps_Total))
      allocate(MR_MetStep_DOY(MR_MetSteps_Total))
      allocate(MR_MetStep_Hour_Of_Day(MR_MetSteps_Total))
      allocate(MR_iwind5_year(MR_MetSteps_Total))

      ! Finally, we need to loop through the steps exactly as above, but this time populate 
      ! the lists just allocated
      if(prestep)THEN
        Found_First_Step = .true.
        istep = 1
        stephour = MR_windfile_starthour(1) + MR_windfile_stephour(1,1) - MR_ForecastInterval
        MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(1)))
        MR_MetStep_findex(istep)          = 1
        MR_MetStep_tindex(istep)          = 1
        MR_MetStep_Hour_since_baseyear(istep) = stephour
        MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
        MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
        MR_iwind5_year(istep)             = MR_Comp_StartYear
      ELSE
        Found_First_Step = .false.
        istep = 0
      ENDIF
      Found_Last_Step  = .false.
      DO iw = 1,MR_iwindfiles
        DO iwstep = 1,MR_windfiles_nt_fullmet(iw)
          stephour = MR_windfile_starthour(iw) + MR_windfile_stephour(iw,iwstep)
          ! Unless the start time is before this step hour, reset the index istep to 1
          !  Otherwise, increment index
          IF(stephour.lt.MR_Comp_StartHour)THEN
            Found_First_Step = .false.
            istep = 1
          ELSE
            !  Otherwise, increment index if we are still in the needed time bracket
            IF(.not.Found_Last_Step)THEN
              Found_First_Step = .true.
              istep = istep + 1
            ENDIF
          ENDIF

          IF(.not.Found_Last_Step)THEN
            MR_MetStep_File(istep)            = trim(adjustl(MR_windfiles(iw)))
            MR_MetStep_findex(istep)          = iw
            MR_MetStep_tindex(istep)          = iwstep
            MR_MetStep_Hour_since_baseyear(istep) = stephour
            MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
            MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
            MR_iwind5_year(istep)             = MR_MetStep_year(istep)
          ENDIF


          ! Check if the current step hour is exclusively after the end of the simulation
          IF(stephour.gt.MR_Comp_StartHour+MR_Comp_Time_in_hours)THEN
            IF(.not.Found_Last_Step)nMetSteps_Comp=istep
            Found_Last_Step = .true.
          ENDIF
        ENDDO
      ENDDO
      IF(.not.Found_Last_Step)THEN
        IF(poststep)THEN
          istep = istep+1
          stephour = MR_MetStep_Hour_since_baseyear(istep-1) + MR_ForecastInterval
          MR_MetStep_File(istep)            = MR_MetStep_File(istep-1)
          MR_MetStep_findex(istep)          = MR_MetStep_findex(istep-1)
          MR_MetStep_tindex(istep)          = MR_MetStep_tindex(istep-1)
          MR_MetStep_Hour_since_baseyear(istep) = stephour
          MR_MetStep_year(istep)            = HS_YearOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_month(istep)           = HS_MonthOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_day(istep)             = HS_DayOfEvent(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_DOY(istep)             = HS_DayOfYear(real(stephour,kind=dp),MR_BaseYear,MR_useLeap)
          MR_MetStep_Hour_Of_Day(istep)     = real(HS_HourOfDay(real(stephour,kind=dp),MR_BaseYear,MR_useLeap),kind=sp)
          MR_iwind5_year(istep)             = MR_MetStep_year(istep)
        ELSE
          write(*,*)"MR ERROR:  Something is wrong.  Could not find the last MetStep needed for"
          write(*,*)"           the simulation."
          stop 1
        ENDIF
      ENDIF

      DO istep = 1,MR_MetSteps_Total-1
        MR_MetStep_Interval(istep) = MR_MetStep_Hour_since_baseyear(istep+1) - &
                                     MR_MetStep_Hour_since_baseyear(istep)
      ENDDO
      MR_MetStep_Interval(MR_MetSteps_Total)=MR_MetStep_Interval(MR_MetSteps_Total-1)

      write(*,*)"Comp start = ",real(MR_Comp_StartHour,kind=4)
      write(*,*)"Comp end   = ",real(MR_Comp_StartHour+MR_Comp_Time_in_hours,kind=4)
      DO istep = 1,MR_MetSteps_Total
        iw       = MR_MetStep_findex(istep)
        iwstep   = MR_MetStep_tindex(istep)
        stephour = MR_MetStep_Hour_since_baseyear(istep)
        write(*,*)istep,trim(adjustl(MR_MetStep_File(istep))),iw,iwstep,real(stephour,kind=4)
      ENDDO

 !     MAKE SURE THE WIND MODEL TIME WINDOW COVERS THE ENTIRE SUMULATION TIME
      write(6,99)
      !write(9,99)
99    format(/,4x,'Making sure the mesoscale model time covers the simulation time . . . ')
      if (MR_MetStep_Hour_since_baseyear(1).gt.MR_Comp_StartHour) then
        write(*,*)"MR ERROR:  Wind data starts after SimStartHour"
        write(*,*)"WindHour(MR_iwindfiles) = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
        write(*,*)"SimStartHour            = ",MR_Comp_StartHour
        write(*,*)"Simtime_in_hours        = ",MR_Comp_Time_in_hours
        write(*,*)"  All steps:"
        DO i = 1,MR_MetSteps_Total
          write(*,*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
        ENDDO
        !write(9,101)
!101     format(4x,'Error.  The model starts after the first eruption. Program stopped.')
        stop 1
      else if (MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total).lt.(MR_Comp_StartHour+MR_Comp_Time_in_hours)) then
        write(*,*)"MR ERROR:  Last met time is earlier than simulation time"
        write(*,*)"MR_MetSteps_Total     = ",MR_MetSteps_Total
        write(*,*)"Last time step       = ",MR_MetStep_Hour_since_baseyear(MR_MetSteps_Total)
        write(*,*)"SimStartHour         = ",MR_Comp_StartHour
        write(*,*)"Simtime_in_hours     = ",MR_Comp_Time_in_hours
        DO i = 1,MR_MetSteps_Total
          write(*,*)"    ",i,MR_MetStep_Hour_since_baseyear(i)
        ENDDO
        write(9,102)
102     format(4x,'Error.  The model ends before the end of the last eruption.  Program stopped.')
        stop 1
      end if
      write(6,103)
      !write(9,103)
103   format(4x,'Good.  It does.',/)

      write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_Met_Times


!##############################################################################
!
!     MR_Read_HGT_arrays
!
!     This subroutine basically does the same thing as Read_3d_MetP_Variable
!     but specificly for the variable HGT (ivar=1).  This is needed because the
!     more general Read_3d_MetP_Variable just returns data in dum3d_metP and
!     expects the calling program to save the results.  The HGT variables,
!     geoH_metP_last and geoH_metP_next need to persist locally, however,  since they are
!     used in QC calculations and are needed in converting PresVertVel (from Pa
!     s to m/s).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: istep :: specified the met step
!
!     Sets: geoH_metP_last, geoH_metP_next
!
!##############################################################################

      subroutine MR_Read_HGT_arrays(istep,reset_first_time)

      implicit none

      integer,intent(in)           :: istep
      logical, optional,intent(in) :: reset_first_time

      integer :: ivar
      logical,save :: first_time = .true.
      integer :: k

      if(present(reset_first_time)) then
        first_time = .true.
      endif

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_HGT_arrays                                    ----------"
      !write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"     Reading HGT array for istep = ",istep

      ivar = 1 ! HGT
      IF(first_time)THEN
        IF(MR_idataFormat.eq.1)THEN
          call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
        ELSEIF(MR_idataFormat.eq.2)THEN
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.3)THEN
#ifdef USEGRIB2
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.4)THEN
#ifdef USEHDF
          !call MR_Read_MetP_Variable_HDF(ivar,istep)
#endif
        ENDIF
        MR_geoH_metP_last(:,:,:) = MR_dum3d_metP(:,:,:)
        first_time = .false.
        Max_geoH_metP_last = maxval(MR_geoH_metP_last(:,:,np_fullmet))
      ELSE
        MR_geoH_metP_last(:,:,:) = MR_geoH_metP_next(:,:,:)
        Max_geoH_metP_last = Max_geoH_metP_next
        Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
      ENDIF
      IF(MR_idataFormat.eq.1)THEN
          call MR_Read_MetP_Variable_ASCII_1d(ivar,istep+1)
      ELSEIF(MR_idataFormat.eq.2)THEN
#ifdef USENETCDF
        call MR_Read_MetP_Variable_netcdf(ivar,istep+1)
#endif
      ELSEIF(MR_idataFormat.eq.3)THEN
#ifdef USEGRIB2
        call MR_Read_MetP_Variable_GRIB(ivar,istep+1)
#endif
      ELSEIF(MR_idataFormat.eq.4)THEN
#ifdef USEHDF
        !call MR_Read_MetP_Variable_HDF(ivar,istep+1)
#endif
      ENDIF
      MR_geoH_metP_next(:,:,:) = MR_dum3d_metP(:,:,:)

      Max_geoH_metP_next = maxval(MR_geoH_metP_next(:,:,np_fullmet))
        ! Now determine the maximum value provided between last and next steps
      Max_geoH_metP = max(Max_geoH_metP_last,Max_geoH_metP_next)
        ! and compare with the maximun value needed by the computational grid
      IF(Max_geoH_metP.lt.z_comp_sp(nz_comp))THEN
        ! use highest needed point (rounded up to 5-km increment) for the
        ! padding height
        Suppl_H = real(ceiling(z_comp_sp(nz_comp)/5.0_sp),kind=sp) * 5.0_sp
      ELSE
        ! otherwise use the height from the met files to determine padding
        Suppl_H = real(ceiling(Max_geoH_metP/5.0_sp),kind=sp) * 5.0_sp
      ENDIF
      IF(MR_iHeightHandler.eq.1)THEN
        If(z_comp_sp(nz_comp).gt.Max_geoH_metP)THEN
          write(*,*)"Computational grid extends higher than met grid"
          write(*,*)"iHeightHandler = 1:  exiting"
          stop 1
        ENDIF
      ENDIF

      !write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_HGT_arrays


!##############################################################################
!
!     MR_Read_3d_MetP_Variable
!
!     This subroutine extracts the variable Met_var_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid on pressure coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets: dum3d_metP
!##############################################################################

      subroutine MR_Read_3d_MetP_Variable(ivar,istep)

      implicit none

      integer,intent(in)        :: ivar,istep

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_3d_MetP_Variable                              ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        call MR_Read_MetP_Variable_ASCII_1d(ivar,istep)
      case(2)
        !call Read_3d_ASCII
      case (3:5)
        IF(MR_idataFormat.eq.2)THEN
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.3)THEN
#ifdef USEGRIB2
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.4)THEN
#ifdef USEHDF
          !call MR_Read_MetP_Variable_HDF(ivar,istep)
#endif
        ENDIF

      case DEFAULT

      end select

      !write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_3d_MetP_Variable

!##############################################################################
!
!     MR_Read_3d_MetH_Variable
!
!     This subroutine extracts the variable Met_var_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Reads : MR_dum3d_metP
!     Sets  : MR_dum3d_metH
!##############################################################################

      subroutine MR_Read_3d_MetH_Variable(ivar,istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)        :: ivar,istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext
      integer :: np_fully_padded

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_3d_MetH_Variable                              ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      np_fully_padded = np_fullmet+1+np_fullmet_pad

      ! First get the variable on the native pressure coordinate
      call MR_Read_3d_MetP_Variable(ivar,istep)

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp
      ! Allocate 1-d arrays, padding one layer to the bottom and np_fullmet_pad
      ! layers on the top
      allocate(  z_col_metP(np_fully_padded))
      allocate(var_col_metP(np_fully_padded))
      allocate(var_col_metH(nz_comp))

!     CREATE 1-D ARRAYS IN P, AND REGRID THEM INTO 1-D ARRAYS IN z
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          IF(istep.eq.MR_iMetStep_Now)THEN
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          ELSE
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          ENDIF
          ! Set first node at lower boundary
          z_col_metP(1)                = 0.0e-4_sp  ! 0.1 m
            ! Sometimes the lowest z from the wind file (k=2) might be
            ! negative.  If so, reassign to just above boundary.
          if (z_col_metP(2).le.2.0e-4_sp) z_col_metP(2) = 2.0e-4_sp
          if (z_col_metP(3).le.3.0e-4_sp) z_col_metP(3) = 3.0e-4_sp

          var_col_metP(1)              = MR_dum3d_metP(i,j,1)
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while ((z_col_metP(knext)<z_col_metP(k-1)).and. &
                        (knext.lt.np_fullmet+1))
                 knext=knext+1
              end do
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              end do
            endif
          enddo
          IF(np_fullmet_pad.gt.0)THEN
            ! Assign top node (Should these indicies be +1?)
            IF(MR_iHeightHandler.eq.1)THEN
              z_col_metP(np_fullmet+2) = z_col_metP(np_fullmet+1) + 10.0_sp
              write(*,*)z_col_metP(np_fullmet+1:np_fullmet+2)
            ELSEIF(MR_iHeightHandler.eq.2)THEN
              z_col_metP(np_fullmet+2) = Suppl_H
              IF(ivar.ne.5)THEN
                ! For all variables except temperature (ivar=5), just copy highest
                ! met node
                var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)
              ELSE
                ! For temperature, use the 1976 US Standard Atmosphere
                var_col_metP(np_fullmet+2) = MR_Temp_US_StdAtm(Suppl_H)
              ENDIF
            ENDIF
          ENDIF

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid
          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,  var_col_metP, & 
                                    nz_comp,       z_comp_sp,  var_col_metH)

          MR_dum3d_metH(i,j,:) = var_col_metH

        enddo ! j
      enddo  ! i

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Read_3d_MetH_Variable

!##############################################################################
!
!     MR_Read_3d_Met_Variable_to_CompGrid
!
!     This subroutine extracts the variable Met_var_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_3d_MetH_Variable                   (gets variable on Met_x, Met_y, Comp_z)
!        -> Read_3d_MetP_Variable               (gets variable on native Met subgrid)
!            -> Read_3d_MetP_Variable_[format]  (direct read of variable in
!                                                whatever format : nc,hdf,grib,ascii)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum3d_compH
!               MR_dum3d_metH and MR_dum3d_metP are also filled in the course of
!               generating MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Read_3d_Met_Variable_to_CompGrid(ivar,istep,IsNext)

      implicit none

      integer,intent(in)           :: ivar,istep
      logical, optional,intent(in) :: IsNext

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_3d_Met_Variable_to_CompGrid                   ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

        ! First get the variable on the height coordinate
      call MR_Read_3d_MetH_Variable(ivar,istep)

      IF(MR_Save_Velocities)THEN
        ! Note: this flag for saving the velocity values is useful in special
        ! cased such where the velocities might be read and used for a local
        ! calculation, but then can be used later.  Variable diffusivity uses
        ! this.
        if(present(IsNext)) then
          ! MR_dum3d_metP still contains the variable just read
          IF(ivar.eq.2)THEN
            IF(IsNext)THEN
              MR_vx_metP_last = MR_vx_metP_next
              MR_vx_metP_next = MR_dum3d_metP
            ELSE
              MR_vx_metP_last = MR_dum3d_metP
            ENDIF
          ELSEIF(ivar.eq.3)THEN
            IF(IsNext)THEN
              MR_vy_metP_last = MR_vy_metP_next
              MR_vy_metP_next = MR_dum3d_metP
            ELSE
              MR_vy_metP_last = MR_dum3d_metP
            ENDIF
          ENDIF
        ENDIF
      ENDIF
        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp))

      DO k=1,nz_comp
        IF(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)THEN
          !HFS This will not work for multi-site sonde data
          MR_dum3d_compH(:,:,k) = MR_dum3d_metH(1,1,k)
          cycle
        ENDIF

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

        DO i = 1,nx_comp
          DO j = 1,ny_comp
            IF(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          ENDDO
        ENDDO
        MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
      ENDDO

      deallocate(tmp_regrid2d_sp)

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Read_3d_Met_Variable_to_CompGrid


!##############################################################################
!
!     MR_Read_2d_Met_Variable
!
!     This subroutine extracts the variable Met_var_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on height coordinates.
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_met or MR_dum2d_met_int
!##############################################################################

      subroutine MR_Read_2d_Met_Variable(ivar,istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)        :: ivar,istep

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_2d_Met_Variable                               ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      select case (MR_iwind)
      case(1)   ! if we're using a 1-D wind sounding
        !call Read_1d_windfile
      case(2)
        !call Read_2d_ASCII
      case (3:5)
        IF(MR_idataFormat.eq.2)THEN
#ifdef USENETCDF
          call MR_Read_MetP_Variable_netcdf(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.3)THEN
#ifdef USEGRIB2
          call MR_Read_MetP_Variable_GRIB(ivar,istep)
#endif
        ELSEIF(MR_idataFormat.eq.4)THEN
#ifdef USEHDF
          !call MR_Read_MetP_Variable_HDF(ivar,istep)
#endif
        ENDIF
      case DEFAULT

      end select

      !write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_2d_Met_Variable

!##############################################################################
!
!     MR_Read_2d_Met_Variable_to_CompGrid
!
!     This subroutine extracts the variable Met_var_names(ivar) from the
!     windfile/timestep given by MetStep_File(istep),MetStep_tindex(istep).
!     The values extracted are just on the needed subgrid of the full met
!     grid remapped on to the computational grid.  This is done by calling:
!       Read_2d_Met_Variable                   (gets variable on Met_x, Met_y)
!          -> Read_Met_Variable_[format]       (direct read of variable in
!                                                whatever format :
!                                                nc,hdf,grib,ascii)
!
!     Takes as input :: ivar  :: specifies which variable to read
!                       istep :: specified the met step
!     Sets  : MR_dum2d_comp
!               MR_dum2d_met and MR_dum2d_met are also filled in the course of
!               generating MR_dum2d_comp
!
!##############################################################################

      subroutine MR_Read_2d_Met_Variable_to_CompGrid(ivar,istep)

      implicit none

      integer,intent(in)        :: ivar,istep

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Read_2d_Met_Variable_to_CompGrid                   ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      call MR_Read_2d_Met_Variable(ivar,istep)

      allocate(tmp_regrid2d_sp(nx_comp,ny_comp))

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

      DO i = 1,nx_comp
        DO j = 1,ny_comp
          IF(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
        ENDDO
      ENDDO
      MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Read_2d_Met_Variable_to_CompGrid

!##############################################################################
!
!     MR_Rotate_UV_GR2ER_Met
!
!     This subroutine reads U and V on the metP grid, rotates to GridRelative
!
!     Takes as input :: istep :: specified the met step
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_GR2ER_Met(istep)

      implicit none

      integer,intent(in)        :: istep

      integer             :: i,j,k

      real(kind=sp) :: vx_old,vy_old
      real(kind=sp) :: vx_new,vy_new
      real(kind=sp) :: rotang

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Rotate_UV_GR2ER_Met                                ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      ! Rotate wind vectors on the MetP grid for NARR cases with Map_Case = 2 (both projected)
      ! or rotate to LL on MetP grid for Map_Case = 4 (Met=proj, Comp=LL)
      write(*,*)" Rotating Earth-relative projected Met winds to grid-relative"
      write(*,*)"  or rotating projected grid-relative Met winds to Earth-relative"

      IF(.not.MR_Save_Velocities)THEN
        write(*,*)"MR WARNING: Velocities not saved"
      ENDIF

      call MR_Read_3d_MetP_Variable(2,istep)
        MR_u_ER_metP = MR_dum3d_metP
      call MR_Read_3d_MetP_Variable(3,istep)
        MR_v_ER_metP = MR_dum3d_metP
      ! Using these earth relative velocities, rotate to grid relative and
      ! copy to MR_dum3d_metP
      DO i=1,nx_submet
        DO j=1,ny_submet
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids_netcdf
          rotang = real(theta_Met(i,j),kind=sp)
          DO k=1,np_fullmet
            vx_old = MR_u_ER_metP(i,j,k)
            vy_old = MR_v_ER_metP(i,j,k)
          ! Project vx_old onto grid
            vx_new = vx_old * cos(rotang)
            vy_new = vx_old * sin(rotang)
          ! Add projection of vel_old
            vx_new = vx_new - vy_old * sin(rotang)
            vy_new = vy_new + vy_old * cos(rotang)
            MR_u_ER_metP(i,j,k) = vx_new
            MR_v_ER_metP(i,j,k) = vy_new
          ENDDO
        ENDDO
      ENDDO
      !MR_dum3d_metP = MR_v_ER_metP
      !call MR_Regrid_MetP_to_CompGrid(istep)
      !MR_dum3d_compH_2 = MR_dum3d_compH
      !MR_dum3d_metP = MR_u_ER_metP
      !call MR_Regrid_MetP_to_CompGrid(istep)

      !write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Rotate_UV_GR2ER_Met


!##############################################################################
!
!     MR_Rotate_UV_ER2GR_Comp
!
!     This subroutine is needed when we have a Earth-relative Met grid (either
!     natively or through MR_Rotate_UV_GR2ER_Met) and need to rotate that into a
!     projected Comp grid.  This corresponds to:
!       Map_Case = 3 : Met=LL,   Comp=proj
!       Map_Case = 5 : Met=proj, Comp=proj (differentn)
!     This subroutine reads U and V on the metP projected grid,
!     then resamples onto the projected compH grid, rotates to grid relative and
!     returns to calling program
!
!     Takes as input :: istep :: specified the met step
!
!     Sets  : MR_dum3d_compH    holds U
!             MR_dum3d_compH_2  holds V
!
!##############################################################################

      subroutine MR_Rotate_UV_ER2GR_Comp(istep)

      implicit none

      integer,intent(in)        :: istep

      integer             :: i,j,k

      real(kind=sp) :: vx_old,vy_old
      real(kind=sp) :: vx_new,vy_new
      real(kind=sp) :: rotang1

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Rotate_UV_ER2GR_Met                                ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      IF(Map_Case.eq.3)THEN
        ! Met grid is natively LL and Comp grid is projected
        ! Fill the y velocities to spare array (comp_2)
        call MR_Read_3d_MetP_Variable(3,istep) ! This fills MR_dum3d_metP
        call MR_Regrid_MetP_to_CompGrid(istep) ! Takes MR_dum3d_metP and fills MR_dum3d_compH
          MR_dum3d_compH_2 = MR_dum3d_compH

        call MR_Read_3d_MetP_Variable(2,istep) 
        call MR_Regrid_MetP_to_CompGrid(istep)
        ! Now compH and compH_2 have vx and vy
      ELSEIF(Map_Case.eq.5)THEN
        ! Assume MR_u_ER_metP and MR_v_ER_metP have already been set by
        ! MR_Rotate_UV_GR2ER_Met
        MR_dum3d_metP = MR_v_ER_metP
        call MR_Regrid_MetP_to_CompGrid(istep)
          MR_dum3d_compH_2 = MR_dum3d_compH

        MR_dum3d_metP = MR_u_ER_metP
        call MR_Regrid_MetP_to_CompGrid(istep)
        ! Now compH and compH_2 have vx and vy
      ELSE
        write(*,*)"Should not be calling this subroutine unless Map_Case=3 or 5"
        stop 1
      ENDIF

      ! Now loop through the comp points and rotate the vectors in place
      DO i=1,nx_comp
        DO j=1,ny_comp
          ! The angle theta for the Earth to Grid conversion was
          ! precalculated in Set_MetComp_Grids_netcdf
          rotang1 = real(theta_Comp(i,j),kind=sp)

          DO k=1,nz_comp
            vx_old = MR_dum3d_compH(i,j,k)
            vy_old = MR_dum3d_compH_2(i,j,k)

            vx_new = vx_old * cos(rotang1) - vy_old * sin(rotang1)
            vy_new = vx_old * sin(rotang1) + vy_old * cos(rotang1)

            MR_dum3d_compH(i,j,k) = vx_new
            MR_dum3d_compH_2(i,j,k) = vy_new
          ENDDO
        ENDDO
      ENDDO

      !write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Rotate_UV_ER2GR_Comp

!##############################################################################
!
!     MR_Regrid_MetP_to_CompGrid
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH by the subroutine Regrid_MetP_to_MetH.
!     MR_dum3d_metH is then regridded onto MR_dum3d_compH
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_CompGrid(istep)

      implicit none

      integer,intent(in)  :: istep

      integer             :: i,j,k
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Regrid_MetP_to_CompGrid                            ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      ! convert MR_dum3d_MetP to MR_dum3d_metH
      call MR_Regrid_MetP_to_MetH(istep)

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp))

      DO k=1,nz_comp
        IF(MR_iwindformat.eq.1.or.MR_iwindformat.eq.2)THEN
          !HFS This will not work for multi-site sonde data
          MR_dum3d_compH(:,:,k) = MR_dum3d_metH(1,1,k)
          cycle
        ENDIF

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum3d_metH(1:nx_submet,1:ny_submet,k),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

        DO i = 1,nx_comp
          DO j = 1,ny_comp
            IF(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          ENDDO
        ENDDO
        MR_dum3d_compH(:,:,k) = tmp_regrid2d_sp(:,:)
      ENDDO

      deallocate(tmp_regrid2d_sp)

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Regrid_MetP_to_CompGrid

!##############################################################################
!
!     MR_Regrid_MetP_to_MetH
!
!     This subroutine expects the calling program to populate MR_dum3d_metP.
!     This is regridded onto MR_dum3d_metH .
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_Regrid_MetP_to_MetH(istep)

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer,intent(in)  :: istep

      real(kind=sp),dimension(:),allocatable :: z_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metP
      real(kind=sp),dimension(:),allocatable :: var_col_metH
      integer :: i,j,k
      integer :: kc,knext

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Regrid_MetP_to_MetH                                ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

      ! Now remap from
      !    nx_submet,ny_submet,np_fullmet to
      !    nx_submet,ny_submet,nz_comp

      allocate(  z_col_metP(np_fullmet+2))
      allocate(var_col_metP(np_fullmet+2))
      allocate(var_col_metH(nz_comp))

!     CREATE 1-D ARRAYS IN P, AND REGRID THEM INTO 1-D ARRAYS IN z
      do i=1,nx_submet
        do j=1,ny_submet
            ! copy the column of z values for this i,j
            ! Note: These are the height values from the windfile plus an extra
            ! point at the surface and an extra point above the wind grid.
          IF(istep.eq.MR_iMetStep_Now)THEN
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_last(i,j,1:np_fullmet)
          ELSE
            z_col_metP(2:np_fullmet+1)  = MR_geoH_metP_next(i,j,1:np_fullmet)
          ENDIF
          ! Set first node at lower boundary
          z_col_metP(1)                = 1.0e-4_sp  ! 0.1 m
            ! Sometimes the lowest z from the wind file (k=2) might be
            ! negative.  If so, reassign to just above boundary.
          if (z_col_metP(2).le.2.0e-4_sp) z_col_metP(2) = 2.0e-4_sp
          if (z_col_metP(3).le.3.0e-4_sp) z_col_metP(3) = 3.0e-4_sp

          var_col_metP(1)              = 0.0_sp
          var_col_metP(2:np_fullmet+1) = MR_dum3d_metP(i,j,1:np_fullmet)
          ! Fix occasional erroneous pressure levels that are non-increasing by
          ! resetting height to average of the neighbors in z, or, if more than
          ! one z value is lower than a previous one, reset all such values to
          ! increase linearly with k
          do k=2,np_fullmet
            if (z_col_metP(k)<z_col_metP(k-1))then
              knext=k       !find the first value of z_col_metP that's above z_col_metP(k-1)
              do while (z_col_metP(knext)<z_col_metP(k-1))
                 knext=knext+1
              end do
              do kc=k,knext-1    !correct all intervening values of z_col_metP
                 z_col_metP(kc)=z_col_metP(k-1)+(real(kc-(k-1),kind=sp)/real(knext-(k-1),kind=sp))* &
                                              (z_col_metP(knext)-z_col_metP(k-1))
              end do
            endif
          enddo
            ! Assign top node (Should these indicies be +1?)
          z_col_metP(np_fullmet+2) = max(z_comp_sp(nz_comp)+1.0_sp,z_col_metP(np_fullmet+1)+1.0_sp)
          var_col_metP(np_fullmet+2) = var_col_metP(np_fullmet+1)

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !   Interpolate these values to a regular grid with
          !   spacing equal to the simulation grid

          call MR_Regrid_P2H_linear(np_fullmet+2, z_col_metP,  var_col_metP, & !vx
                                    nz_comp,       z_comp_sp,  var_col_metH)

          MR_dum3d_metH(i,j,:) = var_col_metH
        enddo ! j
      enddo  ! i

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Regrid_MetP_to_MetH


!##############################################################################
!
!     MR_Regrid_Met2d_to_Comp2d
!
!     This subroutine expects the calling program to populate MR_dum2d_met.
!     MR_dum2d_met is then regridded onto MR_dum2d_comp
!
!     Takes as input :: istep :: specified the met step
!     Sets  : MR_dum3d_compH
!
!##############################################################################

      subroutine MR_Regrid_Met2d_to_Comp2d

      implicit none

      integer             :: i,j
      real(kind=sp),dimension(:,:),allocatable :: tmp_regrid2d_sp

      !write(*,*)"--------------------------------------------------------------------------------"
      !write(*,*)"----------      MR_Regrid_Met2d_to_Comp2d                             ----------"
      !write(*,*)"--------------------------------------------------------------------------------"

        ! Now we have MR_dum3d_metH; interpolate onto computational grid
        !  Since MR_dum3d_metH and MR_dum3d_compH have the same z-coordinate, we only
        !  need to do 2d regridding on each k-slice
      allocate(tmp_regrid2d_sp(nx_comp,ny_comp))

        call MR_Regrid_Met2Comp(nx_submet,ny_submet, MR_dum2d_met(1:nx_submet,1:ny_submet),       &
                                nx_comp,  ny_comp,   tmp_regrid2d_sp(1:nx_comp,1:ny_comp))

        DO i = 1,nx_comp
          DO j = 1,ny_comp
            IF(isnan(tmp_regrid2d_sp(i,j)))tmp_regrid2d_sp(i,j)=0.0_sp
          ENDDO
        ENDDO
        MR_dum2d_comp(:,:) = tmp_regrid2d_sp(:,:)

      deallocate(tmp_regrid2d_sp)

      !write(*,*)"--------------------------------------------------------------------------------"

      return

      end subroutine MR_Regrid_Met2d_to_Comp2d


!##############################################################################
!
!     MR_DelMetP_Dx
!
!     Calculated the x derivitive of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dx

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer :: i
      integer :: lside,rside
      real(kind=sp) :: dx_fac
      real(kind=sp) :: KM_2_M

      KM_2_M = 100.0_sp

      do i=1,nx_submet
        IF(i.eq.1)THEN
          IF(IsPeriodic_CompGrid)THEN
            ! Two-sided derivitive wrapping around
            lside  = nx_submet
            rside  = i+1
            dx_fac = 2.0_sp
          ELSE
            ! One-sided derivitive on left
            lside  = i+0
            rside  = i+1
            dx_fac = 1.0_sp
          ENDIF
        ELSEIF(i.eq.nx_submet)THEN
          IF(IsPeriodic_CompGrid)THEN
            ! Two-sided derivitive wrapping around
            lside  = i-1
            rside  = 1
            dx_fac = 2.0_sp
          ELSE
            ! One-sided derivitive on right
            lside  = i-1
            rside  = i+0
            dx_fac = 1.0_sp
          ENDIF
        ELSE
          ! Two-sided derivitive
          lside  = i-1
          rside  = i+1
          dx_fac = 2.0_sp
        ENDIF
        IF(IsLatLon_MetGrid)THEN
          IF(IsRegular_MetGrid)THEN
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (rdlambda_MetP_sp(i,:,:)   * &
                                    dx_fac*KM_2_M)
          ELSE
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (0.5_sp*(rdlambda_MetP_sp(i,:,:)+rdlambda_MetP_sp(rside,:,:)) * &
                                    dx_fac*KM_2_M)
          ENDIF
        ELSE
          IF(IsRegular_MetGrid)THEN
            MR_dum3d2_metP(i,:,:) = (MR_dum3d_metP(lside,:,:)  + &
                                     MR_dum3d_metP(rside,:,:)) / &
                                    (dx_fac*MR_dx_submet(i)*KM_2_M)
          ELSE
            write(*,*)"Need to fix DelMetP_Dx for non-regular grids."
            stop 1
          ENDIF
        ENDIF
      enddo


      end subroutine MR_DelMetP_Dx


!##############################################################################
!
!     MR_DelMetP_Dy
!
!     Calculated the y derivitive of the variable in MR_dum3d_MetP.
!     This subroutine expects the calling program to populate MR_dum3d2_metP.
!
!     Sets  : MR_dum3d_metH
!
!##############################################################################

      subroutine MR_DelMetP_Dy

      implicit none

      integer, parameter :: sp        = 4 ! single precision

      integer :: i,j
      integer :: lside,rside
      real(kind=sp) :: dy_fac
      real(kind=sp) :: KM_2_M

      KM_2_M = 100.0_sp

      do j=1,ny_submet
        IF(j.eq.1)THEN
          ! One-sided derivitive on left
          lside=j+0
          rside=j+1
          dy_fac = 1.0_sp
        ELSEIF(j.eq.ny_submet)THEN
          ! One-sided derivitive on right
          lside=j-1
          rside=j+0
          dy_fac = 1.0_sp
        ELSE
          ! Two-sided derivitive
          lside=j-1
          rside=j+1
          dy_fac = 2.0_sp
        ENDIF
        IF(IsLatLon_MetGrid)THEN
          IF(IsRegular_MetGrid)THEN
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,lside,:)  + &
                                       MR_dum3d_metP(i,rside,:)) / &
                                      (rdphi_MetP_sp(j,:)   * &
                                       dy_fac*KM_2_M)
            enddo
          ELSE
            do i=1,nx_submet
              MR_dum3d2_metP(i,j,:) = (MR_dum3d_metP(i,lside,:)  + &
                                       MR_dum3d_metP(i,rside,:)) / &
                                      (0.5_sp*(rdphi_MetP_sp(j,:)+rdphi_MetP_sp(rside,:))   * &
                                       dy_fac*KM_2_M)
            enddo
          ENDIF
        ELSE
          IF(IsRegular_MetGrid)THEN
            MR_dum3d2_metP(:,j,:) = (MR_dum3d_metP(:,lside,:)  + &
                                     MR_dum3d_metP(:,rside,:)) / &
                                    (dy_fac*MR_dy_submet(j)*KM_2_M)
          ELSE
            write(*,*)"Need to fix DelMetP_Dy for non-regular grids."
            stop 1
          ENDIF
        ENDIF
      enddo

      end subroutine MR_DelMetP_Dy

!##############################################################################
!
!  Returns the temperature (in K) given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!
!##############################################################################

      function MR_Temp_US_StdAtm(zin)

      implicit none

      real(kind=sp) :: MR_Temp_US_StdAtm
      real(kind=sp) :: zin

      real(kind=sp),dimension(7) :: US_StdAtm_znodes
      real(kind=sp),dimension(7) :: US_StdAtm_Tnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_temp
      integer :: k,kk

      US_StdAtm_znodes = (/ 0.0_sp, 11.0_sp, 20.0_sp, 32.0_sp, &
                           47.0_sp, 51.0_sp, 71.0_sp/)
      US_StdAtm_Tnodes = (/288.15_sp, 216.65_sp, 216.65_sp, 228.65_sp, &
                           270.65_sp, 270.65_sp, 214.65_sp/)

      IF(zin.le.US_StdAtm_znodes(1))THEN
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(1)
      ELSEIf(zin.ge.US_StdAtm_znodes(7))THEN
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(7)
        write(*,*)"Stopping in Temp_US_StdAtm"
        stop 1
      ELSE
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        DO k = 6,1,-1
          If(zin.ge.US_StdAtm_znodes(k).and.zin.lt.US_StdAtm_znodes(k+1)) kk = k
        ENDDO
        frac = (zin-US_StdAtm_znodes(kk)) / &
               (US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk))
        Delta_temp = US_StdAtm_Tnodes(kk+1)-US_StdAtm_Tnodes(kk)
        MR_Temp_US_StdAtm = US_StdAtm_Tnodes(kk) + Delta_temp * frac
      ENDIF

      return

      end function MR_Temp_US_StdAtm

!##############################################################################
!
!  Returns the height in km given the pressure in hPa
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!
!##############################################################################

      function MR_Z_US_StdAtm(pin)

      implicit none

      real(kind=sp) :: MR_Z_US_StdAtm  ! in km
      real(kind=sp) :: pin          ! in mb

      real(kind=sp),dimension(7) :: US_StdAtm_znodes
      real(kind=sp),dimension(7) :: US_StdAtm_pnodes

      real(kind=sp) :: frac
      real(kind=sp) :: Delta_z
      integer :: k,kk

      US_StdAtm_znodes = (/ 0.0_sp, 11.0_sp, 20.0_sp, 32.0_sp, &
                           47.0_sp, 51.0_sp, 71.0_sp/)
      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 0.669389_sp, 0.0395642_sp/)

      IF(pin.ge.US_StdAtm_pnodes(1))THEN
        MR_Z_US_StdAtm = US_StdAtm_pnodes(1)
      ELSEIf(pin.le.US_StdAtm_pnodes(7))THEN
        MR_Z_US_StdAtm = US_StdAtm_pnodes(7)
        write(*,*)"MR ERROR: Pressure provided is below lowest value of US StdAtm."
        write(*,*)"Stopping in Z_US_StdAtm"
        stop 1
      ELSE
        ! interpolate
        ! Start from the top since we assume the requested point is near the top
        kk = 1
        DO k = 6,1,-1
          If(pin.le.US_StdAtm_pnodes(k).and.pin.gt.US_StdAtm_pnodes(k+1)) kk = k
        ENDDO
        frac = (pin-US_StdAtm_pnodes(kk)) / &
               (US_StdAtm_pnodes(kk+1)-US_StdAtm_pnodes(kk))
        Delta_z = US_StdAtm_znodes(kk+1)-US_StdAtm_znodes(kk)
        ! Note: this should really be an exponential fit, not linear
        MR_Z_US_StdAtm = US_StdAtm_znodes(kk) + Delta_z * frac
      ENDIF

      return

      end function MR_Z_US_StdAtm


!##############################################################################
!
!  Returns the pressure in hPa given the height in km
!
!  From http://en.wikipedia.org/wiki/U.S._Standard_Atmosphere
!
!##############################################################################

      function MR_Pres_US_StdAtm(zin)

      implicit none

      real(kind=sp) :: MR_Pres_US_StdAtm
      real(kind=sp) :: zin

      real(kind=sp),dimension(7) :: US_StdAtm_znodes
      real(kind=sp),dimension(7) :: US_StdAtm_pnodes

      real(kind=sp) :: pres0 = 1013.25_sp
      real(kind=sp) :: skinz   = 7.0

      US_StdAtm_znodes = (/ 0.0_sp, 11.0_sp, 20.0_sp, 32.0_sp, &
                           47.0_sp, 51.0_sp, 71.0_sp/)
      US_StdAtm_pnodes = (/1013.25_sp, 226.321_sp, 54.7489_sp, 8.68019_sp, &
                           1.10906_sp, 0.669389_sp, 0.0395642_sp/)

        ! Eq 1.8 of Wallace and Hobbs
      MR_Pres_US_StdAtm = pres0 * exp(-zin/skinz)

      return

      end function MR_Pres_US_StdAtm

!##############################################################################

!##############################################################################
!
!    MR_QC_3dvar
!
!##############################################################################

        subroutine MR_QC_3dvar(              &
                   nx_max,ny_max,nz1_max, &
                   z_array_sp,               &
                   nz2_max,               &
                   dum_array_sp,             &
                   fill_val_sp,              &
                   bc_low_sp, bc_high_sp)


        implicit none

        integer, parameter :: sp        = 4 ! single precision
        integer, parameter :: dp        = 8 ! double precision

        integer,intent(in) :: nx_max,ny_max,nz1_max
        real(kind=sp),intent(in)    :: z_array_sp(nx_max,ny_max,nz1_max)
        integer,intent(in) :: nz2_max
        real(kind=sp),intent(inout) :: dum_array_sp(nx_max,ny_max,nz2_max)
        real(kind=sp)          :: fill_val_sp
        real(kind=sp),optional :: bc_low_sp
        real(kind=sp),optional :: bc_high_sp

        logical,dimension(nz1_max) :: IsFillValue

        integer ::  i,j,k,kk,kkk,klow,khigh

        DO i=1,nx_max
          DO j=1,ny_max
            ! find all fill values in the column
            IsFillValue = .false.
            DO k=1,nz2_max
              IF(isnan(dum_array_sp(i,j,k)).or.      &  ! Some windfiles use NaN's for fill values
                 dum_array_sp(i,j,k).eq.fill_val_sp) &  ! Others have a specific number for fill
                     IsFillValue(k) = .true.
            ENDDO

            ! Set lower BC if requested and needed
            IF(IsFillValue(1).and.present(bc_low_sp)) &
                dum_array_sp(i,j,1) = bc_low_sp

            IF(nz2_max.gt.1)THEN
              ! Set upper BC if requested and needed
              IF(IsFillValue(nz2_max).and.present(bc_high_sp)) &
                  dum_array_sp(i,j,nz2_max) = bc_high_sp

              ! Now find the lowest non-Fill value
              DO klow=1,nz2_max
                If(.not.IsFillValue(klow))exit
              enddo
              ! And the highest non-Fill value
              DO khigh=nz2_max,1,-1
                If(.not.IsFillValue(khigh))exit
              enddo

              ! Set bottom values to the lowest real number
              !IF(klow.gt.1.and..not.present(bc_low_sp))THEN
              IF(klow.gt.1)THEN
                DO k = 1,klow-1
                  dum_array_sp(i,j,k) = dum_array_sp(i,j,klow)
                ENDDO
              ENDIF

              ! Set top values to the highest real number
              !IF(khigh.lt.nz2_max.and..not.present(bc_high_sp))THEN
              IF(khigh.lt.nz2_max)THEN
                DO k = nz2_max,khigh+1,-1
                  dum_array_sp(i,j,k) = dum_array_sp(i,j,khigh)
                ENDDO
              ENDIF

              ! Now check if there are any intermediate FillValues and
              ! interpolate
              DO k=klow+1,khigh-1
                IF(IsFillValue(k))THEN
                   ! linearly interpolate in z
                     ! Find first number bove
                   do kk = k+1,khigh,1
                     IF(.not.IsFillValue(kk))exit
                   enddo
                     ! Find first number beneath
                   do kkk = max(k-1,1),klow,-1
                     IF(.not.IsFillValue(kkk))exit
                   enddo
                   ! Note:  this assumes that nz2_max.le.nz1_max and that the
                   ! pressure values for dum_array_sp(1:nz2_max) are the same as
                   ! the pressure for z_array_sp(nz1_max)
                   dum_array_sp(i,j,k) = dum_array_sp(i,j,kk) + &
                         (dum_array_sp(i,j,kk)-dum_array_sp(i,j,kkk)) * &
                         (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))/ &
                         (z_array_sp(i,j,kk)-z_array_sp(i,j,kkk))
                ENDIF

              ENDDO
            ENDIF
          ENDDO
        ENDDO

        end subroutine MR_QC_3dvar

!##############################################################################

      end module MetReader

