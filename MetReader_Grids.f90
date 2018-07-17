!##############################################################################
!
!     MR_Set_Met_NCEPGeoGrid
!
!     This subroutine reads the variable and dimension IDs, fills the
!     coordinate dimension variables, and calculates xLL_meso, xUR_meso,
!     yLL_meso, and yUR_meso.
!
!     Allocated the dummy arrays for storing met data on met and computational
!     grids:  dum2d_met(nx_submet,ny_submet)
!             dum3d_metP(nx_submet,ny_submet,np_fullmet)
!             dum3d_metH(nx_submet,ny_submet,nz_comp)
!             dum2d_comp(nx_comp,ny_comp)
!             dum3d_compH(nx_comp,ny_comp,nz_comp)
!
!     Note: NCEP grids described here:
!              http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
!
!##############################################################################

      subroutine MR_Set_Met_NCEPGeoGrid(igrid)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: igrid

      integer :: i
      real(kind=dp) :: x_start,y_start
      real(kind=sp),dimension(:)  ,allocatable :: tmp_sp      !

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                          MR_Set_Met_NCEPGeoGrid            ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(igrid.eq.1051)then
        ! Not an NCEP grid
        !  This grid is for the SENAMHI 22 km files
        ! proj +proj=merc  +lat_ts=56.792 +lon_0=274.784 +R=6367.470
        ! 198.475 18.073
        !   0.00    1920.62
        ! 206.131 23.088
        !   800.00  2480.60
        ! 0 5 274.784 56.792 0.933 6367.470 #Proj flags and params

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 5
        Met_lam0          = 274.784027099609_8
        Met_phi0          =  56.7920036315918_8
        Met_k0            =  0.933_8
        Met_Re            =  6367.470_8

        nx_fullmet = 93
        ny_fullmet = 112
        dx_met_const = 22.0_sp
        dy_met_const = 22.0_sp
        !dx_met_const = 11.0_sp
        !dy_met_const = 11.0_sp
        x_start = 0.0_sp
        y_start = -1232.906_sp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1041)then
         ! Not an NCEP grid
         !  This grid is for the NASA Np
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 1152
        ny_fullmet = 721
        dx_met_const = 0.3125_sp
        dy_met_const = 0.25_sp
        x_start = -180.0_dp    ! These are double-precision for GEOS_Np
        y_start = -90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1040.or.igrid.eq.1024)then
         ! Not an NCEP grid
         !  This grid is for the NASA GEOS-5 Cp or MERRA-2

        !Met_dim_names(3) = "lat"        ! y        (-90.0 -> 90.0) 361
        !  Met_dim_IsAvailable(3)=.true.
        !Met_dim_names(4) = "lon"        ! x        (-180.0 -> 179.375) 576

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 576
        ny_fullmet = 361
        dx_met_const = 0.625_sp
        dy_met_const = 0.5_sp
        x_start = -180.0_dp
        y_start = -90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1033)then
         ! Not an NCEP grid
         !  This grid is for the CAM files
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 96
        ny_fullmet = 48
        dx_met_const = 3.75_sp
        dy_met_const = 3.7_sp
        x_start = 0.0_dp
        y_start = -87.1590945558628_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1032)then
         ! Not an NCEP grid
         !  This grid is for the AFWA files
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 1440
        ny_fullmet = 721
        dx_met_const = 0.25_sp
        dy_met_const = 0.25_sp
        x_start = 0.0_dp
        y_start = -90.0_dp
        !dx_met_const = 0.234374_sp
        !dy_met_const = 0.15625_sp
        !x_start = 0.117187_dp
        !y_start = -89.921875_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1031)then
         ! Not an NCEP grid
         !  This grid is for the Catania files
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        nx_fullmet = 94
        ny_fullmet = 98
        dx_met_const = 0.06_sp
        dy_met_const = 0.06_sp
        x_start = 12.5_dp
        y_start = 34.5_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      !elseif(igrid.eq.1024)then
      !   ! Not an NCEP grid
      !   !  This grid is for the NASA MERRA files
      !  IsLatLon_MetGrid  = .true.
      !  IsGlobal_MetGrid  = .true.
      !  IsRegular_MetGrid = .true.
      !  nx_fullmet = 288
      !  ny_fullmet = 144
      !  dx_met_const = 1.25_sp
      !  dy_met_const = 1.25_sp
      !  x_start = -179.375_dp
      !  y_start = 89.375_dp
      !  allocate(x_fullmet_sp(0:nx_fullmet+1))
      !  allocate(y_fullmet_sp(ny_fullmet))
      !  allocate(MR_dx_met(nx_fullmet))
      !  allocate(MR_dy_met(ny_fullmet))
      !  do i = 0,nx_fullmet+1
      !    x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
      !  enddo
      !  do i = 1,ny_fullmet
      !    y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
      !  enddo
      !  do i = 1,nx_fullmet
      !    MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
      !  enddo
      !  do i = 1,ny_fullmet-1
      !    MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
      !  enddo
      !  MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.1027)then
         ! Not an NCEP grid
         !  This grid is for the NOAA Reanalysis
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 180
        ny_fullmet = 91
        dx_met_const = 2.0_sp
        dy_met_const = 2.0_sp
        x_start =  0.0_dp
        y_start = 90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.2)then
         ! Used by NCEP DOE reanalysis, NCEP-1
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 144
        ny_fullmet = 73
        dx_met_const = 2.5_sp
        dy_met_const = 2.5_sp
        x_start =  0.0_dp
        y_start = 90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.3)then
         ! Used by GFS forecast
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 360
        ny_fullmet = 181
        dx_met_const = 1.0_sp
        dy_met_const = 1.0_sp
        x_start =  0.0_dp
        y_start = 90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.4)then
         ! Used by GFS forecast
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 720
        ny_fullmet = 361
        dx_met_const = 0.5_sp
        dy_met_const = 0.5_sp
        x_start =  0.0_dp
        y_start = 90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.91)then
        ! NAM 3-km Polar Sterographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID91
        !   1822145-point (1649x1105) N. Hemisphere Polar Stereographic grid
        !   oriented 150W; 
        !        Polar_Stereographic:grid_mapping_name = "polar_stereographic" ;
        !        Polar_Stereographic:longitude_of_projection_origin = 210. ;
        !        Polar_Stereographic:straight_vertical_longitude_from_pole = 210.;
        !        Polar_Stereographic:scale_factor_at_projection_origin = 0.933 ;
        !        Polar_Stereographic:latitude_of_projection_origin = 90. ;
        !        Polar_Stereographic:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_Dx = 2976.0  ;
        !        Polar_Stereographic:GRIB_param_Dy = 2976.0 ;
        !        Polar_Stereographic:GRIB_param_GDSkey = -1633822368 ;
        !        Polar_Stereographic:GRIB_param_La1 = 40.53 ;
        !        Polar_Stereographic:GRIB_param_LaD = 60. ;
        !        Polar_Stereographic:GRIB_param_Lo1 = 181.42899 ;
        !        Polar_Stereographic:GRIB_param_LoV = 210. ;
        !        Polar_Stereographic:GRIB_param_NpProj = "true" ;
        !        Polar_Stereographic:GRIB_param_Nx = 1649 ;
        !        Polar_Stereographic:GRIB_param_Ny = 1105 ;
        !        Polar_Stereographic:GRIB_param_ProjFlag = 0 ;
        !        Polar_Stereographic:GRIB_param_Quasi = "false" ;
        !        Polar_Stereographic:GRIB_param_ResCompFlag = 8 ;
        !        Polar_Stereographic:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Polar_Stereographic:GRIB_param_Winds = "Relative" ;
        !        Polar_Stereographic:GRIB_param_grid_name = "Polar_Stereographic" ;
        !        Polar_Stereographic:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Polar_Stereographic:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_grid_shape_code = 6 ;
        !        Polar_Stereographic:GRIB_param_grid_type = 20 ;
        !        Polar_Stereographic:GRIB_param_grid_units = "m" ;
        ! proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229
        !   181.42899 40.5301
        !      -2619.36159134661 -4810.03724324973
        !   266.3082 63.9757
        !       2285.91081099714 -1523.98097371848
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -150.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -150.0 90.0 0.933 6371.229    #Proj flags and params
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 1
        Met_lam0          = -150.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 1649
        ny_fullmet = 1105
        dx_met_const = 2.976000_sp
        dy_met_const = 2.976000_sp
        x_start = -2619.397217_dp
        y_start = -4810.102539_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)


      elseif(igrid.eq.104)then
        ! NAM 90-km Polar Sterographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID104
        !   16170-point (147x110) N. Hemisphere Polar Stereographic grid oriented
        !   105W; pole at (75.5,109.5). (NGM Super C grid)
        !   90.75464 km at 60N
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 255. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=255  +lat_0=90 +k_0=0.933 +R=6371.229
        !
        ! -6761.21 -9846.821
        ! 
        !  6489.02 45.47379
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -105.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -105.0 90.0 0.933 6371.229    #Proj flags and params
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 1
        Met_lam0          = -105.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 147
        ny_fullmet = 110
        dx_met_const = 90.75500_sp
        dy_met_const = 90.75500_sp
        x_start = -6761.11795473085_dp
        y_start = -9846.6868574523_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.170)then
        ! Global Gaussian Lat/Lon T170
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID170
        ! This is used by the ERA data

        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .false.
        nx_fullmet = 512
        ny_fullmet = 256
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(tmp_sp(ny_fullmet))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        x_fullmet_sp(1:nx_fullmet) = &
         (/ 0.0_sp,0.7031252_sp, 1.406250_sp, 2.109376_sp, 2.812501_sp, 3.515626_sp, 4.218751_sp,& 
       4.921877_sp, 5.625002_sp, 6.328127_sp, 7.031252_sp, 7.734378_sp, 8.437503_sp, 9.140628_sp,& 
       9.843754_sp, 10.54688_sp, 11.25000_sp, 11.95313_sp, 12.65625_sp, 13.35938_sp, 14.06250_sp,& 
       14.76563_sp, 15.46876_sp, 16.17188_sp, 16.87501_sp, 17.57813_sp, 18.28126_sp, 18.98438_sp,& 
       19.68751_sp, 20.39063_sp, 21.09376_sp, 21.79688_sp, 22.50001_sp, 23.20313_sp, 23.90626_sp,& 
       24.60938_sp, 25.31251_sp, 26.01563_sp, 26.71876_sp, 27.42188_sp, 28.12501_sp, 28.82813_sp,& 
       29.53126_sp, 30.23439_sp, 30.93751_sp, 31.64064_sp, 32.34376_sp, 33.04689_sp, 33.75001_sp,& 
       34.45314_sp, 35.15626_sp, 35.85939_sp, 36.56251_sp, 37.26564_sp, 37.96876_sp, 38.67189_sp,& 
       39.37502_sp, 40.07814_sp, 40.78127_sp, 41.48439_sp, 42.18752_sp, 42.89064_sp, 43.59377_sp,& 
       44.29689_sp, 45.00002_sp, 45.70314_sp, 46.40627_sp, 47.10939_sp, 47.81252_sp, 48.51564_sp,& 
       49.21877_sp, 49.92189_sp, 50.62502_sp, 51.32814_sp, 52.03127_sp, 52.73439_sp, 53.43752_sp,& 
       54.14064_sp, 54.84377_sp, 55.54689_sp, 56.25002_sp, 56.95314_sp, 57.65627_sp, 58.35939_sp,& 
       59.06252_sp, 59.76564_sp, 60.46877_sp, 61.17190_sp, 61.87502_sp, 62.57815_sp, 63.28127_sp,& 
       63.98440_sp, 64.68752_sp, 65.39065_sp, 66.09377_sp, 66.79690_sp, 67.50002_sp, 68.20315_sp,& 
       68.90627_sp, 69.60940_sp, 70.31252_sp, 71.01565_sp, 71.71877_sp, 72.42190_sp, 73.12502_sp,& 
       73.82815_sp, 74.53127_sp, 75.23440_sp, 75.93752_sp, 76.64065_sp, 77.34378_sp, 78.04691_sp,& 
       78.75003_sp, 79.45316_sp, 80.15628_sp, 80.85941_sp, 81.56253_sp, 82.26566_sp, 82.96878_sp,& 
       83.67191_sp, 84.37503_sp, 85.07816_sp, 85.78128_sp, 86.48441_sp, 87.18753_sp, 87.89066_sp,& 
       88.59378_sp, 89.29691_sp, 90.00003_sp, 90.70316_sp, 91.40628_sp, 92.10941_sp, 92.81253_sp,& 
       93.51566_sp, 94.21878_sp, 94.92191_sp, 95.62503_sp, 96.32816_sp, 97.03128_sp, 97.73441_sp,& 
       98.43753_sp, 99.14066_sp, 99.84379_sp, 100.5469_sp, 101.2500_sp, 101.9532_sp, 102.6563_sp,& 
       103.3594_sp, 104.0625_sp, 104.7657_sp, 105.4688_sp, 106.1719_sp, 106.8750_sp, 107.5782_sp,& 
       108.2813_sp, 108.9844_sp, 109.6875_sp, 110.3907_sp, 111.0938_sp, 111.7969_sp, 112.5000_sp,& 
       113.2032_sp, 113.9063_sp, 114.6094_sp, 115.3125_sp, 116.0157_sp, 116.7188_sp, 117.4219_sp,& 
       118.1250_sp, 118.8282_sp, 119.5313_sp, 120.2344_sp, 120.9375_sp, 121.6407_sp, 122.3438_sp,& 
       123.0469_sp, 123.7500_sp, 124.4532_sp, 125.1563_sp, 125.8594_sp, 126.5625_sp, 127.2657_sp,& 
       127.9688_sp, 128.6719_sp, 129.3750_sp, 130.0782_sp, 130.7813_sp, 131.4844_sp, 132.1875_sp,& 
       132.8907_sp, 133.5938_sp, 134.2969_sp, 135.0000_sp, 135.7032_sp, 136.4063_sp, 137.1094_sp,& 
       137.8125_sp, 138.5157_sp, 139.2188_sp, 139.9219_sp, 140.6250_sp, 141.3282_sp, 142.0313_sp,& 
       142.7344_sp, 143.4375_sp, 144.1407_sp, 144.8438_sp, 145.5469_sp, 146.2500_sp, 146.9532_sp,& 
       147.6563_sp, 148.3594_sp, 149.0625_sp, 149.7657_sp, 150.4688_sp, 151.1719_sp, 151.8750_sp,& 
       152.5782_sp, 153.2813_sp, 153.9844_sp, 154.6876_sp, 155.3907_sp, 156.0938_sp, 156.7969_sp,& 
       157.5001_sp, 158.2032_sp, 158.9063_sp, 159.6094_sp, 160.3126_sp, 161.0157_sp, 161.7188_sp,& 
       162.4219_sp, 163.1251_sp, 163.8282_sp, 164.5313_sp, 165.2344_sp, 165.9376_sp, 166.6407_sp,& 
       167.3438_sp, 168.0469_sp, 168.7501_sp, 169.4532_sp, 170.1563_sp, 170.8594_sp, 171.5626_sp,& 
       172.2657_sp, 172.9688_sp, 173.6719_sp, 174.3751_sp, 175.0782_sp, 175.7813_sp, 176.4844_sp,& 
       177.1876_sp, 177.8907_sp, 178.5938_sp, 179.2969_sp, 180.0001_sp, 180.7032_sp, 181.4063_sp,& 
       182.1094_sp, 182.8126_sp, 183.5157_sp, 184.2188_sp, 184.9219_sp, 185.6251_sp, 186.3282_sp,& 
       187.0313_sp, 187.7344_sp, 188.4376_sp, 189.1407_sp, 189.8438_sp, 190.5469_sp, 191.2501_sp,& 
       191.9532_sp, 192.6563_sp, 193.3594_sp, 194.0626_sp, 194.7657_sp, 195.4688_sp, 196.1719_sp,& 
       196.8751_sp, 197.5782_sp, 198.2813_sp, 198.9845_sp, 199.6876_sp, 200.3907_sp, 201.0938_sp,& 
       201.7970_sp, 202.5001_sp, 203.2032_sp, 203.9063_sp, 204.6095_sp, 205.3126_sp, 206.0157_sp,& 
       206.7188_sp, 207.4220_sp, 208.1251_sp, 208.8282_sp, 209.5313_sp, 210.2345_sp, 210.9376_sp,& 
       211.6407_sp, 212.3438_sp, 213.0470_sp, 213.7501_sp, 214.4532_sp, 215.1563_sp, 215.8595_sp,& 
       216.5626_sp, 217.2657_sp, 217.9688_sp, 218.6720_sp, 219.3751_sp, 220.0782_sp, 220.7813_sp,& 
       221.4845_sp, 222.1876_sp, 222.8907_sp, 223.5938_sp, 224.2970_sp, 225.0001_sp, 225.7032_sp,& 
       226.4063_sp, 227.1095_sp, 227.8126_sp, 228.5157_sp, 229.2188_sp, 229.9220_sp, 230.6251_sp,& 
       231.3282_sp, 232.0313_sp, 232.7345_sp, 233.4376_sp, 234.1407_sp, 234.8438_sp, 235.5470_sp,& 
       236.2501_sp, 236.9532_sp, 237.6563_sp, 238.3595_sp, 239.0626_sp, 239.7657_sp, 240.4688_sp,& 
       241.1720_sp, 241.8751_sp, 242.5782_sp, 243.2813_sp, 243.9845_sp, 244.6876_sp, 245.3907_sp,& 
       246.0938_sp, 246.7970_sp, 247.5001_sp, 248.2032_sp, 248.9063_sp, 249.6095_sp, 250.3126_sp,& 
       251.0157_sp, 251.7188_sp, 252.4220_sp, 253.1251_sp, 253.8282_sp, 254.5313_sp, 255.2345_sp,& 
       255.9376_sp, 256.6407_sp, 257.3438_sp, 258.0470_sp, 258.7501_sp, 259.4532_sp, 260.1563_sp,& 
       260.8595_sp, 261.5626_sp, 262.2657_sp, 262.9688_sp, 263.6720_sp, 264.3751_sp, 265.0782_sp,& 
       265.7813_sp, 266.4845_sp, 267.1876_sp, 267.8907_sp, 268.5938_sp, 269.2970_sp, 270.0001_sp,& 
       270.7032_sp, 271.4063_sp, 272.1095_sp, 272.8126_sp, 273.5157_sp, 274.2188_sp, 274.9220_sp,& 
       275.6251_sp, 276.3282_sp, 277.0313_sp, 277.7345_sp, 278.4376_sp, 279.1407_sp, 279.8438_sp,& 
       280.5470_sp, 281.2501_sp, 281.9532_sp, 282.6563_sp, 283.3595_sp, 284.0626_sp, 284.7657_sp,& 
       285.4688_sp, 286.1720_sp, 286.8751_sp, 287.5782_sp, 288.2813_sp, 288.9845_sp, 289.6876_sp,& 
       290.3907_sp, 291.0938_sp, 291.7970_sp, 292.5001_sp, 293.2032_sp, 293.9063_sp, 294.6095_sp,& 
       295.3126_sp, 296.0157_sp, 296.7188_sp, 297.4220_sp, 298.1251_sp, 298.8282_sp, 299.5313_sp,& 
       300.2345_sp, 300.9376_sp, 301.6407_sp, 302.3438_sp, 303.0470_sp, 303.7501_sp, 304.4532_sp,& 
       305.1563_sp, 305.8595_sp, 306.5626_sp, 307.2657_sp, 307.9689_sp, 308.6720_sp, 309.3751_sp,& 
       310.0782_sp, 310.7814_sp, 311.4845_sp, 312.1876_sp, 312.8907_sp, 313.5939_sp, 314.2970_sp,& 
       315.0001_sp, 315.7032_sp, 316.4064_sp, 317.1095_sp, 317.8126_sp, 318.5157_sp, 319.2189_sp,& 
       319.9220_sp, 320.6251_sp, 321.3282_sp, 322.0314_sp, 322.7345_sp, 323.4376_sp, 324.1407_sp,& 
       324.8439_sp, 325.5470_sp, 326.2501_sp, 326.9532_sp, 327.6564_sp, 328.3595_sp, 329.0626_sp,& 
       329.7657_sp, 330.4689_sp, 331.1720_sp, 331.8751_sp, 332.5782_sp, 333.2814_sp, 333.9845_sp,& 
       334.6876_sp, 335.3907_sp, 336.0939_sp, 336.7970_sp, 337.5001_sp, 338.2032_sp, 338.9064_sp,& 
       339.6095_sp, 340.3126_sp, 341.0157_sp, 341.7189_sp, 342.4220_sp, 343.1251_sp, 343.8282_sp,& 
       344.5314_sp, 345.2345_sp, 345.9376_sp, 346.6407_sp, 347.3439_sp, 348.0470_sp, 348.7501_sp,& 
       349.4532_sp, 350.1564_sp, 350.8595_sp, 351.5626_sp, 352.2657_sp, 352.9689_sp, 353.6720_sp,& 
       354.3751_sp, 355.0782_sp, 355.7814_sp, 356.4845_sp, 357.1876_sp, 357.8907_sp, 358.5939_sp,& 
       359.2970_sp /)
        x_fullmet_sp(0)            = x_fullmet_sp(nx_fullmet) - 360.0_sp
        x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(         1) + 360.0_sp

        y_fullmet_sp(1:ny_fullmet) = &
                  (/ 89.46282_sp,  88.76695_sp,  88.06697_sp,   87.36607_sp,  86.66480_sp,  85.96337_sp, & 
       85.26185_sp,  84.56026_sp,  83.85863_sp,   83.15699_sp,  82.45532_sp,  81.75363_sp,  81.05194_sp, & 
       80.35023_sp,  79.64853_sp,  78.94681_sp,   78.24509_sp,  77.54337_sp,  76.84164_sp,  76.13991_sp, & 
       75.43818_sp,  74.73644_sp,  74.03471_sp,   73.33297_sp,  72.63123_sp,  71.92949_sp,  71.22775_sp, & 
       70.52601_sp,  69.82426_sp,  69.12252_sp,   68.42078_sp,  67.71903_sp,  67.01729_sp,  66.31554_sp, & 
       65.61379_sp,  64.91205_sp,  64.21030_sp,   63.50855_sp,  62.80680_sp,  62.10506_sp,  61.40331_sp, & 
       60.70156_sp,  59.99981_sp,  59.29806_sp,   58.59631_sp,  57.89456_sp,  57.19281_sp,  56.49106_sp, & 
       55.78931_sp,  55.08756_sp,  54.38581_sp,   53.68406_sp,  52.98231_sp,  52.28056_sp,  51.57881_sp, & 
       50.87706_sp,  50.17531_sp,  49.47356_sp,   48.77180_sp,  48.07005_sp,  47.36830_sp,  46.66655_sp, & 
       45.96480_sp,  45.26305_sp,  44.56129_sp,   43.85954_sp,  43.15779_sp,  42.45604_sp,  41.75429_sp, & 
       41.05254_sp,  40.35078_sp,  39.64903_sp,   38.94728_sp,  38.24553_sp,  37.54378_sp,  36.84202_sp, & 
       36.14027_sp,  35.43852_sp,  34.73677_sp,   34.03502_sp,  33.33326_sp,  32.63151_sp,  31.92976_sp, & 
       31.22800_sp,  30.52625_sp,  29.82450_sp,   29.12275_sp,  28.42099_sp,  27.71924_sp,  27.01749_sp, & 
       26.31573_sp,  25.61398_sp,  24.91223_sp,   24.21048_sp,  23.50872_sp,  22.80697_sp,  22.10522_sp, & 
       21.40347_sp,  20.70171_sp,  19.99996_sp,   19.29821_sp,  18.59645_sp,  17.89470_sp,  17.19295_sp, & 
       16.49120_sp,  15.78944_sp,  15.08769_sp,   14.38594_sp,  13.68418_sp,  12.98243_sp,  12.28068_sp, & 
       11.57893_sp,  10.87717_sp,  10.17542_sp,   9.473666_sp,  8.771913_sp,  8.070160_sp,  7.368407_sp, & 
       6.666654_sp,  5.964901_sp,  5.263148_sp,   4.561395_sp,  3.859642_sp,  3.157889_sp,  2.456136_sp, & 
       1.754383_sp,  1.052630_sp, 0.3508765_sp, -0.3508765_sp,  -1.05263_sp, -1.754383_sp, -2.456136_sp, & 
      -3.157889_sp, -3.859642_sp, -4.561395_sp,  -5.263148_sp, -5.964901_sp, -6.666654_sp, &
      -7.368407_sp, -8.070160_sp, -8.771913_sp,  -9.473666_sp, -10.17542_sp, -10.87717_sp, &
      -11.57893_sp, -12.28068_sp, -12.98243_sp,  -13.68418_sp, -14.38594_sp, -15.08769_sp, &
      -15.78944_sp, -16.49120_sp, -17.19295_sp,  -17.89470_sp, -18.59645_sp, -19.29821_sp, &
      -19.99996_sp, -20.70171_sp, -21.40347_sp,  -22.10522_sp, -22.80697_sp, -23.50872_sp, &
      -24.21048_sp, -24.91223_sp, -25.61398_sp,  -26.31573_sp, -27.01749_sp, -27.71924_sp, &
      -28.42099_sp, -29.12275_sp, -29.82450_sp,  -30.52625_sp, -31.22800_sp, -31.92976_sp, -32.63151_sp, & 
      -33.33326_sp, -34.03502_sp, -34.73677_sp,  -35.43852_sp, -36.14027_sp, -36.84202_sp, &
      -37.54378_sp, -38.24553_sp, -38.94728_sp,  -39.64903_sp, -40.35078_sp, -41.05254_sp, &
      -41.75429_sp, -42.45604_sp, -43.15779_sp,  -43.85954_sp, -44.56129_sp, -45.26305_sp, &
      -45.96480_sp, -46.66655_sp, -47.36830_sp,  -48.07005_sp, -48.77180_sp, -49.47356_sp, -50.17531_sp, & 
      -50.87706_sp, -51.57881_sp, -52.28056_sp,  -52.98231_sp, -53.68406_sp, -54.38581_sp, &
      -55.08756_sp, -55.78931_sp, -56.49106_sp,  -57.19281_sp, -57.89456_sp, -58.59631_sp, &
      -59.29806_sp, -59.99981_sp, -60.70156_sp,  -61.40331_sp, -62.10506_sp, -62.80680_sp, &
      -63.50855_sp, -64.21030_sp, -64.91205_sp,  -65.61379_sp, -66.31554_sp, -67.01729_sp, &
      -67.71903_sp, -68.42078_sp, -69.12252_sp,  -69.82426_sp, -70.52601_sp, -71.22775_sp, &
      -71.92949_sp, -72.63123_sp, -73.33297_sp,  -74.03471_sp, -74.73644_sp, -75.43818_sp, &
      -76.13991_sp, -76.84164_sp, -77.54337_sp,  -78.24509_sp, -78.94681_sp, -79.64853_sp, &
      -80.35023_sp, -81.05194_sp, -81.75363_sp,  -82.45532_sp, -83.15699_sp, -83.85863_sp, &
      -84.56026_sp, -85.26185_sp, -85.96337_sp,  -86.66480_sp, -87.36607_sp, -88.06697_sp, &
      -88.76695_sp, -89.46282_sp /)
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i)-y_fullmet_sp(i+1)
        enddo
        deallocate(tmp_sp)
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.193)then
         ! Used by GFS forecast (0.25)
        IsLatLon_MetGrid  = .true.
        IsGlobal_MetGrid  = .true.
        IsRegular_MetGrid = .true.
        nx_fullmet = 1440
        ny_fullmet = 721
        dx_met_const = 0.25_sp
        dy_met_const = 0.25_sp
        x_start =  0.0_dp
        y_start = 90.0_dp
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.196)then
        ! HI 2.5-km Mercator
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID196
        !   72225-point (321x225) Mercator
        !        Mercator:grid_mapping_name = "mercator" ;
        !        Mercator:standard_parallel = 20. ;
        !        Mercator:longitude_of_projection_origin = 198.475006103516 ;
        !        Mercator:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Mercator:GRIB_param_BasicAngle = 0 ;
        !        Mercator:GRIB_param_Dx = 2500. ;
        !        Mercator:GRIB_param_Dy = 2500. ;
        !        Mercator:GRIB_param_GDSkey = -1286480248 ;
        !        Mercator:GRIB_param_La1 = 18.073 ;
        !        Mercator:GRIB_param_La2 = 23.088 ;
        !        Mercator:GRIB_param_LaD = 20. ;
        !        Mercator:GRIB_param_Lo1 = 198.475 ;
        !        Mercator:GRIB_param_Lo2 = 206.13101 ;
        !        Mercator:GRIB_param_Nx = 321 ;
        !        Mercator:GRIB_param_Ny = 225 ;
        !        Mercator:GRIB_param_Quasi = "false" ;
        !        Mercator:GRIB_param_ResCompFlag = 56 ;
        !        Mercator:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Mercator:GRIB_param_Winds = "Relative" ;
        !        Mercator:GRIB_param_grid_name = "Mercator" ;
        !        Mercator:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Mercator:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Mercator:GRIB_param_grid_shape_code = 6 ;
        !        Mercator:GRIB_param_grid_type = 10 ;
        !        Mercator:GRIB_param_grid_units = "m" ;
        ! proj +proj=merc  +lat_ts=20.0 +lon_0=198.475 +R=6371.229
        ! 198.475 18.073
        !   0.00    1920.62
        ! 206.131 23.088
        !   800.00  2480.60
        ! 0 5 198.475 20.0 0.933 6371.229    #Proj flags and params
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 5
        Met_lam0          = 198.475_8
        Met_phi0          =  20.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 321
        ny_fullmet = 225
        dx_met_const = 2.5_sp
        dy_met_const = 2.5_sp
        x_start = 0.0_dp
        y_start = 1920.618_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.198)then
        ! NAM 6-km Polar Sterographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID198
        !   456225-point (825x553) N. Hemisphere Polar Stereographic grid
        !   oriented 150W; 
        !        Polar_Stereographic:grid_mapping_name = "polar_stereographic" ;
        !        Polar_Stereographic:longitude_of_projection_origin = 210. ;
        !        Polar_Stereographic:straight_vertical_longitude_from_pole = 210.;
        !        Polar_Stereographic:scale_factor_at_projection_origin = 0.933 ;
        !        Polar_Stereographic:latitude_of_projection_origin = 90. ;
        !        Polar_Stereographic:earth_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_Dx = 5953.0005 ;
        !        Polar_Stereographic:GRIB_param_Dy = 5953.0005 ;
        !        Polar_Stereographic:GRIB_param_GDSkey = -1633822368 ;
        !        Polar_Stereographic:GRIB_param_La1 = 40.53 ;
        !        Polar_Stereographic:GRIB_param_LaD = 60. ;
        !        Polar_Stereographic:GRIB_param_Lo1 = 181.42899 ;
        !        Polar_Stereographic:GRIB_param_LoV = 210. ;
        !        Polar_Stereographic:GRIB_param_NpProj = "true" ;
        !        Polar_Stereographic:GRIB_param_Nx = 825 ;
        !        Polar_Stereographic:GRIB_param_Ny = 553 ;
        !        Polar_Stereographic:GRIB_param_ProjFlag = 0 ;
        !        Polar_Stereographic:GRIB_param_Quasi = "false" ;
        !        Polar_Stereographic:GRIB_param_ResCompFlag = 8 ;
        !        Polar_Stereographic:GRIB_param_VectorComponentFlag = "gridRelative" ;
        !        Polar_Stereographic:GRIB_param_Winds = "Relative" ;
        !        Polar_Stereographic:GRIB_param_grid_name = "Polar_Stereographic" ;
        !        Polar_Stereographic:GRIB_param_grid_radius_spherical_earth = 6371229. ;
        !        Polar_Stereographic:GRIB_param_grid_shape = "Earth spherical with radius of 6371229.0 m" ;
        !        Polar_Stereographic:GRIB_param_grid_shape_code = 6 ;
        !        Polar_Stereographic:GRIB_param_grid_type = 20 ;
        !        Polar_Stereographic:GRIB_param_grid_units = "m" ;
        ! proj +proj=stere  +lon_0=210  +lat_0=90 +k_0=0.933 +R=6371.229
        !   181.42899 40.5301
        !      -2619.36159134661 -4810.03724324973
        !   266.3082 63.9757
        !       2285.91081099714 -1523.98097371848
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -150.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -150.0 90.0 0.933 6371.229    #Proj flags and params
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 1
        Met_lam0          = -150.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 825
        ny_fullmet = 553
        dx_met_const = 5.953000_sp
        dy_met_const = 5.953000_sp
        x_start = -2619.397217_dp
        y_start = -4810.102539_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.212)then
        ! CONUS 40-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID212
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3250.731 4368.582
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params                               

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 185
        ny_fullmet = 129
        dx_met_const = 40.635_sp
        dy_met_const = 40.635_sp
        x_start =  -4226.10860264919_dp
        y_start =  -832.697842685899_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.215)then
        ! CONUS 20-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID215
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3250.916 4368.71
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params  
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 369
        ny_fullmet = 257
        dx_met_const = 20.317625_sp
        dy_met_const = 20.317625_sp
        x_start =  -4226.108_dp
        y_start =  -832.6978_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.216)then
        ! NAM 45-km Polar Sterographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID216
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 225. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229
        ! 187.0 30.0
        !   -4225.928 -5408.941
        ! 297.15 70.111
        !    1984.072 -638.9415
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -135.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -135.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 1
        Met_lam0          = -135.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 139
        ny_fullmet = 107
        dx_met_const = 45.00000_sp
        dy_met_const = 45.00000_sp
        x_start =  -4225.87071154953_dp
        y_start =  -5408.86785597764_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.218)then
        ! CONUS 12-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID218
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 25.;
        !        LambertConformal_Projection:longitude_of_central_meridian = 265.;
        !        LambertConformal_Projection:standard_parallel = 25. ;
        !        LambertConformal_Projection:earth_radius = 6371229. ;
        ! proj +proj=lcc +lon_0=265.0 +lat_0=25.0 +lat_1=25.0 +lat_2=25.0 +R=6371.229
        !   226.541 12.190
        !     -4226.108 -832.6978
        !   310.615 57.290
        !      3246.974 4372.859
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 265.0 25.0 25.0 25.0 6371.229    #Proj flags and params  
        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 4
        Met_lam0          =  265.0_8
        Met_phi0          =  25.0_8
        Met_phi1          =  25.0_8
        Met_phi2          =  25.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 614
        ny_fullmet = 428
        dx_met_const = 12.191000_sp
        dy_met_const = 12.191000_sp
        x_start =  -4226.108_dp
        y_start =  -832.6978_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.221)then
        ! NAM 32-km Lambert Conformal
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID221
        !        LambertConformal_Projection:grid_mapping_name = "lambert_conformal_conic" ;
        !        LambertConformal_Projection:latitude_of_projection_origin = 50. ;
        !        LambertConformal_Projection:longitude_of_central_meridian = 253. ;
        !        LambertConformal_Projection:standard_parallel = 50. ;
        !    NCEP FC        LambertConformal_Projection:earth_radius = 6371229. ;
        !    Reported NARR Reanal    Lambert_Conformal:GRIB_param_grid_radius_spherical_earth = 6367.47 ;
        ! proj +proj=lcc +lon_0=-107.0 +lat_0=50.0 +lat_1=50.0 +lat_2=50.0 +R=6371.229
        !   214.50  1.0
        !     -5632.668 -4612.566
        !   357.43 46.352
        !      5664.457 4347.222
        !  latlonflag  = 0         : projected grid
        !  projflag    = 4         : Lambert conformal conic
        !  lam0     = 265.0     : longitude of projection point
        !  phi0        =  25.0     : latitude of projection point
        !  phi1        =  25.0     : latitude of cone intersection 1
        !  phi2        =  25.0     : latitude of cone intersection 2
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 4 -107.0 50.0 50.0 50.0 6371.229    #Proj flags and params  

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 4
        Met_lam0          =  107.0_8
        Met_phi0          =  50.0_8
        Met_phi1          =  50.0_8
        Met_phi2          =  50.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 349
        ny_fullmet = 277
        dx_met_const = 32.463_sp
        dy_met_const = 32.463_sp
        x_start =  -5632.66705905226_dp
        y_start =  -4612.56764884087_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(igrid.eq.242)then
        ! NAM 11.25-km Polar Sterographic
        ! http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html#GRID242
        !        PolarStereographic_Projection:grid_mapping_name = "stereographic" ;
        !        PolarStereographic_Projection:longitude_of_projection_origin = 225. ;
        !        PolarStereographic_Projection:latitude_of_projection_origin = 90. ;
        !        PolarStereographic_Projection:scale_factor_at_projection_origin = 0.933012701892219 ;
        !        PolarStereographic_Projection:earth_radius = 6371229. ;
        ! proj +proj=stere  +lon_0=225  +lat_0=90 +k_0=0.933 +R=6371.229
        ! 187.0 30.0
        !   -4225.928 -5408.941
        ! 297.15 70.111
        !    1984.072 -638.9415
        !
        !  latlonflag  = 0         : projected grid
        !  projflag    = 1         : polar stereographic projection
        !  lam0     = -135.0    : longitude of projection point
        !  phi0        =  90.0     : latitude of projection point
        !  k0          =  0.933    : scale factor at projection point
        !  radius      =  6371.229 : earth radius for spherical earth
        ! 0 1 -135.0 90.0 0.933 6371.229    #Proj flags and params

        IsLatLon_MetGrid  = .false.
        IsGlobal_MetGrid  = .false.
        IsRegular_MetGrid = .true.
        Met_iprojflag     = 1
        Met_lam0          = -135.0_8
        Met_phi0          =  90.0_8
        Met_k0            =  0.933_8
        Met_Re            =  6371.229_8

        nx_fullmet = 553
        ny_fullmet = 425
        dx_met_const = 11.250000_sp
        dy_met_const = 11.250000_sp
        x_start =  -4225.87071154953_dp
        y_start =  -5408.86785597764_dp
        !allocate(x_fullmet_sp(nx_fullmet))
        allocate(x_fullmet_sp(0:nx_fullmet+1))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))
        do i = 0,nx_fullmet+1
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
        enddo
        do i = 1,nx_fullmet
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      else
        write(MR_global_info,*)"MR ERROR: MR_Set_Met_NCEPGeoGrid called with invalid code."
        stop 1
      endif
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_Met_NCEPGeoGrid

!##############################################################################
!
!     MR_Set_MetComp_Grids
!
!     This subroutine evaluates the full NWP grid for the subgrid needed.
!     If the NWP projection and the computational grid projection differ,
!     then the mapping from each computational grid point onto the NWP grid
!     is calculated through MR_Set_Comp2Met_Map.
!
!     Sets: 
!           n[t,x,y,p]_met     :: sets the size of the dimensions of the sub-met grid
!           [x,y,p]_met_sp     :: arrays holding dimension values of the sub-met grid
!           MR_dum2d_met(nx_submet,ny_submet)
!           MR_dum3d_metP(nx_submet,ny_submet,np_fullmet)
!           MR_dum3d_metH(nx_submet,ny_submet,nz_comp)
!           MR_dum2d_comp(nx_comp,ny_comp)
!           MR_dum3d_compH(nx_comp,ny_comp,nz_comp)
!           CompPoint_on_subMet_idx
!           bilin_map_wgt
!           CompPoint_X_on_Met_sp
!           CompPoint_Y_on_Met_sp
!
!##############################################################################

      subroutine MR_Set_MetComp_Grids

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: i, j
      integer :: ii,jj
      integer :: isubmet, jsubmet

      real(kind=sp) :: xLL,yLL
      real(kind=sp) :: xUR,yUR

      real(kind=sp) :: dely_sp,x_loc,y_loc,xc_sp,yc_sp,xfrac_sp,yfrac_sp
      integer :: ilat,ilon,ix1,ix2,iy1,iy2

      real(kind=dp) :: ptlon,ptlat,xin,yin,de_x,de_y

      real(kind=sp) :: xc,xfrac,yc,yfrac,px,py
      real(kind=sp) :: x_start_sub,y_start_sub

      real(kind=dp) :: xout,yout

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                 MR_Set_MetComp_Grids                       ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      call MR_Set_Comp2Met_Map

      ! Now calculate the indicies of the subgrid containing the
      ! computational grid.  Note, the subgrid of the wind file
      ! will be much coarser than the computational grid
      if(Map_Case.eq.1.or. & !  Both Comp Grid and Met grids are Lat/Lon
         Map_Case.eq.2)then  !  Both Comp Grid and Met grids are the same projection
        xLL = x_comp_sp(1)
        yLL = y_comp_sp(1)
        xUR = x_comp_sp(nx_comp)
        yUR = y_comp_sp(ny_comp)
      else
        write(MR_global_info,*)"Met and comp grids differ:"
        write(MR_global_info,2504)
        write(MR_global_info,2505)x_comp_sp(1),&
                     y_comp_sp(1),&
                     CompPoint_X_on_Met_sp(1,1),&
                     CompPoint_Y_on_Met_sp(1,1)
        write(MR_global_info,2505)x_comp_sp(nx_comp),&
                     y_comp_sp(1),&
                     CompPoint_X_on_Met_sp(nx_comp,1),&
                     CompPoint_Y_on_Met_sp(nx_comp,1)
        write(MR_global_info,2505)x_comp_sp(nx_comp),&
                     y_comp_sp(ny_comp),&
                     CompPoint_X_on_Met_sp(nx_comp,ny_comp),&
                     CompPoint_Y_on_Met_sp(nx_comp,ny_comp)
        write(MR_global_info,2505)x_comp_sp(1),&
                     y_comp_sp(ny_comp),&
                     CompPoint_X_on_Met_sp(1,ny_comp),&
                     CompPoint_Y_on_Met_sp(1,ny_comp)
        write(MR_global_info,*)" "

          ! This the branch for when Met and Comp grids differ
        xLL = minval(CompPoint_X_on_Met_sp(:,:))
        yLL = minval(CompPoint_Y_on_Met_sp(:,:))
        xUR = maxval(CompPoint_X_on_Met_sp(:,:))
        yUR = maxval(CompPoint_Y_on_Met_sp(:,:))
 2504   format(8x,'Comp grid corner',28x,'Met grid corner')
 2505   format(4x,'(',f10.4,',',f10.4,')',8x,'--->',8x,'(',f10.4,',',f10.4,')')
      endif

      if(IsLatLon_MetGrid)then
        if(xLL.gt.x_fullmet_sp(nx_fullmet).and.x_fullmet_sp(nx_fullmet).le.180.0_sp)then
          ! If the comp grid starts in the western himisphere (xLL>180) and if
          ! the global Met grid only extends to 180, then shift the comp grid
          ! into the domain of the met grid
          xLL=xLL-360.0_sp  ! This should only be true western hemisphere (xLL>180)
          xUR=xUR-360.0_sp  ! cases using MERRA
          x_comp_sp = x_comp_sp-360.0_sp
        endif
      endif

      write(MR_global_info,*)"Region of Met grid required by comp grid (in Met coordinates):"

      write(MR_global_info,2501)
      write(MR_global_info,2502)xLL,yUR,xUR,yUR
      write(MR_global_info,2503)
      write(MR_global_info,2503)
      write(MR_global_info,2502)xLL,yLL,xUR,yLL
      write(MR_global_info,2501)

 2501 format(4x,'----------------------------------------------------------------------')
 2502 format(4x,'| (',f10.4,',',f10.4,')',20x,'(',f10.4,',',f10.4,') |')
 2503 format(4x,'|                                                                    |')

      if(IsPeriodic_CompGrid)then
          ! If the domain is periodic, use the whole x-range of the wind file
          ! including a periodic mapping at either end
        nx_submet = nx_fullmet
        istart = 1
        iend = nx_fullmet
        write(MR_global_info,*) "Computational domain is periodic"
      else
        if(x_fullmet_sp(1).le.xLL)then
           ! Make sure the start of the comp grid is not below the domain of the
           ! met files
          if(IsRegular_MetGrid)then
            istart = max(1,floor((abs(x_fullmet_sp(1)-xLL)/dx_met_const)+1))
          else
            ! Need to actually march through all the values of x_fullmet_sp and find
            ! which interval is needed
            istart = -1
            do i = 2,nx_fullmet
              if(x_fullmet_sp(i).ge.xLL.and.x_fullmet_sp(i-1).lt.xLL)istart=i-1
            enddo
            if(istart.lt.0)then
              write(MR_global_info,*)"MR ERROR: Could not find iend"
              write(MR_global_info,*)"  ",xLL
              stop 1
            endif
          endif
        else
          write(MR_global_info,*)"MR ERROR: xLL < x_fullmet_sp(1)"
          write(MR_global_info,*)"     x_fullmet_sp(1) = ",x_fullmet_sp(1)
          write(MR_global_info,*)"     xLL             = ",xLL
          stop 1
        endif
        if(IsRegular_MetGrid)then
          iend = max(2,floor((abs(x_fullmet_sp(1)-xUR)/dx_met_const))+2)
        else
          ! Need to actually march through all the values of x_fullmet_sp and find
          ! which interval is needed
          iend = -1
          do i = 2,nx_fullmet
            if(x_fullmet_sp(i).ge.xUR.and.x_fullmet_sp(i-1).lt.xUR)iend=i
          enddo
          if(iend.lt.0)then
            write(MR_global_info,*)"MR ERROR: Could not find iend"
            write(MR_global_info,*)"  ",xUR
            stop 1
          endif
        endif
        nx_submet = iend-istart+1
        write(MR_global_info,*) "Domain is NOT periodic"
      endif

      !SEE if COMPUTATIONAL REGION STRADDLES THE BREAK IN THE WIND FILE
      !  (EITHER THE PRIME OR ANTI-MERIDIAN)
      if(iend.le.nx_fullmet)then        !yes
        wrapgrid = .false.
        write(MR_global_info,*)"Comp grid maps within a contiguous region of the Met grid"
        write(MR_global_info,*)"           wrapgrid = ",wrapgrid
        write(MR_global_info,*)"Met Sub grid specifications:"
        write(MR_global_info,*)"             istart = ",istart
        write(MR_global_info,*)"               iend = ",iend
        write(MR_global_info,*)"          nx_submet = ",nx_submet
        write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(istart)
        write(MR_global_info,*)"                xLL = ",xLL
        write(MR_global_info,*)"                xUR = ",xUR
        write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(iend)
      else                            !no
        if(IsGlobal_MetGrid)then
          wrapgrid = .true.

          ilhalf_fm_l = istart                        ! start index of left half on full met grid
          ilhalf_fm_r = nx_fullmet                    ! end index of left half on full met grid
          ilhalf_nx   = ilhalf_fm_r - ilhalf_fm_l +1  ! width of left half
          irhalf_fm_l = 1                             ! start index of right half on full met grid
          irhalf_fm_r = nx_submet - ilhalf_nx         ! end index of right half on full met grid
          irhalf_nx   = irhalf_fm_r - irhalf_fm_l +1  ! width of right half

          write(MR_global_info,*)"Comp grid span beyond the upper end of the Met grid"
          write(MR_global_info,*)"           wrapgrid = ",wrapgrid
          write(MR_global_info,*)"Met Sub grid specifications:"
          write(MR_global_info,*)"        ilhalf_fm_l = ",ilhalf_fm_l  ! start index of left half on full met grid
          write(MR_global_info,*)"        ilhalf_fm_r = ",ilhalf_fm_r  ! end index of left half on full met grid
          write(MR_global_info,*)"          ilhalf_nx = ",ilhalf_nx    ! width of left half
          write(MR_global_info,*)"        irhalf_fm_l = ",irhalf_fm_l  ! start index of right half on full met grid
          write(MR_global_info,*)"        irhalf_fm_r = ",irhalf_fm_r  ! end index of right half on full met grid
          write(MR_global_info,*)"          irhalf_nx = ",irhalf_nx    ! width of right half

          write(MR_global_info,*)"          nx_submet = ",nx_submet
          write(MR_global_info,*)"ilhalf_nx+irhalf_nx = ",ilhalf_nx+irhalf_nx
          write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(ilhalf_fm_l)
          write(MR_global_info,*)"                xLL = ",xLL
          write(MR_global_info,*)"                xUR = ",xUR
          write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(irhalf_fm_r)
        else
          write(MR_global_info,*)"MR ERROR: Comp grid requirements extend beyond Met grid"
          write(MR_global_info,*)"                xLL = ",xLL
          write(MR_global_info,*)"                xUR = ",xUR
          write(MR_global_info,*)"             istart = ",istart
          write(MR_global_info,*)"               iend = ",iend
          write(MR_global_info,*)"         xsubMetMin = ",x_fullmet_sp(1),x_fullmet_sp(istart)
          write(MR_global_info,*)"         xsubMetMax = ",x_fullmet_sp(iend),x_fullmet_sp(nx_fullmet)
          stop 1
        endif
      endif
      write(MR_global_info,*)"-------------"

      !SEE IF THE MODEL DOMAIN EXTENDS NORTH OR SOUTH OF THE MESOSCALE DOMAIN
      If(UseFullMetGrid)then
          ! This is the special case where the comp grid equals the Met grid
        jstart = 1
        jend = ny_fullmet
      else
        if(y_inverted)then
          ! Find start index
          if(y_fullmet_sp(1).ge.yUR)then
              ! This is the normal case where the UR of the comp grid is within
              ! the lat values of the wind file
            if(IsRegular_MetGrid)then
             jstart = max(1,floor((abs(y_fullmet_sp(1)-yUR)/dy_met_const))+1)
            else
               ! Need to actually march through all the values of y_fullmet_sp and find
               ! which interval is needed
              jstart = -1
              do j = 2,ny_fullmet
                if(y_fullmet_sp(j).lt.yUR.and.y_fullmet_sp(j-1).ge.yUR) jstart = j-1
              enddo
              if(jstart.lt.0)then
                write(MR_global_info,*)"MR ERROR: Could not find jstart"
                write(MR_global_info,*)"  ",yUR
                stop 1
              endif
            endif
          elseif(IsGlobal_MetGrid)then
              ! There are some special cases where the met grid is global, but do not
              ! have values at the poles (e.g. ERA).  There are occasional instances where we need
              ! values between the extreme lat value and the pole
            jstart = 1
            y_pad_North = .true.
          else
             write(MR_global_info,*)"MR ERROR: yUR > y_fullmet_sp(1)"
             write(MR_global_info,*)"     y_fullmet_sp(1).gt.yUR", &
                        y_fullmet_sp(1),yUR
             stop 1
          endif

          ! Find end index
          if(y_fullmet_sp(ny_fullmet).le.yLL)then
              ! Again, this is the normal case where the LL of the comp grid is within
              ! the lat values of the wind file
            if(IsRegular_MetGrid)then
              jend = max(2,floor((abs(y_fullmet_sp(1)-yLL)/dy_met_const))+2)
                ! Make sure jend doesn't exceed the met grid
              jend = min(jend,ny_fullmet)
            else
               ! Need to actually march through all the values of y_fullmet_sp and find
               ! which interval is needed
              jend = -1
              do j = 2,ny_fullmet
                if(y_fullmet_sp(j).lt.yLL.and.y_fullmet_sp(j-1).ge.yLL) jend = j
              enddo
              if(jend.lt.0)then
                write(MR_global_info,*)"MR ERROR: Could not find jend"
                write(MR_global_info,*)"  ",yLL
                stop 1
              endif
            endif
          elseif(IsGlobal_MetGrid)then
              ! Here is the same special case as above, but for the southern boundary
            jend = ny_fullmet
            y_pad_South = .true.
          else
             write(MR_global_info,*)"MR ERROR: y_fullmet_sp(ny_fullmet).lt.yLL",&
                               y_fullmet_sp(ny_fullmet),yLL
             stop 1
          endif

        else ! .not.y_inverted
          ! y values go from - to +
          if(y_fullmet_sp(1).le.yLL)then
            if(IsRegular_MetGrid)then
              jstart = max(1,floor((abs(y_fullmet_sp(1)-yLL)/dy_met_const)+1))
            else
              ! Need to actually march through all the values of y_fullmet_sp and find
              ! which interval is needed
              jstart = -1
              do j = 2,ny_fullmet
                if(y_fullmet_sp(j).ge.yLL.and.y_fullmet_sp(j-1).lt.yLL)jstart=j-1
              enddo
              if(jstart.lt.0)then
                write(MR_global_info,*)"MR ERROR: Could not find jend"
                write(MR_global_info,*)"  ",yLL
              endif
            endif
          elseif(IsGlobal_MetGrid)then
            jstart = 1
          else
            write(MR_global_info,*)"MR ERROR: yLL < y_fullmet_sp(1)"
            write(MR_global_info,*)"y_fullmet_sp(1),yLL",y_fullmet_sp(1),yLL,&
                      y_fullmet_sp(1).lt.yLL,y_fullmet_sp(1)-yLL
            stop 1
          endif
          if(y_fullmet_sp(ny_fullmet).ge.yUR)then
            if(IsRegular_MetGrid)then
              jend = max(2,floor((abs(y_fullmet_sp(1)-yUR)/dy_met_const))+2)
                ! Make sure jend doesn't exceed the met grid
              jend = min(jend,ny_fullmet)
            else
              ! Need to actually march through all the values of y_fullmet_sp and find
              ! which interval is needed
              jend = -1
              do j = 2,ny_fullmet
                if(y_fullmet_sp(j).ge.yUR.and.y_fullmet_sp(j-1).lt.yUR)jend=j
              enddo
              if(jend.lt.0)then
                write(MR_global_info,*)"MR ERROR: Could not find jend"
                write(MR_global_info,*)"  ",yUR
              endif
            endif
          elseif(IsGlobal_MetGrid)then
            jend = ny_fullmet
            y_pad_South = .true.
          else
            write(MR_global_info,*)"MR ERROR: yUR > y_fullmet_sp(ny_fullmet)"
            write(MR_global_info,*)"y_fullmet_sp(my_fullmet),yUR",y_fullmet_sp(ny_fullmet),yUr
            stop 1
          endif
        endif
      endif

      ! Calculate size of arrays that will hold the relavent section of
      ! the mesoscale model
      ny_submet = jend-jstart+1
      write(MR_global_info,*)"-------------"
      write(MR_global_info,*)"             jstart =" ,jstart
      write(MR_global_info,*)"               jend =" ,jend
      write(MR_global_info,*)"          ny_submet =" ,ny_submet
      write(MR_global_info,*)"         ysubMetMin =" ,y_fullmet_sp(jstart)
      write(MR_global_info,*)"                yLL =" ,yLL
      write(MR_global_info,*)"                yUR =" ,yUR
      write(MR_global_info,*)"         ysubMetMax = ",y_fullmet_sp(jend)
      write(MR_global_info,*)"-------------"

      if(IsPeriodic_CompGrid)then
        allocate( x_submet_sp(0:nx_submet+1))
        allocate( MR_dx_submet(0:nx_submet+1))
      else
        allocate( x_submet_sp(1:nx_submet))
        allocate( MR_dx_submet(1:nx_submet))
      endif
      allocate( y_submet_sp(1:ny_submet))
      allocate( MR_dy_submet(1:ny_submet))

      ! POPULATE THE X AND Y ARRAYS (Z ARRAY WILL BE POPULATED DIFFERENTLY
      ! AT EACH POINT, IN READ4DWINDARRAY.F90)
      if (wrapgrid) then
        x_submet_sp(          1:ilhalf_nx) = x_fullmet_sp(ilhalf_fm_l:ilhalf_fm_l+ilhalf_nx-1)
        x_submet_sp(ilhalf_nx+1:nx_submet) = x_fullmet_sp(irhalf_fm_l:irhalf_fm_l+irhalf_nx-1) + 360.0_sp
        MR_dx_submet(          1:ilhalf_nx) = MR_dx_met(ilhalf_fm_l:ilhalf_fm_l+ilhalf_nx-1)
        MR_dx_submet(ilhalf_nx+1:nx_submet) = MR_dx_met(irhalf_fm_l:irhalf_fm_l+irhalf_nx-1) + 360.0_sp
        if(IsPeriodic_CompGrid)then
          x_submet_sp(0)           = x_fullmet_sp(nx_submet  ) - 360.0_sp
          x_submet_sp(nx_submet+1) = x_fullmet_sp(nx_submet+1) + 360.0_sp
          MR_dx_submet(0)           = MR_dx_met(nx_submet  ) - 360.0_sp
          MR_dx_submet(nx_submet+1) = MR_dx_met(nx_submet+1) + 360.0_sp
        endif
      else
        x_submet_sp(1:nx_submet) = x_fullmet_sp(istart:iend)
        MR_dx_submet(1:nx_submet) = MR_dx_met(istart:iend)
        if(IsPeriodic_CompGrid)then
          x_submet_sp(0)           = x_fullmet_sp(nx_submet) - 360.0_sp
          x_submet_sp(nx_submet+1) = x_fullmet_sp(1     ) + 360.0_sp
          MR_dx_submet(0)           = MR_dx_met(nx_submet) - 360.0_sp
          MR_dx_submet(nx_submet+1) = MR_dx_met(1     ) + 360.0_sp
        endif
      endif

      do j=jstart,jend
        if(y_inverted)then
          y_submet_sp(jend-j+1) = y_fullmet_sp(j)
          MR_dy_submet(jend-j+1) = -MR_dy_met(j) ! Note that we need to negated dy so that
                                                 ! MR_dy_submet is always +
        else
          y_submet_sp(j-jstart+1) = y_fullmet_sp(j)
          MR_dy_submet(j-jstart+1) = MR_dy_met(j)
        endif
      enddo

      ! Set up for interpolation if needed
      write(MR_global_info,*)" Calculating mapping of comp "
      x_start_sub = x_submet_sp(1)
      y_start_sub = y_submet_sp(1)
      do i=1,nx_comp
        do j=1,ny_comp
          px = CompPoint_X_on_Met_sp(i,j)
          py = CompPoint_Y_on_Met_sp(i,j)
          if(IsLatLon_MetGrid)then
            ! When the met grid is lat/lon, we need to make sure that px
            ! map to the correct domain of the met file (-180->180 v.s. 0->360)
            if(IsGlobal_MetGrid)then
              if(IsPeriodic_CompGrid)then
                ! For global Lon/Lat Met data, allow points up to
                ! x_fullmet_sp(nx_fullmet)+dx
                if(px.ge.x_fullmet_sp(1)+360.0_sp)then
                  px=px-360.0_sp
                endif
              else
                ! For global Lon/Lat, but not a periodic comp grid
                if(px.gt.x_submet_sp(nx_submet))then
                  px=px-360.0_sp
                endif
                if(px.lt.x_submet_sp(1))then
                  px=px+360.0_sp
                endif
              endif
            else
              ! For non-global Met data, require values to be strictly within
              ! limits of the SUB-grid (i.e. might be >360)
              if(px.gt.x_submet_sp(nx_submet))then
                px=px-360.0_sp
              endif
              if(px.lt.x_submet_sp(1))then
                px=px+360.0_sp
              endif
            endif
          endif
          if(.not.IsPeriodic_CompGrid)then
            if(px.lt.x_start_sub.or.px.gt.x_submet_sp(nx_submet))then
              write(MR_global_info,*)"MR ERROR: Comp point maps out of sub_Met in x."
              write(MR_global_info,*)"Comp i,j, x      :",i,j,px
              write(MR_global_info,*)"sub_Met xmin,xmax:",x_start_sub,x_submet_sp(nx_submet)
              stop 1
            endif
            if((py.lt.y_start_sub           .and..not.y_pad_South).or.&
               (py.gt.y_submet_sp(ny_submet).and..not.y_pad_North))then
              write(MR_global_info,*)"MR ERROR: Comp point maps out of sub_Met in y."
              write(MR_global_info,*)"Comp i,j, y      :",i,j,py
              write(MR_global_info,*)"sub_Met ymin,ymax:",y_start_sub,y_submet_sp(ny_submet)
              stop 1
            endif
          endif

          ! Get the sub_Met index of LL corner of cell containing px,py
          If(IsRegular_MetGrid)then
            isubmet = int((px-x_start_sub)/dx_met_const) + 1
          else
            isubmet = -1
            do ii = 1,nx_submet
              if(px.ge.x_submet_sp(ii).and.px.lt.x_submet_sp(ii+1))then
                isubmet = ii
                exit
              endif
            enddo
          endif
          CompPoint_on_subMet_idx(i,j,1) = isubmet

            ! Check if the point is within the upper and lower bounds
          if(py.lt.y_submet_sp(ny_submet).and.py.ge.y_submet_sp(1))then
            If(IsRegular_MetGrid)then
              jsubmet = int((py-y_start_sub)/dy_met_const) + 1
            else
              jsubmet = -1
              do jj = 1,ny_submet-1
                if(py.ge.y_submet_sp(jj).and.py.lt.y_submet_sp(jj+1))then
                  jsubmet = jj
                  exit
                endif
              enddo
            endif
          elseif(abs(py-y_submet_sp(ny_submet)).lt.1.0e-7_sp)then
              ! This is to fix the occasional instances where the top comp point
              ! maps almost right on the top submet point
            If(IsRegular_MetGrid)then
              jsubmet = int((py-y_start_sub)/dy_met_const)
            else
              jsubmet = -1
              do jj = 1,ny_submet-1
                if(py.ge.y_submet_sp(jj).and.py.lt.y_submet_sp(jj+1))then
                  jsubmet = jj
                  exit
                endif
              enddo
            endif
          elseif(py.gt.y_submet_sp(ny_submet).and.IsGlobal_MetGrid.and.y_pad_North)then
            jsubmet = ny_submet
          elseif(py.lt.y_submet_sp(1)        .and.IsGlobal_MetGrid.and.y_pad_South)then
            jsubmet = 0
          endif
          CompPoint_on_subMet_idx(i,j,2) = jsubmet

          ! Get fractional position of comp point in met cell
          xfrac=(px-x_submet_sp(isubmet))/MR_dx_submet(isubmet)
          if(py.gt.y_submet_sp(ny_submet).and.IsGlobal_MetGrid.and.y_pad_North)then
              ! If comp point is above all met points
            yfrac=(py-y_submet_sp(ny_submet)+MR_dy_submet(ny_submet))/MR_dy_submet(ny_submet)
          elseif(py.lt.y_submet_sp(1).and.IsGlobal_MetGrid.and.y_pad_South)then
              ! If comp point is below all met points
            yfrac=(py-y_submet_sp(1)-MR_dy_submet(1))/MR_dy_submet(1)
          else
              ! Normal case where comp point is strictly within the met grid
            yfrac=(py-y_submet_sp(jsubmet))/MR_dy_submet(jsubmet)
          endif
          xc = 1.0_sp-xfrac
          yc = 1.0_sp-yfrac

          bilin_map_wgt(i,j,1)=xc*yc
          bilin_map_wgt(i,j,2)=xfrac*yc
          bilin_map_wgt(i,j,3)=xfrac*yfrac
          bilin_map_wgt(i,j,4)=yfrac*xc

        enddo
      enddo
      ! Now we have all the dimension sizes needed for allocating all grids

      ! We might need to rotate the wind vectors on the met grid in place if
      !   we need to convert ER to GR for the same grid or
      !   we need ER vectors from a projected Met grid
      if(.not.isGridRelative.or. &  ! We are dealing with NARR data
               Map_Case.eq.4.or. &  ! Met Grid is projected and Comp grid is Lat/Lon
               Map_Case.eq.5)then   ! Met Grid and Comp grids have different projections
        write(MR_global_info,*)"  Setting up for arrays for rotating vectors on Met grid."
        allocate(MR_u_ER_metP(nx_submet,ny_submet,np_fullmet))
        allocate(MR_v_ER_metP(nx_submet,ny_submet,np_fullmet))
        allocate(theta_Met(nx_submet,ny_submet))  ! This holds the angle between the projected
                                                  ! Met grid and the earth grid; used for rotating
                                                  ! grid velocities to Earth-Relative, or in the
                                                  ! special NARR case, rotating ER to GR

          !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + ip*(nc_nxmax*nc_nymax)
        do i=1,nx_submet
          do j=1,ny_submet
              ! Get lon/lat of point in question
            xin = real(x_submet_sp(i),kind=dp)
            yin = real(y_submet_sp(j),kind=dp)
            call PJ_proj_inv(xin,yin, Met_iprojflag, &
                          Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           ptlon,ptlat)
              ! Get projected coordinate of de at the currrent point
            call PJ_proj_for(ptlon+1.0_dp/60.0_dp,ptlat, Met_iprojflag, &
                       Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       xout,yout)
            de_x = xout-xin
            de_y = yout-yin
              ! Now recover the angle between de and x (of map grid)
            theta_Met(i,j) = atan(de_y/de_x)
          enddo
        enddo
      endif

      ! theta_Met ensures that we have Met data that is Earth-Relative, even if the
      ! underlying wind data are projected.  We might, however, need to rotate these
      ! ER values to a projected computational grid.  So we set up another rotation
      ! to map ER values that were interpolated onto a computational grid to GR
      if(Map_Case.eq.3.or. & ! Met is Lat/Lon, but Comp is projected
         Map_Case.eq.5)then  ! Met Grid and Comp grids have different projections
        write(MR_global_info,*)"  Setting up for arrays for rotating vectors on comp grid."
        allocate(MR_dum3d_compH_2(nx_comp,ny_comp,nz_comp))
        allocate(theta_Comp(nx_comp,ny_comp))
        ! If met and comp grids differ, first get Met grid winds as Earth-relative
        ! Note: This is only needed if Met grid is projected, ie for Map_Case = 4 or 5
        do i=1,nx_comp
          do j=1,ny_comp
              ! Get lon/lat of point in question
            xin = real(x_comp_sp(i),kind=dp)
            yin = real(y_comp_sp(j),kind=dp)
            call PJ_proj_inv(xin,yin, Comp_iprojflag, &
                          Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                           ptlon,ptlat)
              ! Get projected coordinate of de at the currrent point
            call PJ_proj_for(ptlon+1.0_dp/60.0_dp,ptlat, Comp_iprojflag, &
                       Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                       xout,yout)
            de_x = xout-xin
            de_y = yout-yin
              ! Again, we want the angle between de and x (but now, of comp grid)
            theta_Comp(i,j) = atan(de_y/de_x)
          enddo
        enddo
      endif

      allocate(MR_dum2d_met_int(nx_submet,ny_submet))
      allocate(MR_dum2d_met(nx_submet,ny_submet))
      allocate(MR_dum3d_metP(nx_submet,ny_submet,np_fullmet))
      allocate(MR_dum3d2_metP(nx_submet,ny_submet,np_fullmet))
      allocate(MR_dum3d_metH(nx_submet,ny_submet,nz_comp))
      allocate(MR_dum2d_comp_int(nx_comp,ny_comp))
      allocate(MR_dum2d_comp(nx_comp,ny_comp))
      allocate(MR_dum3d_compH(nx_comp,ny_comp,nz_comp))
      !  The only 3d met data that persists locally is Geopotential Height.
      !  This is needed for interpolating the other variables onto a catersian
      !  grid.
      allocate(MR_geoH_metP_last(nx_submet,ny_submet,np_fullmet))
      allocate(MR_geoH_metP_next(nx_submet,ny_submet,np_fullmet))

      if(MR_iwindformat.eq.25)then
        ! The following is only needed if we are reading 2d variables from NCEP
        ! wind files, but we will set up the grids regardless as a precaution
          ! Here are the weights starting with LL then counter-clockwise
        allocate(amap_iwf25(nx_submet,ny_submet,4))
        amap_iwf25 = 0.0_sp
        !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + sp*(as_nxmax*as_nymax*4)
          ! Here is the x1,x2,y1,y2 indicies
        allocate(imap_iwf25(nx_submet,ny_submet,4))
        imap_iwf25 = 0
        !Tot_Bytes_on_Heap = Tot_Bytes_on_Heap + sp*(as_nxmax*as_nymax*4)

        do ilon = 1,nx_submet
          x_loc = max(0.0_sp,x_submet_sp(ilon))
          do i = 1,191
            if(max(0.0_sp,x_in_iwf25_sp(i)).le.x_loc.and.x_in_iwf25_sp(i+1).gt.x_loc)then
              ix1 = i
              ix2 = i+1
              xfrac_sp = (x_loc - x_in_iwf25_sp(i))/1.875_sp
              xc_sp = 1.0_sp - xfrac_sp
              exit ! leave do loop
            endif
          enddo

          do ilat = 1,ny_submet
            y_loc = y_submet_sp(ilat)
            do j = 94,2,-1
              if(y_in_iwf25_sp(j).le.y_loc.and.y_in_iwf25_sp(j-1).gt.y_loc)then
                iy1 = j
                iy2 = j-1
                dely_sp = y_in_iwf25_sp(j-1)-y_in_iwf25_sp(j)
                yfrac_sp = (y_loc - y_in_iwf25_sp(j))/dely_sp
                yc_sp = 1.0_sp - yfrac_sp
                exit ! leave do loop
              endif
            enddo

            imap_iwf25(ilon,ilat,1)=ix1
            imap_iwf25(ilon,ilat,2)=ix2
            imap_iwf25(ilon,ilat,3)=iy1
            imap_iwf25(ilon,ilat,4)=iy2
            amap_iwf25(ilon,ilat,1)=xc_sp*yc_sp
            amap_iwf25(ilon,ilat,2)=xfrac_sp*yc_sp
            amap_iwf25(ilon,ilat,3)=xfrac_sp*yfrac_sp
            amap_iwf25(ilon,ilat,4)=yfrac_sp*xc_sp

          enddo
        enddo

      endif

      write(MR_global_info,*)"Finished MR_Set_MetComp_Grids"
      write(MR_global_info,*)" "
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_MetComp_Grids

!##############################################################################


!##############################################################################
!
!     MR_Set_Comp2Met_Map
!
!     This subroutine generates the mapping between the computational grid and
!     the Met grid.  
!
!     Allocates the mapping arrays:
!
!        CompPoint_on_subMet_idx
!        bilin_map_wgt
!        CompPoint_X_on_Met_sp
!        CompPoint_Y_on_Met_sp
!
!##############################################################################


      subroutine MR_Set_Comp2Met_Map

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=dp), parameter :: tol = 1.0e-3_dp

      real(kind=dp) :: dum1,dum2,dum3,dum4,dum5
      integer :: i,j
      real(kind=dp) :: x_in ,y_in
      real(kind=dp) :: x_out,y_out

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                          MR_Set_Comp2Met_Map               ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      ! We now have the full definition of the Met grid and the Comp grid
      ! Figure out if we need to do any remapping
      !Possibilities are:
      !(1) Both Comp Grid and Met grids are Lat/Lon
      !(2) Both Comp Grid and Met grids are the same projection
      !(3) Met Grid is Lat/Lon and Comp grid is projected
      !(4) Met Grid is projected and Comp grid is Lat/Lon
      !(5) Met Grid and Comp grids have different projections
      !      For 3,4,5, we need to map the comp points to the Met grid
      !      Set up CompPoint_Metx, CompPoint_Mety (these are the xy (lon/lat) of each comp point)
      !             CompPoint_Meti, CompPoint_Metj (these are the i,j indices of the sub-Met grid)
      !             CompPoint_Met_Wgt (weights given to the four surrounding points
      if(IsLatLon_MetGrid.and.IsLatLon_CompGrid)then
        ! Both Comp and Met are in Lat/Lon 
        Map_Case = 1
      elseif(IsLatLon_MetGrid.and..not.IsLatLon_CompGrid)then
        ! Met is Lat/Lon, but Comp is projected
        Map_Case = 3
      elseif(.not.IsLatLon_MetGrid.and.IsLatLon_CompGrid)then
        ! Met is projected, but Comp is Lat/Lon
        Map_Case = 4
      else
        ! Both Met and Comp are projected.
        ! Test if the projections are the same.
        if(Met_iprojflag.ne.Comp_iprojflag)then
          ! Met and Comp are completely different projection types
          Map_Case = 5
        else
          ! Projections are the same type, test individual parameters
          if(Comp_iprojflag.eq.1)then
            ! Polar stereographic
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_k0 - Met_k0)
            dum4 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.2)then
            ! Albers Equal Area
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_phi1 - Met_phi1)
            dum4 = abs(Comp_phi2 - Met_phi2)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.3)then
            ! UTM
            stop 1
          elseif(Comp_iprojflag.eq.4)then
            ! Lambert conformal conic (NARR, NAM218, NAM221)
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_phi1 - Met_phi1)
            dum4 = abs(Comp_phi2 - Met_phi2)
            dum5 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol.and.&
               dum4.gt.tol.and.&
               dum5.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          elseif(Comp_iprojflag.eq.5)then
            ! Mercator (NAM196)
            dum1 = abs(Comp_lam0 - Met_lam0)
            dum2 = abs(Comp_phi0 - Met_phi0)
            dum3 = abs(Comp_Re - Met_Re)
            if(dum1.gt.tol.and.&
               dum2.gt.tol.and.&
               dum3.gt.tol)then
              Map_Case = 5
            else
              Map_Case = 2
            endif
          endif !Comp_iprojflag
        endif !Met_iprojflag.ne.Comp_iprojflag
      endif ! Met and Comp projected

      if(Map_Case.eq.1)then
        write(MR_global_info,*)"Map_case = ",Map_case
        write(MR_global_info,*)"  Both Comp Grid and Met grids are Lat/Lon"
      elseif(Map_Case.eq.2)then
        write(MR_global_info,*)"Map_case = ",Map_case
        write(MR_global_info,*)"  Both Comp Grid and Met grids are the same projection"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      elseif(Map_Case.eq.3)then
        write(MR_global_info,*)"Map_case = ",Map_case
        write(MR_global_info,*)"  Met Grid is Lat/Lon and Comp grid is projected"
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      elseif(Map_Case.eq.4)then
        write(MR_global_info,*)"Map_case = ",Map_case
        write(MR_global_info,*)"  Met Grid is projected and Comp grid is Lat/Lon"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
      elseif(Map_Case.eq.5)then
        write(MR_global_info,*)"Map_case = ",Map_case
        write(MR_global_info,*)"  Met Grid and Comp grids have different projections"
        write(MR_global_info,*)"   Met_iprojflag",Met_iprojflag
        write(MR_global_info,*)"   Met_lam0",real(Met_lam0,kind=sp)
        write(MR_global_info,*)"   Met_phi0",real(Met_phi0,kind=sp)
        write(MR_global_info,*)"   Met_phi1",real(Met_phi1,kind=sp)
        write(MR_global_info,*)"   Met_phi2",real(Met_phi2,kind=sp)
        write(MR_global_info,*)"   Met_Re  ",real(Met_Re,kind=sp)
        write(MR_global_info,*)"   Met_k0  ",real(Met_k0,kind=sp)
        write(MR_global_info,*)"   Comp_iprojflag",Comp_iprojflag
        write(MR_global_info,*)"   Comp_lam0",real(Comp_lam0,kind=sp)
        write(MR_global_info,*)"   Comp_phi0",real(Comp_phi0,kind=sp)
        write(MR_global_info,*)"   Comp_phi1",real(Comp_phi1,kind=sp)
        write(MR_global_info,*)"   Comp_phi2",real(Comp_phi2,kind=sp)
        write(MR_global_info,*)"   Comp_Re  ",real(Comp_Re,kind=sp)
        write(MR_global_info,*)"   Comp_k0  ",real(Comp_k0,kind=sp)
      endif

      allocate(CompPoint_on_subMet_idx(nx_comp,ny_comp,2))
      allocate(bilin_map_wgt(nx_comp,ny_comp,4))
      allocate(CompPoint_X_on_Met_sp(nx_comp,ny_comp))
      allocate(CompPoint_Y_on_Met_sp(nx_comp,ny_comp))

      if(Map_Case.eq.1.or.Map_Case.eq.2)then
        ! Map and Comp are on same grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            CompPoint_X_on_Met_sp(i,j) = real(x_in,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_in,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.3)then
          ! We just need to map the projected comp grid to the Lon/Lat Met grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            call PJ_proj_inv(x_in, y_in, Comp_iprojflag, &
                           Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.4)then
          ! We just need to map the Lon/Lat comp grid to the projected Met grid
        do i=1,nx_comp
          x_in = x_comp_sp(i)
          do j=1,ny_comp
            y_in = y_comp_sp(j)
            call PJ_proj_for(x_in, y_in, Met_iprojflag, &
                           Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      elseif(Map_Case.eq.5)then
          ! Here, we need to map the projected comp grid to a Lon/Lat grid, then
          ! map to the projected Met grid
        do i=1,nx_comp
          do j=1,ny_comp
            x_in = x_comp_sp(i)
            y_in = y_comp_sp(j)
            call PJ_proj_inv(x_in, y_in, Comp_iprojflag, &
                           Comp_lam0,Comp_phi0,Comp_phi1,Comp_phi2,Comp_k0,Comp_Re, &
                           x_out,y_out)
            x_in = x_out
            y_in = y_out
            call PJ_proj_for(x_in, y_in, Met_iprojflag, &
                           Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                           x_out,y_out)
            CompPoint_X_on_Met_sp(i,j) = real(x_out,kind=sp)
            CompPoint_Y_on_Met_sp(i,j) = real(y_out,kind=sp)
          enddo
        enddo
      endif ! Map_Case

      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_Comp2Met_Map

!##############################################################################
!
!     MR_Regrid_Met2Comp
!
!     This subroutine does the 2-D regridding using the mapping
!     arrays filled in MR_Set_Comp2Met_Map
!
!##############################################################################


      subroutine MR_Regrid_Met2Comp(nx1,ny1,wrk_met,nx2,ny2,wrk_comp)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer                         ,intent(in)  :: nx1,ny1
      real(kind=sp),dimension(nx1,ny1),intent(in)  :: wrk_met
      integer                         ,intent(in)  :: nx2,ny2
      real(kind=sp),dimension(nx2,ny2),intent(out) :: wrk_comp

      integer :: i,j,ii,jj
      integer :: nx_max
      real(kind=sp) :: a1,a2,a3,a4
      real(kind=sp) :: tmp

      real(kind=sp),dimension(:,:),allocatable :: wrk_loc

      if(IsPeriodic_CompGrid)then
        nx_max = nx1+1
        allocate(wrk_loc(0:nx1+1,0:ny1+1))
        if(y_pad_North)then
          tmp = sum(wrk_met(1:nx1,ny1))/real(nx1,kind=sp)
          wrk_loc(0:nx1+1,ny1+1) = tmp
        else
          wrk_loc(0:nx1+1,ny1+1) = 0.0_sp
        endif
        if(y_pad_South)then
          tmp = sum(wrk_met(1:nx1,1))/real(nx1,kind=sp)
          wrk_loc(0:nx1+1,0) = tmp
        else
          wrk_loc(0:nx1+1,0) = 0.0_sp
        endif
        wrk_loc(1:nx1,1:ny1) = wrk_met(1:nx1,1:ny1) 
        wrk_loc(0    ,1:ny1) = wrk_met(  nx1,1:ny1)
        wrk_loc(nx1+1,1:ny1) = wrk_met(1    ,1:ny1)
      else
        nx_max = nx1
        allocate(wrk_loc(nx1,0:ny1+1))
        if(y_pad_North)then
          tmp = sum(wrk_met(1:nx1,ny1))/real(nx1,kind=sp)
          wrk_loc(1:nx1,ny1+1) = tmp
        else
          wrk_loc(1:nx1,ny1+1) = 0.0_sp
        endif
        wrk_loc(1:nx1,1:ny1) = wrk_met(1:nx1,1:ny1)
      endif

      ! Loop over all comp points
      do i = 1,nx2
        do j = 1,ny2
          ! Get the Met cell id this comp point maps to
          ii = CompPoint_on_subMet_idx(i,j,1)
          jj = CompPoint_on_subMet_idx(i,j,2)
          if(ii.lt.1.or.ii.gt.nx_max-1)then
            write(MR_global_info,*)"MR ERROR: ii maps out of grid: ",ii
            stop 1
          endif
          if(jj.lt.0.or.jj.gt.ny1)then
            write(MR_global_info,*)"MR ERROR: jj maps out of grid: ",jj,ny1
            stop 1
          endif
          ! Look up this comp points weights
          a1 = bilin_map_wgt(i,j,1)
          a2 = bilin_map_wgt(i,j,2)
          a3 = bilin_map_wgt(i,j,3)
          a4 = bilin_map_wgt(i,j,4)
          ! Now interpolate from Met to Comp
          wrk_comp(i,j) = a1*wrk_loc(ii  ,jj  ) + &
                          a2*wrk_loc(ii+1,jj  ) + &
                          a3*wrk_loc(ii+1,jj+1) + &
                          a4*wrk_loc(ii  ,jj+1)

        enddo
      enddo


      end subroutine MR_Regrid_Met2Comp

!##############################################################################
!
!     MR_Regrid_P2H_linear
!
!     This subroutine does the 1-D regridding (replaces rgrd1d) linearly
!     interpolating in z.
!
!##############################################################################


      subroutine MR_Regrid_P2H_linear(nzm,z_met ,var_met, &
                                      nzc,z_comp,var_comp)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer                     ,intent(in)  :: nzm
      real(kind=sp),dimension(nzm),intent(in)  :: z_met
      real(kind=sp),dimension(nzm),intent(in)  :: var_met
      integer                     ,intent(in)  :: nzc
      real(kind=sp),dimension(nzc),intent(in)  :: z_comp
      real(kind=sp),dimension(nzc),intent(out) :: var_comp

      integer :: km,kc,km_interv
      real(kind=sp) :: a1,a2,dz,z1

      logical :: found_interv

      var_comp = -9999.0_sp
      ! Loop over all comp points
      km_interv = 1
      do kc = 1,nzc
        found_interv = .false.
        z1 = z_comp(kc)
        ! For each comp point, check which met interval it is in, starting from
        ! the last interval found
        do km = km_interv,nzm-1
        !do km = 1,nzm-1
          if(z1.ge.z_met(km).and.z1.le.z_met(km+1))then
            found_interv = .true.
            km_interv = km
            dz = z_met(km_interv+1)-z_met(km_interv)
            a1 = (z1-z_met(km_interv))/dz
            a2 = 1.0_sp-a1
            var_comp(kc) = var_met(km_interv  ) * a2 + &
                           var_met(km_interv+1) * a1
            exit
          else
            cycle
          endif
        enddo
        ! Check that interval was found
        if(.not.found_interv)then
          write(MR_global_info,*)"MR ERROR:  Did not find interval in vertical 1-D interpolation."
          write(MR_global_info,*)"z_met = "
          write(MR_global_info,*)z_met
          write(MR_global_info,*)" "
          write(MR_global_info,*)"z_comp = "
          write(MR_global_info,*)z_comp
          stop 1
        endif
      enddo

      return

      end subroutine MR_Regrid_P2H_linear

!##############################################################################
!
!     MR_Read_Met_Template
!
!     This subroutine reads an auxillary file that specifies the custom windfile
!     structure
!
!##############################################################################


      subroutine MR_Read_Met_Template

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      logical      :: IsThere
      character(len=130) :: lllinebuffer
      character(len=80)  :: Met_projection_line

      integer :: ioerr
      character(len=3)   :: useLeap_str
      integer :: ndims_custom, nvars_custom
      integer :: i
      character     :: dv_char
      integer       :: dimID,varID,zindx,vndim
      real(kind=sp) :: fac
      character(len=30) :: dname
      character(len=71) :: vname
      character(len=5)  :: vname_WMO
      real(kind=8)      :: StepInterval

      write(MR_global_info,*)"Inside MR_Read_Met_Template"
      inquire( file=adjustl(trim(MR_iwf_template)), exist=IsThere )
      if(.not.IsThere)then
        write(MR_global_info,*)"MR ERROR: Could not find NWP template file ",&
                   adjustl(trim(MR_iwf_template))
        write(MR_global_info,*)"          Make sure the calling program sets MR_iwf_template"
        write(MR_global_info,*)"          and that it is linked to the cwd."
        stop 1
      endif

      open(unit=27,file=adjustl(trim(MR_iwf_template)),status='unknown')
      read(27,'(a130)')lllinebuffer
      Met_projection_line = lllinebuffer(1:80)
      call PJ_Set_Proj_Params(Met_projection_line)

      if(PJ_ilatlonflag.eq.0)then
        IsLatLon_MetGrid = .false.
        IsGlobal_MetGrid = .false.
      else
        IsLatLon_MetGrid = .true.
        if(PJ_iprojflag.eq.0)then
          IsGlobal_MetGrid = .false.
        else
          IsGlobal_MetGrid = .true.
        endif
      endif
      Met_iprojflag = PJ_iprojflag
      Met_lam0      = PJ_lam0
      Met_lam1      = PJ_lam1
      Met_lam2      = PJ_lam2
      Met_phi0      = PJ_phi0
      Met_phi1      = PJ_phi1
      Met_phi2      = PJ_phi2
      Met_k0        = PJ_k0
      Met_Re        = PJ_radius_earth

      read(27,'(a130)')lllinebuffer
      read(lllinebuffer,*,iostat=ioerr)StepInterval,useLeap_str

      if (ioerr.eq.0)then
        ! Two values read, process useLeap_str to determine T or F
        if(useLeap_str(1:1).eq.'F'.or.useLeap_str(1:1).eq.'f')then
          MR_useLeap = .false.
          write(MR_global_info,*)"This windfile template specifies that leap years are NOT to"
          write(MR_global_info,*)"be used.  Resetting MR_useLeap = .false."
        endif
      else
        ! default will be whatever is set in the host program or in
        ! MetReader.f90 if the calling program doesn't specify
        read(lllinebuffer,*,iostat=ioerr)StepInterval
      endif
      read(27,'(a130)')lllinebuffer
      read(lllinebuffer,*,err=2002)ndims_custom,nvars_custom

      write(MR_global_info,*)"  Reading dimensions: ",ndims_custom
      do i = 1,ndims_custom
        read(27,'(a130)')lllinebuffer
        read(lllinebuffer,1501)dv_char,dimID,fac,dname
        if(dv_char.ne.'d')then
          write(MR_global_info,*)"MR ERROR : Trying to read variable into dimension"
          write(MR_global_info,*)"dv_char = ",dv_char
          write(MR_global_info,*)"dimID   = ",dimID
          write(MR_global_info,*)"fac     = ",fac
          write(MR_global_info,*)"dname   = ",dname
          stop 1
        endif
        if(dimID.le.9)then
          Met_dim_IsAvailable(dimID) = .true.
          Met_dim_names(dimID)       = adjustl(trim(dname))
          Met_dim_fac(i)             = fac
        else
          write(MR_global_info,*)"MR ERROR: dimID too large",dimID
          stop 1
        endif
      enddo
      write(MR_global_info,*)"  Reading variables: ",nvars_custom
      do i = 1,nvars_custom
        read(27,'(a130)')lllinebuffer
        read(lllinebuffer,1511)dv_char,vndim,zindx,varID, &
                                fac,vname_WMO,vname
        if(dv_char.ne.'v')then
          write(MR_global_info,*)"MR ERROR : Trying to read dimension into variable"
          write(MR_global_info,*)"dv_char   = ",dv_char
          write(MR_global_info,*)"vndim     = ",vndim
          write(MR_global_info,*)"zindx     = ",zindx
          write(MR_global_info,*)"varID     = ",varID
          write(MR_global_info,*)"fac       = ",fac
          write(MR_global_info,*)"vname_WMO = ",vname_WMO
          write(MR_global_info,*)"vname     = ",vname
          stop 1
        endif

        if(varID.le.50)then
          Met_var_IsAvailable(varID)       = .true.
          Met_var_names(varID)             = adjustl(trim(vname))
          Met_var_names_WMO(varID)         = adjustl(trim(vname_WMO))
          Met_var_ndim(varID)              = vndim
          Met_var_zdimID(varID)            = zindx
          Met_var_conversion_factor(varID) = fac
          write(MR_global_info,*)varID,Met_var_names_WMO(varID),' ',Met_var_names(varID)
        else
          write(MR_global_info,*)"MR ERROR: varID too large",varID
          stop 1
        endif
      enddo

      close(27)

      return

1501  format(a1,i9     ,f9.2,a30)
1511  format(a1,i3,i3,i3,f9.2,a7,a71)

!2001  write(MR_global_info,*)  'error reading ForecastInterval'
!      write(MR_global_info,*)lllinebuffer
!      stop 1
2002  write(MR_global_info,*)  'error reading number of custom dims and vars.'
      write(MR_global_info,*)lllinebuffer
      stop 1

      end subroutine MR_Read_Met_Template

