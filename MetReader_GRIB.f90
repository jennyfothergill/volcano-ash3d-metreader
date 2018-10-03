!##############################################################################
!##############################################################################
!##############################################################################


!##############################################################################
!
!     MR_Read_Met_DimVars_GRIB
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_GRIB_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_GRIB

      use MetReader
      use grib_api
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: i, j, k
      real(kind=sp) :: xLL_fullmet
      real(kind=sp) :: yLL_fullmet
      real(kind=sp) :: xUR_fullmet
      real(kind=sp) :: yUR_fullmet

      integer            :: ifile
      integer            :: iret
      integer            :: igrib
      integer            :: iw

      integer          :: ivar,iivar
      integer          :: idx
      integer          :: maxdimlen

      character(len=130) :: grib_file_path
      integer(kind=4)  :: iv_discpl
      integer(kind=4)  :: iv_paramC
      integer(kind=4)  :: iv_paramN
      integer(kind=4)  :: iv_typeSf
      character(len=7) :: iv_marsParam
      character(len=3) :: iv_typeSfc

      integer(kind=4)  :: numberOfPoints
      integer(kind=4)  :: Ni,Nx
      integer(kind=4)  :: Nj,Ny
      integer(kind=4)  :: dum_int
      real(kind=8)     :: dum_dp
      character(len=19) :: dum_str
      real(kind=dp) :: x_start,y_start
      real(kind=dp) :: Lon_start,Lat_start
      real(kind=dp),dimension(:),allocatable     :: lats,lons,values

      integer(kind=4)  :: typeOfFirstFixedSurface
      integer            :: count1=0
        ! Stores values of keys read from grib file
      character(len=7) :: grb_marsParam
      character(len=4) :: grb_typeSfc
      integer(kind=4)  :: grb_discipline
      integer(kind=4)  :: grb_parameterCategory
      integer(kind=4)  :: grb_parameterNumber
      integer(kind=4)  :: grb_level
      integer(kind=4)  :: grb_scaledValueOfFirstFixedSurface

      integer(kind=4),dimension(MR_MAXVARS,100) :: zlev_dum  ! This will hold the z-levels, up to 100
      integer(kind=4),dimension(MR_MAXVARS)     :: zcount    ! This will hold the length of the z-coord
      logical :: Check
      logical :: FoundOldDim
      logical :: IsTruncatedDim
      logical :: ReadGrid

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Read_Met_DimVars_GRIB                    ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

        !---------------------------------------------------------------------------------
        ! Checking for dimension length and values for x,y,t,p
        !   Assume all files have the same format

      if(MR_iwindformat.eq.27)then
        stop 5
      else
        write(*,*)"Opening grib file to find version number"
        iw = 1
        call grib_open_file(ifile,trim(ADJUSTL(MR_windfiles(iw))),'R')
        call grib_new_from_file(ifile,igrib,iret)

        call grib_get(igrib,'editionNumber',MR_GRIB_Version)
        call grib_release(igrib)
        call grib_close_file(ifile)
      endif
      write(*,*)"Grib version = ",MR_GRIB_Version
      !---------------------------------------------------------------------------------
      ! Checking for dimension length and values for x,y,t,p
      !   Assume all files have the same format
      maxdimlen = 0
      if(MR_GRIB_Version.eq.1)then

      elseif(MR_GRIB_Version.eq.2)then
        ! Loop through all the grib messages,
        ! If we find a message that matches a variable criteria, then log the level to 
        !  a dummy array.
        ! Finally sort the pressure values and evaluate the distinct pressure coordinates
        grib_file_path  = adjustl(trim(MR_windfiles(1)))

        call grib_open_file(ifile,grib_file_path,'R')
        ! Loop on all the messages in a file.
        call grib_new_from_file(ifile,igrib,iret)
        count1=0
        zcount(:)     = 0
        zlev_dum(:,:) = 0
        do while (iret/=GRIB_END_OF_FILE)
          count1=count1+1
          if(count1.eq.1)then
            ! For the first record, get the x,y grid info
            ReadGrid = .false.
            call grib_get(igrib,'Nx',nx_fullmet)
            call grib_get(igrib,'Ny',ny_fullmet)
            allocate(x_fullmet_sp(0:nx_fullmet+1))
            allocate(y_fullmet_sp(ny_fullmet))
            allocate(MR_dx_met(nx_fullmet))
            allocate(MR_dy_met(ny_fullmet))

            call grib_get(igrib,'gridType',dum_str)
            call grib_get(igrib,'latitudeOfFirstGridPointInDegrees',dum_dp)
            Lat_start = dum_dp
            call grib_get(igrib,'longitudeOfFirstGridPointInDegrees',dum_dp)
            Lon_start = dum_dp

            dum_int = 0
            Met_Re =  6371.229_8
            call grib_get(igrib,'shapeOfTheEarth',dum_int)
            if (dum_int.eq.0)then
                ! 0  Earth assumed spherical with radius = 6,367,470.0 m
              Met_Re =  6367.470_8
            elseif(dum_int.eq.1)then
                ! 1  Earth assumed spherical with radius specified by data producer
                !  Try to read the radius of earth
                !  For now, just assign the default
              Met_Re =  6371.229_8
            elseif(dum_int.eq.2)then
                ! 2  Earth assumed oblate spheroid with size as determined by IAU in 1965
                !    (major axis = 6,378,160.0 m, minor axis = 6,356,775.0 m, f = 1/297.0)
              Met_Re =  6371.229_8
            elseif(dum_int.eq.3)then
                ! 3  Earth assumed oblate spheroid with major and minor axes specified by data producer
              Met_Re =  6371.229_8
            elseif(dum_int.eq.4)then
                ! 4  Earth assumed oblate spheroid as defined in IAG-GRS80 model 
                !    (major axis = 6,378,137.0 m, minor axis = 6,356,752.314 m, f = 1/298.257222101)
              Met_Re =  6371.229_8
            elseif(dum_int.eq.5)then
                ! 5  Earth assumed represented by WGS84 (as used by ICAO since 1998)
              Met_Re =  6371.229_8
            elseif(dum_int.eq.6)then
                ! 6  Earth assumed spherical with radius of 6,371,229.0 m
              Met_Re =  6371.229_8
            else
                ! 7-191 Reserved
                ! 192- 254  Reserved for local use
              Met_Re =  6371.229_8
            endif

            if(index(dum_str,'regular_ll').ne.0)then
              IsLatLon_MetGrid = .true.
              Lat_start = y_start
              Lon_start = x_start

              call grib_get(igrib,'numberOfPoints',numberOfPoints)
              allocate(lats(numberOfPoints))
              allocate(lons(numberOfPoints))
              allocate(values(numberOfPoints))
              call grib_get_data(igrib,lats,lons,values)
              do j=1,ny_fullmet
                do i=1,nx_fullmet
                  idx = (j-1)*nx_fullmet + i
                  x_fullmet_sp(i) = lons(idx)
                  y_fullmet_sp(j) = lats(idx)
                enddo
              enddo
              ReadGrid = .true.
              deallocate(lats)
              deallocate(lons)
              deallocate(values)
              call grib_get(igrib,'latitudeOfLastGridPointInDegrees',dum_dp)
              if(Lat_start.gt.dum_dp)then
                y_inverted = .true.
              else
                y_inverted = .false.
              endif
              call grib_get(igrib,'longitudeOfLastGridPointInDegrees',dum_dp)
              if(Lon_start.gt.dum_dp)then
                x_inverted = .true.
              else
                x_inverted = .false.
              endif
              call grib_get(igrib,'iDirectionIncrementInDegrees',dum_dp)
              dx_met_const = real(dum_dp,kind=4)
              call grib_get(igrib,'jDirectionIncrementInDegrees',dum_dp)
              dy_met_const = real(dum_dp,kind=4)
              x_fullmet_sp(0) = x_fullmet_sp(1)-dx_met_const
              x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(nx_fullmet)+dx_met_const
            elseif(index(dum_str,'regular_gg').ne.0)then
              IsLatLon_MetGrid = .true.
              Lat_start = y_start
              Lon_start = x_start

              call grib_get(igrib,'numberOfPoints',numberOfPoints)
              allocate(lats(numberOfPoints))
              allocate(lons(numberOfPoints))
              allocate(values(numberOfPoints))
              call grib_get_data(igrib,lats,lons,values)
              do j=1,ny_fullmet
                do i=1,nx_fullmet
                  idx = (j-1)*nx_fullmet + i
                  x_fullmet_sp(i) = lons(idx)
                  y_fullmet_sp(j) = lats(idx)
                enddo
              enddo
              ReadGrid = .true.
              deallocate(lats)
              deallocate(lons)
              deallocate(values)
              call grib_get(igrib,'latitudeOfLastGridPointInDegrees',dum_dp)
              if(Lat_start.gt.dum_dp)then
                y_inverted = .true.
              else
                y_inverted = .false.
              endif
              call grib_get(igrib,'longitudeOfLastGridPointInDegrees',dum_dp)
              if(Lon_start.gt.dum_dp)then
                x_inverted = .true.
              else
                x_inverted = .false.
              endif
              call grib_get(igrib,'iDirectionIncrementInDegrees',dum_dp)
              dx_met_const = real(dum_dp,kind=4)
              call grib_get(igrib,'jDirectionIncrementInDegrees',dum_dp)
              dy_met_const = real(dum_dp,kind=4)
              x_fullmet_sp(0) = x_fullmet_sp(1)-dx_met_const
              x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(nx_fullmet)+dx_met_const

            elseif(index(dum_str,'polar_stereographic').ne.0)then
              IsLatLon_MetGrid = .false.
              Met_iprojflag     = 1
              write(*,*)"polar_stereographic not implemented"
              stop 1

            elseif(index(dum_str,'albers').ne.0)then
              IsLatLon_MetGrid = .false.
              Met_iprojflag     = 2
              write(*,*)"Alber Equal Area not implemented"
              stop 1

            elseif(index(dum_str,'UTM').ne.0)then
              IsLatLon_MetGrid = .false.
              Met_iprojflag     = 3
              write(*,*)"UTM not implemented"
              stop 1
            elseif(index(dum_str,'lambert').ne.0)then
              IsLatLon_MetGrid = .false.
              Met_iprojflag     = 4
              call grib_get(igrib,'LoVInDegrees',dum_dp)
              Met_lam0 = dum_dp
              call grib_get(igrib,'LaDInDegrees',dum_dp)
              Met_phi0 = dum_dp
              call grib_get(igrib,'Latin1InDegrees',dum_dp)
              Met_phi1 = dum_dp
              call grib_get(igrib,'Latin2InDegrees',dum_dp)
              Met_phi2 = dum_dp
              Met_k0   =  0.933_8
              Met_Re   =  6371.229_8
            elseif(index(dum_str,'mercator').ne.0)then
              IsLatLon_MetGrid = .false.
              Met_iprojflag     = 5
              Met_lam0 = Lon_start
              call grib_get(igrib,'LaDInDegrees',dum_dp)
              Met_phi0 = dum_dp
              Met_k0   =  0.933_8
              Met_Re   =  6371.229_8

            else
              write(MR_global_error,*)'MR ERROR: Cannot determine the projection from the GRIB file.'
              stop 1
            endif

            if(.not.IsLatLon_MetGrid)then
              call grib_get(igrib,'DxInMetres',dum_int)
              dx_met_const = real(dum_int,kind=4)/1000.0_sp
              call grib_get(igrib,'DyInMetres',dum_int)
              dy_met_const = real(dum_int,kind=4)/1000.0_sp

              call PJ_proj_for(Lon_start,Lat_start, Met_iprojflag, &
                       Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
            endif

            if(.not.ReadGrid)then
              do i = 0,nx_fullmet+1
                x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
              enddo
              if(y_inverted)then
                do i = 1,ny_fullmet
                  y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
                enddo
              else
                do i = 1,ny_fullmet
                  y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=sp)
                enddo
              endif
              ReadGrid = .true.
            endif
            do i = 1,nx_fullmet
              MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
            enddo
            do i = 1,ny_fullmet-1
              MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
            enddo
            MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)
              
          endif

          call grib_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface)
          ! for populating z-levels, we are only concerned with specific level types
          if(typeOfFirstFixedSurface.eq.100.or. & ! Isobaric surface  (Pa)
             typeOfFirstFixedSurface.eq.103.or. & ! Specified height level above ground  (m)
             typeOfFirstFixedSurface.eq.106)then  ! Depth below land surface  (m)
            call grib_get(igrib,'discipline',              grb_discipline)
            call grib_get(igrib,'parameterCategory',       grb_parameterCategory)
            call grib_get(igrib,'parameterNumber',         grb_parameterNumber)

            ! Loop through all the variables and see if we have a match with this grib record
            do ivar = 1,MR_MAXVARS
              if (.not.Met_var_IsAvailable(ivar)) cycle
              iv_discpl = Met_var_GRIB2_DPcPnSt(ivar,1)
              iv_paramC = Met_var_GRIB2_DPcPnSt(ivar,2)
              iv_paramN = Met_var_GRIB2_DPcPnSt(ivar,3)
              iv_typeSf = Met_var_GRIB2_DPcPnSt(ivar,4)
              if(iv_discpl.eq.grb_discipline.and.       &
                 iv_paramC.eq.grb_parameterCategory.and.&
                 iv_paramN.eq.grb_parameterNumber.and.  &
                 iv_typeSf.eq.typeOfFirstFixedSurface)then
                ! This is one we are tracking, log the level
                call grib_get(igrib,'scaledValueOfFirstFixedSurface',grb_scaledValueOfFirstFixedSurface)
                zcount(ivar) = zcount(ivar) + 1
                zlev_dum(ivar,zcount(ivar)) = grb_scaledValueOfFirstFixedSurface
              endif
            enddo
          endif
          call grib_release(igrib)
          call grib_new_from_file(ifile,igrib,iret)
        enddo
        call grib_release(igrib)
        call grib_close_file(ifile)

        maxdimlen = maxval(zcount(:))
        nlev_coords_detected = 1
        Met_var_zdim_idx(1) = 1
        do ivar = 2,MR_MAXVARS
          FoundOldDim = .false.
          do iivar = 1,ivar-1
            if (zcount(iivar).eq.zcount(ivar))then  ! This check for a different coordinate is
                                                    ! soley on the size of the dimension.
              FoundOldDim = .true.
              Met_var_zdim_idx(ivar)  = Met_var_zdim_idx(iivar)
              exit
            endif
          enddo
          if(.not.FoundOldDim)then
            nlev_coords_detected = nlev_coords_detected + 1
            Met_var_zdim_idx(ivar)  = nlev_coords_detected
          endif
        enddo
        ! The V part of velocity is typically the second of a multi-component record and the above code
        ! does not catch it.  Assume V has the same characteristics as U and copy U
        ! V @ isobaric
        zcount(3) = zcount(2)
        zlev_dum(3,1:zcount(3)) = zlev_dum(2,1:zcount(2))
        Met_var_zdim_idx(3) = Met_var_zdim_idx(2)
        ! V @ 10 m
        zcount(12) = zcount(11)
        zlev_dum(12,1:zcount(12)) = zlev_dum(11,1:zcount(11))
        Met_var_zdim_idx(12) = Met_var_zdim_idx(11)

        ! We have all the level dimension names and dim_ids; now we need to get the sizes
        allocate(nlevs_fullmet(nlev_coords_detected))
        allocate(levs_code(nlev_coords_detected))
        allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
        do ivar = 1,MR_MAXVARS
          if (.not.Met_var_IsAvailable(ivar)) cycle
          idx = Met_var_zdim_idx(ivar)
          nlevs_fullmet(idx) = zcount(ivar)
          ! Check that the pressure level is monotonically increasing
          Check=.true.
          do k=1,zcount(ivar)-1
            if(zlev_dum(ivar,k).gt.zlev_dum(ivar,k+1)) Check=.false.
          enddo
          if (Check)then
            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = &
              real(zlev_dum(ivar,1:nlevs_fullmet(idx)),kind=sp)
            z_inverted = .true.
            ! Pressure is expected to be from the bottom up so invert
            allocate(p_fullmet_sp(nlevs_fullmet(idx)))
            do k=1,zcount(ivar)
              p_fullmet_sp(k) = levs_fullmet_sp(idx,nlevs_fullmet(idx)+1-k)
            enddo
            levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = p_fullmet_sp(1:nlevs_fullmet(idx))
            deallocate(p_fullmet_sp)
          else
            ! Need to first sort pressure variable here
            stop 6
          endif
        enddo

        ! Now log all pressure coordinates as one-to-one, truncated, or interupted
        levs_code(1:nlev_coords_detected) = 0
        levs_code(1) = 1                       ! The first var checked (GPH) should have a one-to-one mapping
        ! Check how each of the pressure coordinates map onto the GPH grid
        if (nlev_coords_detected.gt.1)then
          ! Only bother if there are multiple pressure coordinates
          do idx = 2,nlev_coords_detected
            if (nlevs_fullmet(idx).gt.nlevs_fullmet(1))then
              ! This coordinate has more values than the GPH pressure coordinate
              levs_code(idx) = 4
            elseif (nlevs_fullmet(idx).lt.nlevs_fullmet(1))then
              ! It there are fewer levels, check if this is a truncated coordiante (code = 2)
              ! or one with missing levels that requires interpolation (code = 3)
              IsTruncatedDim = .true.
              do i=1,nlevs_fullmet(idx)
                if(abs(levs_fullmet_sp(idx,i)-levs_fullmet_sp(1,i)).gt.MR_EPS_SMALL)then
                  IsTruncatedDim = .false.
                  exit
                endif
              enddo
              if(IsTruncatedDim)then
                levs_code(idx) = 2
              else
                levs_code(idx) = 3
              endif
            else
              ! This coordinate has the same dimension as the GPH pressure coordinate.
              ! They are probably the same
              levs_code(idx) = 1
            endif
          enddo
        endif

        write(MR_global_production,*)" Found these levels"
        write(MR_global_production,*)"  VaribleID    LevelIdx       dimID"
        do ivar = 1,MR_MAXVARS
          if (Met_var_IsAvailable(ivar)) &
              write(MR_global_production,*)ivar,Met_var_zdim_idx(ivar)
        enddo

        ! Now assign these levels to the working arrays
        ! Geopotential is the first variable checked, use this for np_fullmet
        nt_fullmet = 1
        np_fullmet = nlevs_fullmet(Met_var_zdim_idx(1))  ! Assign fullmet the length of H,U,V
        allocate(p_fullmet_sp(np_fullmet))
        idx = Met_var_zdim_idx(1)
        p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp)
      else
        write(MR_global_error,*)"MR ERROR : MR_GRIB_Version must be 1 or 2, not  ",&
                                MR_GRIB_Version
        stop 1
      endif

      allocate(z_approx(np_fullmet))
      do k=1,np_fullmet
        ! Calculate heights for US Std Atmos while pressures are still in mbars
        ! or hPa
        z_approx(k) = MR_Z_US_StdAtm(p_fullmet_sp(k))
      enddo

      write(MR_global_info,*)"Dimension info:"
      write(MR_global_info,*)"  record (time): ",nt_fullmet
      write(MR_global_info,*)"  level  (z)   : ",np_fullmet
      write(MR_global_info,*)"  y            : ",ny_fullmet
      write(MR_global_info,*)"  x            : ",nx_fullmet

      !************************************************************************
      ! assign boundaries of mesoscale model
      if(x_inverted)then
          ! I know of no windfiles with x-coordinate reversed
        xLL_fullmet = x_fullmet_sp(nx_fullmet)
        xUR_fullmet = x_fullmet_sp(1)
      else
        xLL_fullmet = x_fullmet_sp(1)
        xUR_fullmet = x_fullmet_sp(nx_fullmet)
      endif

      if(y_inverted)then
          ! Most lon/lat grids have y reversed
        yLL_fullmet = y_fullmet_sp(ny_fullmet)
        yUR_fullmet = y_fullmet_sp(1)
      else
          ! Projected grids have y not reversed
        yLL_fullmet = y_fullmet_sp(1)
        yUR_fullmet = y_fullmet_sp(ny_fullmet)
      endif

      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_Met_DimVars_GRIB

!##############################################################################

!##############################################################################
!
!     MR_Read_Met_Times_GRIB
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine opens each GRIB file and determine the time of each
!     time step of each file in the number of hours since MR_BaseYear.
!     In most cases, the length of the time variable (nt_fullmet) will be 
!     read directly from the file and overwritten (is was set in MR_Read_Met_DimVars_GRIB
!     above).
!
!     After this subroutine completes, the following variables will be set:
!       MR_windfile_starthour(MR_iwindfiles)
!       MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
!
!##############################################################################

      subroutine MR_Read_Met_Times_GRIB

      use MetReader
      use grib_api

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: iw,iws
      integer :: itstart_year,itstart_month
      integer :: itstart_day
      real(kind=sp) :: filestart_hour

      integer :: itstart_hour,itstart_min,itstart_sec

      real(kind=8)       :: HS_hours_since_baseyear
      character(len=130) :: dumstr
      integer            :: iwstep

      integer            :: dataDate
      integer            :: dataTime
      integer            :: forecastTime
      integer            :: ifile
      integer            :: iret
      integer            :: igrib

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Read_Met_Times_GRIB                      ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(.not.Met_dim_IsAvailable(1))then
        write(MR_global_error,*)"MR ERROR: Time dimension is required and not listed"
        write(MR_global_error,*)"          in custom windfile specification file."
        stop 1
      endif

      allocate(MR_windfile_starthour(MR_iwindfiles))
      if(MR_iwindformat.eq.27)then
        ! GRIB1 reader not yet working!!
        write(MR_global_error,*)"MR ERROR: iwf=27 is a GRIB1 format."
        write(MR_global_error,*)"       The GRIB1 reader is not yet working"
        stop 1
        ! Here the branch for when MR_iwindformat = 27
        ! First copy path read in to slot 2
        !if(MR_runAsForecast)then
        !  write(MR_global_error,*)"MR ERROR: iwf=27 cannot be used for forecast runs."
        !  write(MR_global_error,*)"          These are reanalysis files."
        !  stop 1
        !endif
        dumstr = MR_windfiles(1)
 110    format(a50,a1,i4,a1)
        write(MR_windfiles(1),110)trim(ADJUSTL(dumstr)),'/', &
                                   MR_Comp_StartYear,'/'
        write(MR_windfiles(2),110)trim(ADJUSTL(dumstr)),'/', &
                                   MR_Comp_StartYear+1,'/'
        MR_windfile_starthour(1) = real(HS_hours_since_baseyear( &
                                    MR_Comp_StartYear,1,1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
        MR_windfile_starthour(2) = real(HS_hours_since_baseyear( &
                                    MR_Comp_StartYear+1,1,1,0.0_8,MR_BaseYear,MR_useLeap),kind=sp)
        if  ((mod(MR_Comp_StartYear,4).eq.0)     .and.                     &
             (mod(MR_Comp_StartYear,100).ne.0).or.(mod(MR_Comp_StartYear,400).eq.0))then
          nt_fullmet = 1464     ! Leap year
        else
          nt_fullmet = 1460     ! Not a leap year
        endif
        MR_windfiles_nt_fullmet(1)=nt_fullmet
        MR_windfiles_nt_fullmet(2)=nt_fullmet  ! Note: we don't care if the next
                                               !       year is a leap year since
                                               !       the simulation will never
                                               !       be long enough to matter.

        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

          ! the interval for iwf27 is 6 hours
        do iwstep = 1,nt_fullmet
          MR_windfile_stephour(:,iwstep) = (iwstep-1)*6.0_4
        enddo
      else
        ! For all other formats, try to read the first grib message and get
        ! dataDate, dataTime and forecastTime
        ! Loop through all the windfiles
        do iw = 1,MR_iwindfiles

          ! Each wind file needs a ref-time which in almost all cases is given
          ! in the 'units' attribute of the time variable
          write(MR_global_info,*)iw,trim(ADJUSTL(MR_windfiles(iw)))

          if(iw.eq.1)then
            ! For now, assume one time step per file
            nt_fullmet = 1
            write(MR_global_info,*)"  Assuming all NWP files have the same number of steps."
            write(MR_global_info,*)"   For grib, assume one time step per file."
            write(MR_global_info,*)"   Allocating time arrays for ",MR_iwindfiles,"files"
            write(MR_global_info,*)"                              ",nt_fullmet,"step(s) each"
            allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
          endif

          call grib_open_file(ifile,trim(ADJUSTL(MR_windfiles(iw))),'R')
          call grib_new_from_file(ifile,igrib,iret)

          if(iw.eq.1)call grib_get(igrib,'editionNumber',MR_GRIB_Version)

          call grib_get(igrib,'dataDate',dataDate)
          call grib_get(igrib,'dataTime',dataTime)
          if(MR_GRIB_Version.eq.1)then
            ! The only grib1 files we deal with are reanalysis files with no FC time
            forecastTime = 0
          else
            call grib_get(igrib,'forecastTime',forecastTime)
          endif

          itstart_year  = int(dataDate/10000)
          itstart_month = int((dataDate-10000*itstart_year)/100)
          itstart_day   = mod(dataDate,100)
          itstart_hour  = int(dataTime/100)
          itstart_min   = mod(dataTime,100)
          itstart_sec   = 0

          write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                     itstart_hour,itstart_min,itstart_sec

          call grib_release(igrib)
          call grib_close_file(ifile)

          filestart_hour = real(itstart_hour,kind=sp) + &
                           real(itstart_min,kind=sp)/60.0_sp      + &
                           real(itstart_sec,kind=sp)/3600.0_sp

          MR_windfiles_nt_fullmet(iw)=nt_fullmet
          MR_windfile_starthour(iw) =  real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                         itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
          MR_windfile_stephour(iw,1) = real(forecastTime,kind=4)

        enddo
      endif
2100  format(20x,a11,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)

      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)

      write(MR_global_info,*)"File, step, Ref, Offset, HoursSince"
      do iw = 1,MR_iwindfiles
        do iws = 1,nt_fullmet
          write(MR_global_info,*)iw,iws,real(MR_windfile_starthour(iw),kind=4),&
                           real(MR_windfile_stephour(iw,iws),kind=4),&
                           real(MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws),kind=4)
        enddo
      enddo

      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_Met_Times_GRIB
!##############################################################################

!##############################################################################
!
!     MR_Read_MetP_Variable_GRIB
!
!     Called from Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_GRIB(ivar,istep)

      use MetReader
      use grib_api

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: ivar
      integer,intent(in) :: istep

      integer :: iw,iwstep
      integer :: np_met_loc
      character(len=71)  :: invar
      character(len=130) :: index_file
      character(len=130) :: grib_file_path
      character(len=130) :: grib_file

      real(kind=sp) :: del_H,del_P,dpdz

      integer :: i,j,k
      integer :: kk
      integer :: kkk,itmp
      integer :: ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      integer :: Dimension_of_Variable
      logical :: IsCatagorical

      integer,dimension(np_fullmet) :: p_met_loc

      integer            :: ifile
      integer            :: igrib
      integer            :: idx
      integer            :: iret
      integer            :: l,m,t
      integer            :: count1=0
      integer            :: rstrt, rend
      real(kind=8),dimension(:),allocatable     :: values
      real(kind=8),dimension(:,:),allocatable   :: slice
      integer(kind=4)  :: numberOfPoints
      integer(kind=4)  :: Ni
      integer(kind=4)  :: Nj
      integer(kind=4)  :: typeOfFirstFixedSurface
        ! Stores values of keys read from grib file
      character(len=7) :: grb_marsParam
      character(len=4) :: grb_typeSfc
      integer(kind=4)  :: grb_discipline
      integer(kind=4)  :: grb_parameterCategory
      integer(kind=4)  :: grb_parameterNumber
      integer(kind=4)  :: grb_level
      integer(kind=4)  :: grb_scaledValueOfFirstFixedSurface
      !character(len=9)   :: sName
      !integer(kind=4)  :: forecastTime
      integer(kind=4),dimension(:),allocatable :: marsParam_idx
      integer(kind=4),dimension(:),allocatable :: discipline_idx
      integer(kind=4),dimension(:),allocatable :: parameterCategory_idx
      integer(kind=4),dimension(:),allocatable :: parameterNumber_idx
      integer(kind=4),dimension(:),allocatable :: level_idx
      integer(kind=4),dimension(:),allocatable :: forecastTime_idx

      integer(kind=4)  :: marsParamSize
      integer(kind=4)  :: disciplineSize
      integer(kind=4)  :: parameterCategorySize
      integer(kind=4)  :: parameterNumberSize
      integer(kind=4)  :: levelSize
      integer(kind=4)  :: forecastTimeSize
        ! temporary versions of those stored in Met_var_GRIB[...](ivar)
      integer(kind=4)  :: iv_discpl
      integer(kind=4)  :: iv_paramC
      integer(kind=4)  :: iv_paramN
      integer(kind=4)  :: iv_typeSf
      character(len=7) :: iv_marsParam
      character(len=3) :: iv_typeSfc

      real(kind=sp) :: Z_top, T_top
      real(kind=sp),dimension(:,:,:),allocatable :: full_values

      logical :: Use_GRIB_Index = .false.
      integer :: fn_idx

      if(.not.Met_var_IsAvailable(ivar))then
        write(MR_global_error,*)"MR ERROR:  Variable not available for this windfile"
        write(MR_global_error,*)"             ivar = ",ivar
        write(MR_global_error,*)"            vname = ",Met_var_GRIB_names(ivar)
        write(MR_global_error,*)"             iwf  = ",MR_iwindformat
        stop 1
      endif

      if(ivar.eq.3 .or. &   ! V_isobaric
         ivar.eq.12)then    ! V_height_above_ground
        ! The v components might be not reachable with the indexing if they are the second
        ! component of a multi-variable message
        Use_GRIB_Index = .false.
!      else
!        Use_GRIB_Index = .true.
      endif

      if(MR_GRIB_Version.eq.1)then
        ! grib1 uses an index based on marsParam
        iv_marsParam = Met_var_GRIB1_MARS(ivar)
        iv_typeSfc   = Met_var_GRIB1_St(ivar)
      elseif(MR_GRIB_Version.eq.2)then
        ! Get the variable discipline, Parameter Catagory, Parameter Number, and
        ! level type for this variable
        iv_discpl = Met_var_GRIB2_DPcPnSt(ivar,1)
        iv_paramC = Met_var_GRIB2_DPcPnSt(ivar,2)
        iv_paramN = Met_var_GRIB2_DPcPnSt(ivar,3)
        iv_typeSf = Met_var_GRIB2_DPcPnSt(ivar,4)
      else
        write(MR_global_error,*)"MR ERROR:  GRIB type not determined"
        stop 1
      endif

      iw     = MR_MetStep_findex(istep)
      iwstep = MR_MetStep_tindex(istep)

      if(Met_var_GRIB_names(ivar).eq."")then
        write(MR_global_error,*)"Variable ",ivar," not available for MR_iwindformat = ",&
                  MR_iwindformat
        stop 1
      endif

      ! Get the dimension of the variable requested (either 2 or 3-D)
      if(ivar.eq.1 ) Dimension_of_Variable = 3 ! Geopotential Height
      if(ivar.eq.2 ) Dimension_of_Variable = 3 ! Vx
      if(ivar.eq.3 ) Dimension_of_Variable = 3 ! Vy
      if(ivar.eq.4 ) Dimension_of_Variable = 3 ! Vz
      if(ivar.eq.5 ) Dimension_of_Variable = 3 ! Temperature
      if(ivar.eq.6 ) Dimension_of_Variable = 3 ! Pressure (only for WRF or other eta-level files)

      if(ivar.eq.10) Dimension_of_Variable = 2 ! Planetary Boundary Layer Height
      if(ivar.eq.11) Dimension_of_Variable = 2 ! U @ 10m
      if(ivar.eq.12) Dimension_of_Variable = 2 ! V @ 10m
      if(ivar.eq.13) Dimension_of_Variable = 2 ! Friction velocity
      if(ivar.eq.14) Dimension_of_Variable = 2 ! Displacement Height
      if(ivar.eq.15) Dimension_of_Variable = 2 ! Snow cover
      if(ivar.eq.16) Dimension_of_Variable = 2 ! Soil moisture
      if(ivar.eq.17) Dimension_of_Variable = 2 ! Surface roughness
      if(ivar.eq.18) Dimension_of_Variable = 2 ! Wind_speed_gust_surface

      if(ivar.eq.20) Dimension_of_Variable = 2 ! pressure at lower cloud base
      if(ivar.eq.21) Dimension_of_Variable = 2 ! pressure at lower cloud top
      if(ivar.eq.22) Dimension_of_Variable = 2 ! temperature at lower cloud top
      if(ivar.eq.23) Dimension_of_Variable = 2 ! Total Cloud cover
      if(ivar.eq.24) Dimension_of_Variable = 2 ! Cloud cover (low)
      if(ivar.eq.25) Dimension_of_Variable = 2 ! Cloud cover (convective)

      if(ivar.eq.30) Dimension_of_Variable = 3 ! Rel. Hum
      if(ivar.eq.31) Dimension_of_Variable = 3 ! QV (specific humidity)
      if(ivar.eq.32) Dimension_of_Variable = 3 ! QL (liquid)
      if(ivar.eq.33) Dimension_of_Variable = 3 ! QI (ice)

      if(ivar.eq.40) Dimension_of_Variable = 2 ! Categorical rain
      if(ivar.eq.41) Dimension_of_Variable = 2 ! Categorical snow
      if(ivar.eq.42) Dimension_of_Variable = 2 ! Categorical frozen rain
      if(ivar.eq.43) Dimension_of_Variable = 2 ! Categorical ice
      if(ivar.eq.44) Dimension_of_Variable = 2 ! Precipitation rate large-scale (liquid)
      if(ivar.eq.45) Dimension_of_Variable = 2 ! Precipitation rate convective (liquid)
      if(ivar.eq.46) Dimension_of_Variable = 3 ! Precipitation rate large-scale (ice)
      if(ivar.eq.47) Dimension_of_Variable = 3 ! Precipitation rate convective (ice)

      if(ivar.eq.40.or.&
         ivar.eq.41.or.&
         ivar.eq.42.or.&
         ivar.eq.43)then
          ! Catagorical variables are integers and need special interpolation
        IsCatagorical = .true.
      else
          ! The default is to read floating point values
        IsCatagorical = .false.
      endif

      if(MR_iwindformat.eq.27)then
        ! Get correct GRIB1 file
        write(MR_global_error,*)"MR ERROR: iwf27 not working for GRIB1"
        stop 1
        if(ivar.eq.1)then
          write(index_file,125)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_HGT_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.2)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_UGRD_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.3)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_VGRD_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.4)then
          write(index_file,126)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_VVEL_pres.nc"
          write(*,*)"+++++++++++++++++++++++++++++++++++++++++++"
          write(*,*)"  NEED TO FIX THIS Vz"
          np_met_loc = np_fullmet
        elseif(ivar.eq.5)then
          write(index_file,125)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_TMP_pres.nc"
          np_met_loc = np_fullmet
        elseif(ivar.eq.10)then
          write(index_file,128)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_HPBL_sfc.nc"
        elseif(ivar.eq.22)then
          write(index_file,130)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_TMP_low-cldtop.nc"
        elseif(ivar.eq.23)then
          write(index_file,131)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_TCDC_low-cldlay.nc"
        elseif(ivar.eq.30)then
          write(index_file,127)trim(adjustl(MR_MetStep_File(istep))), &
                           "pgrbanl_mean_",MR_iwind5_year(istep), &
                           "_RH_pres.nc"
          write(*,*)"+++++++++++++++++++++++++++++++++++++++++++"
          write(*,*)"  NEED TO FIX THIS : RH"
          np_met_loc = np_fullmet
        elseif(ivar.eq.44)then
          write(index_file,129)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_PRATE_sfc.nc"
        elseif(ivar.eq.45)then
          write(index_file,129)trim(adjustl(MR_MetStep_File(istep))), &
                           "sflxgrbfg_mean_",MR_iwind5_year(istep), &
                           "_CPRAT_sfc.nc"
        else
          write(MR_global_error,*)"MR ERROR : Requested variable not available."
          stop 1
        endif
        index_file = trim(adjustl(index_file))

 125      format(a50,a13,i4,a12)
 126      format(a50,a13,i4,a13)
 127      format(a50,a13,i4,a11)
 128      format(a50,a15,i4,a12)
 129      format(a50,a15,i4,a13)
 130      format(a50,a15,i4,a18)
 131      format(a50,a15,i4,a19)
      else  ! all other cases besides iwf27
        ! Set up pressure level index that we will search for
        p_met_loc = 0
        idx = Met_var_zdim_idx(ivar)
        if(ivar.eq.4)then      ! Vertical_velocity_pressure_isobaric
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        elseif(ivar.eq.10)then ! Planetary_Boundary_Layer_Height_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.11)then ! u-component_of_wind_height_above_ground
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 10
        elseif(ivar.eq.12)then ! v-component_of_wind_height_above_ground
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 10
        elseif(ivar.eq.13)then ! Frictional_Velocity_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.15)then ! Snow_depth_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.16)then ! Volumetric_Soil_Moisture_Content_depth_below_surface_layer
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.17)then ! Surface_roughness_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.18)then ! Wind_speed_gust_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.20)then ! Pressure_cloud_base
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.21)then ! Pressure_cloud_topw
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.23)then ! Total_cloud_cover_entire_atmosphere
           ! Something is wrong reading this
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.30)then ! Relative_humidity_isobaric
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        elseif(ivar.eq.40.or.ivar.eq.41.or.ivar.eq.42.or.ivar.eq.43)then ! categorical precip
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.44)then ! Precipitation_rate_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        else
          np_met_loc = nlevs_fullmet(idx)
          if(MR_GRIB_Version.eq.1)then
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx))/100.0_sp)
          else
            p_met_loc(1:np_met_loc)  = int(levs_fullmet_sp(idx,1:nlevs_fullmet(idx)))
          endif
        endif
        allocate(full_values(nx_fullmet,ny_fullmet,np_met_loc))
          ! Files are listed directly, not through directories (as in MR_iwindformat=25,27)
        grib_file_path  = trim(adjustl(MR_MetStep_File(istep)))
        fn_idx = index(MR_MetStep_File(istep), '/' , BACK=.true.)
        grib_file = MR_MetStep_File(istep)(fn_idx+1:)
        index_file      = trim(adjustl(MR_MetStep_File(istep))) // ".index"
      endif
      invar = Met_var_GRIB_names(ivar)

      ! Load data variables for just the subgrid defined above
      if (wrapgrid) then
        ict        = 2
          ! index on the sub-met
        ileft(1)   = 1;         ileft(2)   = ilhalf_nx+1
        iright(1)  = ilhalf_nx; iright(2)  = nx_submet
          ! indes on the full-met
        iistart(1) = ilhalf_fm_l; iistart(2) = irhalf_fm_l
        iicount(1) = ilhalf_nx  ; iicount(2) = irhalf_nx
      else
        ict        = 1
        ileft(1)   = 1
        iright(1)  = nx_submet
        iistart(1) = istart
        iicount(1) = nx_submet
      endif

      !----------------------------------------------------
      ! This is the part where we actually read the files
      ! 4 branches total
      !  Grib1 with index file
      !  Grib1 without index file
      !  Grib2 with index file
      !  Grib2 without index file
      if(MR_GRIB_Version.eq.1)then
        if(Use_GRIB_Index)then
          write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                    trim(adjustl(index_file))

          call grib_index_read(idx,index_file)
          call grib_multi_support_on()

            ! get the number of distinct values of all the keys in the index
          call grib_index_get_size(idx,'marsParam',marsParamSize)
          call grib_index_get_size(idx,'level',levelSize)

            ! allocate the array to contain the list of distinct values
          allocate(marsParam_idx(marsParamSize))
          allocate(level_idx(levelSize))

            ! get the list of distinct key values from the index
          call grib_index_get(idx,'marsParam',marsParam_idx)
          call grib_index_get(idx,'level',level_idx)

          ! Start marching throught the index file and look for the match with the 
          ! keys
          count1=0
          do l=1,marsParamSize
            call grib_index_select(idx,'marsParam',marsParam_idx(l))
            do i=1,levelSize
              call grib_index_select(idx,'level',level_idx(i))
              call grib_new_from_index(idx,igrib, iret)
              do while (iret /= GRIB_END_OF_INDEX)
                count1=count1+1
                call grib_get(igrib,'indicatorOfTypeOfLevel',grb_typeSfc)

                if( grb_marsParam           .eq. iv_marsParam.and. &
                    grb_typeSfc(1:3)        .eq. iv_typeSfc) then
                  call grib_get(igrib,'numberOfPoints',numberOfPoints)
                  call grib_get(igrib,'Ni',Ni)
                  call grib_get(igrib,'Nj',Nj)
                  if(nx_fullmet.ne.Ni)then
                    write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                    write(MR_global_error,*)"nx_fullmet = ",nx_fullmet
                    write(MR_global_error,*)"Ni         = ",Ni
                    stop 1
                  endif
                  if(ny_fullmet.ne.Nj)then
                    write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                    write(MR_global_error,*)"ny_fullmet = ",ny_fullmet
                    write(MR_global_error,*)"Nj         = ",Nj
                    stop 1
                  endif
                  allocate(values(numberOfPoints))
                  allocate(slice(Ni,Nj))
                  call grib_get(igrib,'values',values)
                  do m = 1,Nj
                    rstrt = (m-1)*Ni + 1
                    rend  = m*Ni
                    slice(1:Ni,m) = values(rstrt:rend)
                  enddo
                  deallocate(values)

                  ! There is no guarentee that grib levels are in order so...
                  ! Now loop through the pressure values for this variable and put
                  ! this slice at the correct level.
                  call grib_get(igrib,'level',grb_level)
                  do kk = 1,np_met_loc
                    if(p_met_loc(kk).eq.grb_level)then
                      full_values(:,:,kk) = real(slice(:,:),kind=sp)
                      exit
                    endif
                  enddo
                  deallocate(slice)
                endif

                call grib_release(igrib)
                call grib_new_from_index(idx,igrib, iret)
              enddo ! while
              call grib_release(igrib)

            enddo ! loop on level
          enddo ! loop on marsParam

          call grib_index_release(idx)

        else ! Non-index Grib1 case

          ! We don't have/(can't make) the index file so scan all messages of the
          ! grib1 file
          write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                    trim(adjustl(grib_file_path))!,nx_submet,ny_submet,np_met_loc
          ifile=5
          call grib_open_file(ifile,grib_file_path,'R')

          !     turn on support for multi fields messages */
          call grib_multi_support_on()

          ! Loop on all the messages in a file.
          call grib_new_from_file(ifile,igrib,iret)
          count1=0

          do while (iret/=GRIB_END_OF_FILE)
            count1=count1+1

            call grib_get(igrib,'marsParam',grb_marsParam)
            call grib_get(igrib,'indicatorOfTypeOfLevel',grb_typeSfc)

            if ( grb_marsParam           .eq. iv_marsParam.and. &
                 grb_typeSfc(1:3)        .eq. iv_typeSfc) then
              call grib_get(igrib,'numberOfPoints',numberOfPoints)
              call grib_get(igrib,'Ni',Ni)
              call grib_get(igrib,'Nj',Nj)
              if(nx_fullmet.ne.Ni)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"nx_fullmet = ",nx_fullmet
                write(MR_global_error,*)"Ni         = ",Ni
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"ny_fullmet = ",ny_fullmet
                write(MR_global_error,*)"Nj         = ",Nj
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
              call grib_get(igrib,'values',values)
              do m = 1,Nj
                rstrt = (m-1)*Ni + 1
                rend  = m*Ni
                slice(1:Ni,m) = values(rstrt:rend)
              enddo
              deallocate(values)

               ! There is no guarentee that grib levels are in order so...
               ! Now loop through the pressure values for this variable and put
               ! this slice at the correct level.
               call grib_get(igrib,'level',grb_level)
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.grb_level)then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
            endif
            call grib_release(igrib)
            call grib_new_from_file(ifile,igrib, iret)

          enddo
          call grib_close_file(ifile)

        endif

      elseif(MR_GRIB_Version.eq.2)then
        if(Use_GRIB_Index)then
          write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                    trim(adjustl(index_file))
          call grib_index_read(idx,index_file)
          call grib_multi_support_on()
    
            ! get the number of distinct values of all the keys in the index
          call grib_index_get_size(idx,'discipline',disciplineSize)
          call grib_index_get_size(idx,'parameterCategory',parameterCategorySize)
          call grib_index_get_size(idx,'parameterNumber',parameterNumberSize)
          call grib_index_get_size(idx,'scaledValueOfFirstFixedSurface',levelSize)
          call grib_index_get_size(idx,'forecastTime',forecastTimeSize)
          
            ! allocate the array to contain the list of distinct values
          allocate(discipline_idx(disciplineSize))
          allocate(parameterCategory_idx(parameterCategorySize))
          allocate(parameterNumber_idx(parameterNumberSize))
          allocate(level_idx(levelSize))
          allocate(forecastTime_idx(forecastTimeSize))
          
            ! get the list of distinct key values from the index
          call grib_index_get(idx,'discipline',discipline_idx)
          call grib_index_get(idx,'parameterCategory',parameterCategory_idx)
          call grib_index_get(idx,'parameterNumber',parameterNumber_idx)
          call grib_index_get(idx,'scaledValueOfFirstFixedSurface',level_idx)
          call grib_index_get(idx,'forecastTime',forecastTime_idx)
    
          ! Start marching throught the index file and look for the match with the 
          ! keys
          count1=0
          do l=1,disciplineSize
            call grib_index_select(idx,'discipline',discipline_idx(l))
        
            do j=1,parameterCategorySize
              call grib_index_select(idx,'parameterCategory',parameterCategory_idx(j))
        
              do k=1,parameterNumberSize
                call grib_index_select(idx,'parameterNumber',parameterNumber_idx(k))
        
                do i=1,levelSize
                  call grib_index_select(idx,'level',level_idx(i))
                  call grib_index_select(idx,'scaledValueOfFirstFixedSurface',level_idx(i))
   
                  do t=1,forecastTimeSize
                    call grib_index_select(idx,'forecastTime',forecastTime_idx(t))
                    call grib_new_from_index(idx,igrib, iret)
                    do while (iret /= GRIB_END_OF_INDEX)
                      count1=count1+1
    
            call grib_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface)

            if ( discipline_idx(l)       .eq. iv_discpl .and. &
                 parameterCategory_idx(j).eq. iv_paramC .and. &
                 parameterNumber_idx(k)  .eq. iv_paramN .and. &
                 typeOfFirstFixedSurface .eq. iv_typeSf) then
  
              call grib_get(igrib,'numberOfPoints',numberOfPoints)
              call grib_get(igrib,'Ni',Ni)
              call grib_get(igrib,'Nj',Nj)
              if(nx_fullmet.ne.Ni)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"nx_fullmet = ",nx_fullmet
                write(MR_global_error,*)"Ni         = ",Ni
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"ny_fullmet = ",ny_fullmet
                write(MR_global_error,*)"Nj         = ",Nj
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
                call grib_get(igrib,'values',values)
                do m = 1,Nj
                  rstrt = (m-1)*Ni + 1
                  rend  = m*Ni
                  slice(1:Ni,m) = values(rstrt:rend)
                enddo
                deallocate(values)
        
               ! There is no guarentee that grib levels are in order so...
               ! Now loop through the pressure values for this variable and put this
               ! slice at the correct level.
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.level_idx(i))then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
             endif
        
                      call grib_release(igrib)
                      call grib_new_from_index(idx,igrib, iret)
                    enddo
                    call grib_release(igrib)
        
                  enddo ! loop on forecastTime
                enddo ! loop on level
              enddo ! loop on parameterNumber
            enddo ! loop on parameterCategory
          enddo ! loop on discipline
        
          call grib_index_release(idx)
        else
          ! We don't have/(can't make) the index file so scan all messages of the
          ! grib2 file
          write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                    trim(adjustl(grib_file_path))
          ifile=5
          call grib_open_file(ifile,grib_file_path,'R')
        
          !     turn on support for multi fields messages */
          call grib_multi_support_on()
        
          ! Loop on all the messages in a file.
          call grib_new_from_file(ifile,igrib,iret)
          count1=0
          do while (iret/=GRIB_END_OF_FILE)
            count1=count1+1
            call grib_get(igrib,'discipline',              grb_discipline)
            call grib_get(igrib,'parameterCategory',       grb_parameterCategory)
            call grib_get(igrib,'parameterNumber',         grb_parameterNumber)
            call grib_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface)
            call grib_get(igrib,'scaledValueOfFirstFixedSurface',grb_level)
            call grib_get(igrib,'scaledValueOfFirstFixedSurface',grb_scaledValueOfFirstFixedSurface)
 
            if ( grb_discipline              .eq. iv_discpl .and. &
                 grb_parameterCategory       .eq. iv_paramC .and. &
                 grb_parameterNumber         .eq. iv_paramN .and. &
                 typeOfFirstFixedSurface     .eq. iv_typeSf) then
              call grib_get(igrib,'numberOfPoints',numberOfPoints)
              call grib_get(igrib,'Ni',Ni)
              call grib_get(igrib,'Nj',Nj)
              if(nx_fullmet.ne.Ni)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"nx_fullmet = ",nx_fullmet
                write(MR_global_error,*)"Ni         = ",Ni
                stop 1
              endif
              if(ny_fullmet.ne.Nj)then
                write(MR_global_error,*)"MR ERROR:  Grid is not the expected size"
                write(MR_global_error,*)"ny_fullmet = ",ny_fullmet
                write(MR_global_error,*)"Nj         = ",Nj
                stop 1
              endif
              allocate(values(numberOfPoints))
              allocate(slice(Ni,Nj))
              call grib_get(igrib,'values',values)
              do m = 1,Nj
                rstrt = (m-1)*Ni + 1
                rend  = m*Ni
                slice(1:Ni,m) = values(rstrt:rend)
              enddo
              deallocate(values)
  
               ! There is no guarentee that grib levels are in order so...
               ! Now loop through the pressure values for this variable and put
               ! this slice at the correct level.
               if(ivar.eq.16)grb_level = grb_scaledValueOfFirstFixedSurface
               do kk = 1,np_met_loc
                 if(p_met_loc(kk).eq.grb_level)then
                   full_values(:,:,kk) = real(slice(:,:),kind=sp)
                   exit
                 endif
               enddo
               deallocate(slice)
            endif
            call grib_release(igrib)
            call grib_new_from_file(ifile,igrib, iret)
          enddo
          call grib_close_file(ifile)
        endif ! Use_GRIB_index
      endif ! MR_GRIB_Version eq 1 or 2

      if(Dimension_of_Variable.eq.3)then
        MR_dum3d_metP = 0.0_sp
        allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))

        do i=1,ict        !read subgrid at current time step
            ! for any other 3d variable (non-WRF, non-NCEP)
          temp3d_sp(ileft(i):iright(i)              ,1:ny_submet            ,1:np_met_loc,1) = &
          full_values(iistart(i):iistart(i)+iicount(i)-1,jstart:jstart+ny_submet-1,1:np_met_loc)
        enddo

          do j=1,ny_submet
            itmp = ny_submet-j+1
            !reverse the j indices (since they increment from N to S)
            if(y_inverted)then
              MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,itmp,1:np_met_loc,1)
            else
              MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,j,1:np_met_loc,1)
            endif
          enddo

        deallocate(temp3d_sp)

      elseif(Dimension_of_Variable.eq.2)then
!        if(IsCatagorical)then
!        else
          allocate(temp2d_sp(nx_submet,ny_submet,1))
          if(ivar.eq.11.or.ivar.eq.12)then
              ! Surface winds usually have a z coordinate as well
            allocate(temp3d_sp(nx_submet,ny_submet,1,1))
          endif
  
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then

            else
              ! 2d variables for iwf .ne. 25
!              if(ivar.eq.11.or.ivar.eq.12)then
!                ! Surface velocities do have a z dimension
!                nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
!                         start = (/iistart(i),jstart,1,iwstep/),       &
!                         count = (/iicount(i),ny_submet,1,1/))
!                if(nSTAT.ne.0)then
!                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   stop 1
!                endif
!                do j=1,ny_submet
!                  itmp = ny_submet-j+1
!                  if(y_inverted)then
!                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,itmp,1,1)
!                  else
!                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,j,1,1)
!                  endif
!                enddo
!              else
!                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(ileft(i):iright(i),:,:), &
!                         start = (/iistart(i),jstart,iwstep/),       &
!                         count = (/iicount(i),ny_submet,1/))
!                if(nSTAT.ne.0)then
!                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
!                   stop 1
!                endif
          temp2d_sp(ileft(i):iright(i)              ,1:ny_submet,1) = &
        full_values(iistart(i):iistart(i)+iicount(i)-1,jstart:jstart+ny_submet-1,1)

                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,itmp,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,j,1)
                  endif
                enddo
!              endif
            endif
          enddo
          deallocate(temp2d_sp)
          if(ivar.eq.11.or.ivar.eq.12) deallocate(temp3d_sp)
!        endif ! IsCatagorical
      endif ! Dimension_of_Variable.eq.2

      if(ivar.eq.1)then
        ! If this is filling HGT, then we need to do a special QC check
        !if(MR_iwindformat.eq.24)then
          ! It seems like only NASA has NaNs for pressures greater than surface
          ! pressure
          do i=1,nx_submet
            do j=1,ny_submet
              do k=1,np_met_loc
                if(MR_dum3d_metP(i,j,k).gt.1.0e10_sp)then
                   ! linearly interpolate in z
                   ! find the first non NaN above k
                   do kk = k+1,np_met_loc,1
                     if(MR_dum3d_metP(i,j,kk).lt.1.0e10_sp)exit
                   enddo
                   if(kk.eq.np_met_loc+1)then
                     kk=np_met_loc
                     MR_dum3d_metP(i,j,kk) = 0.0_sp
                   endif
                   ! find the first non NaN below k if k!=1
                   do kkk = max(k-1,1),1,-1
                     if(MR_dum3d_metP(i,j,kkk).lt.1.0e10_sp)exit
                   enddo
                   if(kkk.eq.0)then
                     kkk=1
                     MR_dum3d_metP(i,j,kkk) = 0.0_sp
                   endif
                   MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,kkk) + &
                         (MR_dum3d_metP(i,j,kk)-MR_dum3d_metP(i,j,kkk)) * &
                         real(k-kkk,kind=sp)/real(kk-kkk,kind=sp)
                endif
              enddo
            enddo
          enddo
        !endif
        ! convert m to km
        MR_dum3d_metP = MR_dum3d_metP / 1000.0_sp
      elseif(Dimension_of_Variable.eq.3)then
        ! Do QC checking of all other 3d variables
        if(ivar.eq.2.or.ivar.eq.3.or.ivar.eq.4)then
          ! taper winds (vx,vy,vz) to zero at ground surface
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=0.0_sp)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=0.0_sp)
          endif
        elseif(ivar.eq.5)then
          ! set ground and top-level conditions for temperature
          Z_top = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/real(100.0,kind=sp))
          T_top = MR_Temp_US_StdAtm(Z_top)
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          endif
        else
          ! For other variables, use the top and bottom non-fill values
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat))
          else
            call MR_QC_3dvar(ivar,nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat))
          endif
        endif
      endif

      if(ivar.eq.4)then
          ! For pressure vertical velocity, convert from Pa s to m/s by dividing
          ! by pressure gradient
        idx = Met_var_zdim_idx(ivar)
        do k=1,np_met_loc
          do i=1,nx_submet
            do j=1,ny_submet
              if(k.eq.1)then
                ! Use one-sided gradients for bottom
                !del_P = p_fullmet_Vz_sp(2)-p_fullmet_Vz_sp(1)
                del_p = levs_fullmet_sp(idx,2) - levs_fullmet_sp(idx,1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,2) - MR_geoH_metP_last(i,j,1)
                else
                  del_H = MR_geoH_metP_next(i,j,2) - MR_geoH_metP_next(i,j,1)
                endif
              elseif(k.eq.np_met_loc)then
                ! Use one-sided gradients for top
                !del_P = p_fullmet_Vz_sp(np_met_loc) - &
                !         p_fullmet_Vz_sp(np_met_loc-1)
                del_p = levs_fullmet_sp(idx,nlevs_fullmet(idx)) - &
                        levs_fullmet_sp(idx,nlevs_fullmet(idx)-1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,np_met_loc) - &
                           MR_geoH_metP_last(i,j,np_met_loc-1)
                else
                  del_H = MR_geoH_metP_next(i,j,np_met_loc) - &
                           MR_geoH_metP_next(i,j,np_met_loc-1)
                endif
              else
                ! otherwise, two-sided calculation
                !del_P = p_fullmet_Vz_sp(k+1)-p_fullmet_Vz_sp(k-1)
                del_p = levs_fullmet_sp(idx,k+1) - &
                        levs_fullmet_sp(idx,k)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,k+1) - MR_geoH_metP_last(i,j,k-1)
                else
                  del_H = MR_geoH_metP_next(i,j,k+1) - MR_geoH_metP_next(i,j,k-1)
                endif
              endif
              del_h = del_H * 1000.0_sp ! convert to m
              if(abs(del_H).gt.MR_EPS_SMALL)then
                dpdz  = del_P/del_H
              else
                write(MR_global_error,*)'MR ERROR: failed to calculate dpdz'
                write(MR_global_error,*)i,j,k,del_P,del_H
                write(MR_global_error,*)MR_geoH_metP_last(i,j,:)
                stop 1
              endif
              MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,k) / dpdz
            enddo
          enddo
        enddo
      endif
      MR_dum3d_metP = MR_dum3d_metP * Met_var_conversion_factor(ivar)

      end subroutine MR_Read_MetP_Variable_GRIB

