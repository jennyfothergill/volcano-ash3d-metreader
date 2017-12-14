!      subroutine MR_Set_Gen_Index_GRIB(grib2_file)
!
!      use MetReader
!      use grib_api
!
!      implicit none
!
!      character(len=130)  :: grib2_file
!
!      character(len=130)  :: index_file
!
!      integer            :: idx
!      integer            :: iret
!      integer            :: i,j,k,l,t
!      integer            :: count1=0
!        ! Used for keys
!      integer(kind=4),dimension(:),allocatable :: discipline_idx
!      integer(kind=4),dimension(:),allocatable :: parameterCategory_idx
!      integer(kind=4),dimension(:),allocatable :: parameterNumber_idx
!      integer(kind=4),dimension(:),allocatable :: level_idx
!      integer(kind=4),dimension(:),allocatable :: forecastTime_idx
!      integer(kind=4)  :: disciplineSize
!      integer(kind=4)  :: parameterCategorySize
!      integer(kind=4)  :: parameterNumberSize
!      integer(kind=4)  :: levelSize
!      integer(kind=4)  :: forecastTimeSize
!
!      ! Message identifier.
!      integer            :: igrib
!
!! discipline
!! 0 Meteorological products
!! 1 Hydrological products
!! 2 Land surface products
!!10 Oceanographic products
!
!! parameterCategory
!! 0 Temperature
!! 1 Moisture
!! 2 Momentum
!! 3 Mass
!! 4 Short-wave Radiation
!! 5 Long-wave Radiation
!! 6 Cloud
!! 7 Thermodynamic Stability indices
!
!! parameterNumber
!!  This is a sub-category of the discipline and parameterCategory
!
!! typeOfFirstFixedSurface
!!  1 Ground or water surface
!!  2 Cloud base level
!!  3 Level of cloud tops
!!100 Isobaric surface  (Pa)
!!103 Specified height level above ground  (m)
!!106 Depth below land surface  (m)
!!200 Unknown code table entry
!! http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table4-1.shtml
!!    shortName 
!!            discipline 
!!              parameterCategory 
!!                  parameterNumber
!! 1  gh      0 3   5 100   HGT   Geopotential_height_isobaric
!! 2  u       0 2   2 100  UGRD   u-component_of_wind_isobaric
!! 3  v       0 2   3 100  VGRD   v-component_of_wind_isobaric
!! 4  w       0 2   8 100  VVEL   Vertical_velocity_pressure_isobaric
!! 5  t       0 0   0 100   TMP   Temperature_isobaric
!! 10 hpbl    0 3 196   1  HPBL   Planetary_Boundary_Layer_Height_surface
!! 11 u       0 2   2 103  UGRD   u-component_of_wind_height_above_ground
!! 12 v       0 2   3 103  VGRD   v-component_of_wind_height_above_ground
!! 13 fricv   0 2 197   1 FRICV   Frictional_Velocity_surface
!! 15 sd      0 1  11   1  SNOD   Snow_depth_surface
!! 16 soilw   2 0 192 106 SOILW   Volumetric_Soil_Moisture_Content_depth_below_surface_layer
!! 17 sr      2 0   1   1  SFCR   Surface_roughness_surface
!! 18 gust    0 2  22   1  GUST   Wind_speed_gust_surface
!! 20 pres    0 3   0   2  PRES   Pressure_cloud_base
!! 21 pres    0 3   0   3  PRES   Pressure_cloud_tops
!! 23 tcc     0 6   1 200  TCDC   Total_cloud_cover_entire_atmosphere
!! 30 r       0 1   1 100    RH   Relative_humidity_isobaric
!! 31 q       0 1   0 100  SPFH   Specific_humidity_isobaric
!! 32 clwmr   0 1  22 100 CLWMR   Cloud_mixing_ratio_isobaric
!! 33 snmr    0 1  25 100  SNMR   Snow_mixing_ratio_isobaric
!! 40 crain   0 1 192   1 CRAIN   Categorical_Rain_surface
!! 41 csnow   0 1 195   1 CSNOW   Categorical_Snow_surface
!! 42 cfrzr   0 1 193   1 CFRZR   Categorical_Freezing_Rain_surface
!! 43 cicep   0 1 194   1 CICEP   Categorical_Ice_Pellets_surface
!! 44 prate   0 1   7   1 PRATE   Precipitation_rate_surface
!! 45                   1 CPRAT   Convective_Precipitation_Rate_surface
!! 45 tp      0 1   8   APCP   Total_precipitation_surface_0_Hour_Accumulation
!! 46                  ACPCP
!! Convective_precipitation_surface_0_Hour_Accumulation
!! 47 ncpcp   0 1   9  NCPCP
!! Large-scale_precipitation_non-convective_surface_0_Hour_Accumulation
!
!      index_file = adjustl(trim(grib2_file)) // ".index"
!      write(MR_global_info,*)"Generating index file: ",index_file
!
!          ! create an index from a grib file using some keys
!        call grib_index_create(idx,adjustl(trim(grib2_file)),&
!              'discipline,parameterCategory,parameterNumber,level,forecastTime')
!        !  call grib_index_read(idx,index_file)
!        call grib_multi_support_on()
!      
!          ! get the number of distinct values of all the keys in the index
!        call grib_index_get_size(idx,'discipline',disciplineSize)
!        call grib_index_get_size(idx,'parameterCategory',parameterCategorySize)
!        call grib_index_get_size(idx,'parameterNumber',parameterNumberSize)
!        call grib_index_get_size(idx,'level',levelSize)
!        call grib_index_get_size(idx,'forecastTime',forecastTimeSize)
!      
!          ! allocate the arry to contain the list of distinct values
!        allocate(discipline_idx(disciplineSize))
!        allocate(parameterCategory_idx(parameterCategorySize))
!        allocate(parameterNumber_idx(parameterNumberSize))
!        allocate(level_idx(levelSize))
!        allocate(forecastTime_idx(forecastTimeSize))
!      
!          ! get the list of distinct key values from the index
!        call grib_index_get(idx,'discipline',discipline_idx)
!        call grib_index_get(idx,'parameterCategory',parameterCategory_idx)
!        call grib_index_get(idx,'parameterNumber',parameterNumber_idx)
!        call grib_index_get(idx,'level',level_idx)
!        call grib_index_get(idx,'forecastTime',forecastTime_idx)
!      
!        count1=0
!        do l=1,disciplineSize
!          call grib_index_select(idx,'discipline',discipline_idx(l))
!      
!          do j=1,parameterCategorySize
!            call grib_index_select(idx,'parameterCategory',parameterCategory_idx(j))
!      
!            do k=1,parameterNumberSize
!              call grib_index_select(idx,'parameterNumber',parameterNumber_idx(k))
!      
!              do i=1,levelSize
!                call grib_index_select(idx,'level',level_idx(i))
!      
!                do t=1,forecastTimeSize
!                  call grib_index_select(idx,'forecastTime',forecastTime_idx(t))
!      
!      
!                call grib_new_from_index(idx,igrib, iret)
!                do while (iret /= GRIB_END_OF_INDEX)
!                   count1=count1+1
!                   !call grib_get(igrib,'shortName',sName)
!                   !write(MR_global_info,*)count1,sName
!                   call grib_release(igrib)
!                   call grib_new_from_index(idx,igrib, iret)
!                enddo
!                call grib_release(igrib)
!      
!                enddo ! loop on forecastTime
!              enddo ! loop on level
!            enddo ! loop on parameterNumber
!          enddo ! loop on parameterCategory
!        enddo ! loop on discipline
!      
!        call grib_index_write(idx,adjustl(trim(index_file)))
!      
!        call grib_index_release(idx)
!
!      end subroutine MR_Set_Gen_Index_GRIB


!##############################################################################
!##############################################################################
!##############################################################################
!    call grib_index_read(idx,index_file)
!
!    ! get the number of distinct values of all the keys in the index
!  call grib_index_get_size(idx,'discipline',disciplineSize)
!  call grib_index_get_size(idx,'parameterCategory',parameterCategorySize)
!  call grib_index_get_size(idx,'parameterNumber',parameterNumberSize)
!  call grib_index_get_size(idx,'level',levelSize)
!  call grib_index_get_size(idx,'forecastTime',forecastTimeSize)
!
!    ! allocate the arry to contain the list of distinct values
!  allocate(discipline_idx(disciplineSize))
!  allocate(parameterCategory_idx(parameterCategorySize))
!  allocate(parameterNumber_idx(parameterNumberSize))
!  allocate(level_idx(levelSize))
!  allocate(forecastTime_idx(forecastTimeSize))
!
!    ! get the list of distinct key values from the index
!  call grib_index_get(idx,'discipline',discipline_idx)
!  call grib_index_get(idx,'parameterCategory',parameterCategory_idx)
!  call grib_index_get(idx,'parameterNumber',parameterNumber_idx)
!  call grib_index_get(idx,'level',level_idx)
!  call grib_index_get(idx,'forecastTime',forecastTime_idx)
!
!  count1=0
!  do l=1,disciplineSize
!    call grib_index_select(idx,'discipline',discipline_idx(l))
!
!    do j=1,parameterCategorySize
!      call grib_index_select(idx,'parameterCategory',parameterCategory_idx(j))
!
!      do k=1,parameterNumberSize
!        call grib_index_select(idx,'parameterNumber',parameterNumber_idx(k))
!
!        do i=1,levelSize
!          call grib_index_select(idx,'level',level_idx(i))
!
!          do t=1,forecastTimeSize
!            call grib_index_select(idx,'forecastTime',forecastTime_idx(t))
!
!
!          call grib_new_from_index(idx,igrib, iret)
!          do while (iret /= GRIB_END_OF_INDEX)
!             count1=count1+1
!
!    call grib_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface)
!
!    iv_discpl = 0
!    iv_paramC = 3
!    iv_paramN = 5
!    iv_typeSf = 100
!    if ( discipline_idx(l)              .eq. iv_discpl .and. &
!         parameterCategory_idx(j)       .eq. iv_paramC .and. &
!         parameterNumber_idx(k)         .eq. iv_paramN .and. &
!         typeOfFirstFixedSurface .eq. iv_typeSf) then
!      call grib_get(igrib,'numberOfPoints',numberOfPoints)
!      call grib_get(igrib,'Ni',Ni)
!      call grib_get(igrib,'Nj',Nj)
!      allocate(values(numberOfPoints))
!      allocate(slice(Ni,Nj))
!
!      call grib_get(igrib,'shortName',sName)
!      call grib_get(igrib,'level',lev)
!      write(MR_global_info,*)sName,lev
!      if(lev.eq.1000)then
!        call grib_get(igrib,'values',values)
!        do m = 1,Nj
!          istrt = (m-1)*Ni + 1
!          iend  = m*Ni
!          slice(1:Ni,m) = values(istrt:iend)
!        enddo
!        write(MR_global_info,*)slice(:,1)
!      endif
!      deallocate(values)
!      deallocate(slice)
!    endif
!             call grib_release(igrib)
!             call grib_new_from_index(idx,igrib, iret)
!          enddo
!          call grib_release(igrib)
!
!          enddo ! loop on forecastTime
!        enddo ! loop on level
!      enddo ! loop on parameterNumber
!    enddo ! loop on parameterCategory
!  enddo ! loop on discipline
!
!  call grib_index_release(idx)



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
!       Met_dim_names, Met_var_names, Met_var_conversion_factor, Met_var_IsAvailable
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

      integer :: i, k
      real(kind=sp) :: xLL_fullmet
      real(kind=sp) :: yLL_fullmet
      real(kind=sp) :: xUR_fullmet
      real(kind=sp) :: yUR_fullmet

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Read_Met_DimVars_GRIB                    ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(MR_iwindformat.ne.0)then
        Met_dim_names = ""
        Met_var_names = ""
        Met_var_conversion_factor(:)  = 1.0_sp
      endif

      ! Initialize dimension and variable names
      if(MR_iwindformat.eq.0)then
          ! This expects that MR_iwf_template has been filled by the calling program
        call MR_Read_Met_Template
      elseif(MR_iwindformat.eq.2)then
        ! This is reserved for reading Radiosonde data
      elseif(MR_iwindformat.eq.3)then
          ! NARR3D NAM221 50 km North America files
          ! Note that winds are "earth-relative" and must be rotated!
          ! See
          ! http://www.emc.ncep.noaa.gov/mmb/rreanl/faq.html#eta-winds
      elseif(MR_iwindformat.eq.4)then
          ! NARR3D NAM221 50 km North America files (RAW : assumes full set of
          ! variables)
          ! Note that winds are "earth-relative" and must be rotated!
          ! See
          ! http://www.emc.ncep.noaa.gov/mmb/rreanl/faq.html#eta-winds
      elseif(MR_iwindformat.eq.5)then
          ! NAM216 AK 45km
      elseif(MR_iwindformat.eq.6)then
          ! NAM Regional 90 km grid 104
      elseif(MR_iwindformat.eq.7)then
          ! CONUS 212 40km
      elseif(MR_iwindformat.eq.8)then
          ! CONUS 218 (12km)
          ! wget
          ! ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.20121231/nam.t00z.awphys${fh[$i]}.grb2.tm00
      elseif(MR_iwindformat.eq.9)then
         !Unassigned
      elseif(MR_iwindformat.eq.11)then
          ! NAM 196 2.5 km HI

        Met_dim_IsAvailable=.false.
        Met_var_IsAvailable=.false.

        Met_dim_names(1) = "time"      ! time
          Met_dim_IsAvailable(1)=.true. 
        Met_dim_names(2) = "isobaric2"  ! pressure
          Met_dim_IsAvailable(2)=.true. 
        Met_dim_names(3) = "y"         ! y
          Met_dim_IsAvailable(3)=.true. 
        Met_dim_names(4) = "x"         ! x
          Met_dim_IsAvailable(4)=.true. 
        Met_dim_names(5) = "isobaric2"  ! pressure coordinate for Vz
          Met_dim_IsAvailable(5)=.true. 
        Met_dim_names(6) = "isobaric2" ! pressure coordinate for RH
          Met_dim_IsAvailable(6)=.true. 

        ! Mechanical / State variables
        Met_var_names(1)               = "gh"
          Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/)
          Met_var_IsAvailable(1)       = .true.
        Met_var_names(2)               = "u"
          Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
          Met_var_IsAvailable(2)       = .true.
        Met_var_names(3)               = "v"
          Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
          Met_var_IsAvailable(3)       = .true.
        Met_var_names(4)               = "w"
          Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
          Met_var_IsAvailable(4)       = .true.
        Met_var_names(5)               = "t"
          Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
          Met_var_IsAvailable(5)       =.true.

        ! Surface
        Met_var_names(10)              = "hpbl"
          Met_var_GRIB2_DPcPnSt(10,1:4)= (/0, 3, 196, 1/)
          Met_var_IsAvailable(10)      = .true.
        Met_var_names(11)              = "u"
          Met_var_GRIB2_DPcPnSt(11,1:4)= (/0, 2, 2, 103/)
          Met_var_IsAvailable(11)      = .true.
        Met_var_names(12)              = "v"
          Met_var_GRIB2_DPcPnSt(12,1:4)= (/0, 2, 3, 103/)
          Met_var_IsAvailable(12)      = .true.
        Met_var_names(13)              = "fricv"
          Met_var_GRIB2_DPcPnSt(13,1:4)= (/0, 2, 197, 1/)
          Met_var_IsAvailable(13)      = .true.
        !  14 = Displacement Height
        Met_var_names(15)              = "sd"
          Met_var_GRIB2_DPcPnSt(15,1:4)= (/0, 1, 11, 1/)
          Met_var_IsAvailable(15)      = .true.
        Met_var_names(16)              = "soilw"
          Met_var_GRIB2_DPcPnSt(16,1:4)= (/2, 0, 192, 106/)
          Met_var_IsAvailable(16)      = .true.
        Met_var_names(17)              = "sr"
          Met_var_GRIB2_DPcPnSt(17,1:4)= (/2, 0, 1, 1/)
          Met_var_IsAvailable(17)      = .true.
        Met_var_names(18)              = "gust"
          Met_var_GRIB2_DPcPnSt(18,1:4)= (/0, 2, 22, 1/)
          Met_var_IsAvailable(18)      = .true.
        ! Atmospheric Structure
        Met_var_names(20)              = "pres"
          Met_var_GRIB2_DPcPnSt(20,1:4)= (/0, 3, 0, 2/)
          Met_var_IsAvailable(20)      = .true.
        Met_var_names(21)              = "pres"
          Met_var_GRIB2_DPcPnSt(21,1:4)= (/0, 3, 0, 3/)
          Met_var_IsAvailable(21)      = .true.
        Met_var_names(23)              = "tcc"
          Met_var_GRIB2_DPcPnSt(23,1:4)= (/0, 6, 1, 200/)
          Met_var_IsAvailable(23)      = .true.
        ! Moisture
        Met_var_names(30)              = "r"
          Met_var_GRIB2_DPcPnSt(30,1:4)= (/0, 1, 1, 100/)
          Met_var_IsAvailable(30)      = .true.
        Met_var_names(31)              = "q"
          Met_var_GRIB2_DPcPnSt(31,1:4) = (/0, 1, 0, 100/)
          Met_var_IsAvailable(31)      = .true.
        Met_var_names(32)              = "clwmr"
          Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
          Met_var_IsAvailable(32)      = .true.
        Met_var_names(33)              = "snmr"
          Met_var_GRIB2_DPcPnSt(33,1:4)= (/0, 1, 25, 100/)
          Met_var_IsAvailable(33)      = .true.
        ! Precipitation
        Met_var_names(40)              = "crain"
          Met_var_GRIB2_DPcPnSt(40,1:4)= (/0, 1, 192, 1/)
          Met_var_IsAvailable(40)      = .true.
        Met_var_names(41)              = "csnow"
          Met_var_GRIB2_DPcPnSt(41,1:4)= (/0, 1, 195, 1/)
          Met_var_IsAvailable(41)      = .true.
        Met_var_names(42)              = "cfrzr"
          Met_var_GRIB2_DPcPnSt(42,1:4)= (/0, 1, 193, 1/)
          Met_var_IsAvailable(42)      = .true.
        Met_var_names(43)              = "cicep"
          Met_var_GRIB2_DPcPnSt(43,1:4)= (/0, 1, 194, 1/)
          Met_var_IsAvailable(43)      = .true.
        Met_var_names(44)              = "prate"
          Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 7, 1/)
          Met_var_IsAvailable(44)      = .true.

        fill_value_sp(MR_iwindformat) = -9999.0_sp ! actually NaNf
        MR_ForecastInterval = 1.0_4

      elseif(MR_iwindformat.eq.12)then
          ! NAM 198 5.953 km AK

        Met_dim_IsAvailable=.false.
        Met_var_IsAvailable=.false.  

        Met_dim_names(1) = "time"      ! time
          Met_dim_IsAvailable(1)=.true.
        Met_dim_names(2) = "isobaric2"  ! pressure
          Met_dim_IsAvailable(2)=.true.
        Met_dim_names(3) = "y"         ! y
          Met_dim_IsAvailable(3)=.true.
        Met_dim_names(4) = "x"         ! x
          Met_dim_IsAvailable(4)=.true.
        Met_dim_names(5) = "isobaric2"  ! pressure coordinate for Vz
          Met_dim_IsAvailable(5)=.true.
        Met_dim_names(6) = "isobaric2" ! pressure coordinate for RH
          Met_dim_IsAvailable(6)=.true.

        ! Mechanical / State variables
        Met_var_names(1)               = "gh"
          Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/)
          Met_var_IsAvailable(1)       = .true.
        Met_var_names(2)               = "u"
          Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
          Met_var_IsAvailable(2)       = .true.
        Met_var_names(3)               = "v"
          Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
          Met_var_IsAvailable(3)       = .true.
        Met_var_names(4)               = "w"
          Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
          Met_var_IsAvailable(4)       = .true.
        Met_var_names(5)               = "t"
          Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
          Met_var_IsAvailable(5)       =.true.

        ! Surface
        Met_var_names(10)              = "hpbl"
          Met_var_GRIB2_DPcPnSt(10,1:4)= (/0, 3, 196, 1/)
          Met_var_IsAvailable(10)      = .true.
        Met_var_names(11)              = "u"
          Met_var_GRIB2_DPcPnSt(11,1:4)= (/0, 2, 2, 103/)
          Met_var_IsAvailable(11)      = .true.
        Met_var_names(12)              = "v"
          Met_var_GRIB2_DPcPnSt(12,1:4)= (/0, 2, 3, 103/)
          Met_var_IsAvailable(12)      = .true.
        Met_var_names(13)              = "fricv"
          Met_var_GRIB2_DPcPnSt(13,1:4)= (/0, 2, 197, 1/)
          Met_var_IsAvailable(13)      = .true.
        !  14 = Displacement Height
        Met_var_names(15)              = "sd"
          Met_var_GRIB2_DPcPnSt(15,1:4)= (/0, 1, 11, 1/)
          Met_var_IsAvailable(15)      = .true.
        Met_var_names(16)              = "soilw"
          Met_var_GRIB2_DPcPnSt(16,1:4)= (/2, 0, 192, 106/)
          Met_var_IsAvailable(16)      = .true.
        Met_var_names(17)              = "sr"
          Met_var_GRIB2_DPcPnSt(17,1:4)= (/2, 0, 1, 1/)
          Met_var_IsAvailable(17)      = .true.
        Met_var_names(18)              = "gust"
          Met_var_GRIB2_DPcPnSt(18,1:4)= (/0, 2, 22, 1/)
          Met_var_IsAvailable(18)      = .true.
        ! Atmospheric Structure
        Met_var_names(20)              = "pres"
          Met_var_GRIB2_DPcPnSt(20,1:4)= (/0, 3, 0, 2/)
          Met_var_IsAvailable(20)      = .true.
        Met_var_names(21)              = "pres"
          Met_var_GRIB2_DPcPnSt(21,1:4)= (/0, 3, 0, 3/)
          Met_var_IsAvailable(21)      = .true.
        Met_var_names(23)              = "tcc"
          Met_var_GRIB2_DPcPnSt(23,1:4)= (/0, 6, 1, 200/)
          Met_var_IsAvailable(23)      = .true.
        ! Moisture
        Met_var_names(30)              = "r"
          Met_var_GRIB2_DPcPnSt(30,1:4)= (/0, 1, 1, 100/)
          Met_var_IsAvailable(30)      = .true.
        Met_var_names(31)              = "q"
          Met_var_GRIB2_DPcPnSt(31,1:4) = (/0, 1, 0, 100/)
          Met_var_IsAvailable(31)      = .true.
        Met_var_names(32)              = "clwmr"
          Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
          Met_var_IsAvailable(32)      = .true.
        Met_var_names(33)              = "snmr"
          Met_var_GRIB2_DPcPnSt(33,1:4)= (/0, 1, 25, 100/)
          Met_var_IsAvailable(33)      = .true.
        ! Precipitation
        Met_var_names(40)              = "crain"
          Met_var_GRIB2_DPcPnSt(40,1:4)= (/0, 1, 192, 1/)
          Met_var_IsAvailable(40)      = .true.
        Met_var_names(41)              = "csnow"
          Met_var_GRIB2_DPcPnSt(41,1:4)= (/0, 1, 195, 1/)
          Met_var_IsAvailable(41)      = .true.
        Met_var_names(42)              = "cfrzr"
          Met_var_GRIB2_DPcPnSt(42,1:4)= (/0, 1, 193, 1/)
          Met_var_IsAvailable(42)      = .true.
        Met_var_names(43)              = "cicep"
          Met_var_GRIB2_DPcPnSt(43,1:4)= (/0, 1, 194, 1/)
          Met_var_IsAvailable(43)      = .true.
        Met_var_names(44)              = "prate"
          Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 7, 1/)
          Met_var_IsAvailable(44)      = .true.

        fill_value_sp(MR_iwindformat) = -9999.0_sp ! actually NaNf
        MR_ForecastInterval = 1.0_4

      elseif(MR_iwindformat.eq.13)then
          ! NAM 91 2.976 km AK (nam198 at twice the resolution)
          !    all the variables below should be the same as iwf12=n198

        Met_dim_IsAvailable=.false.
        Met_var_IsAvailable=.false.

        Met_dim_names(1) = "time"      ! time
          Met_dim_IsAvailable(1)=.true.
        Met_dim_names(2) = "isobaric"  ! pressure
          Met_dim_IsAvailable(2)=.true.
        Met_dim_names(3) = "y"         ! y
          Met_dim_IsAvailable(3)=.true.
        Met_dim_names(4) = "x"         ! x
          Met_dim_IsAvailable(4)=.true.
        Met_dim_names(5) = "isobaric"  ! pressure coordinate for Vz
          Met_dim_IsAvailable(5)=.true.
        Met_dim_names(6) = "isobaric" ! pressure coordinate for RH
          Met_dim_IsAvailable(6)=.true.

        ! Mechanical / State variables
        Met_var_names(1)               = "gh"
          Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/)
          Met_var_IsAvailable(1)       = .true.
        Met_var_names(2)               = "u"
          Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
          Met_var_IsAvailable(2)       = .true.
        Met_var_names(3)               = "v"
          Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
          Met_var_IsAvailable(3)       = .true.
        Met_var_names(4)               = "w"
          Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
          Met_var_IsAvailable(4)       = .true.
        Met_var_names(5)               = "t"
          Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
          Met_var_IsAvailable(5)       =.true.

        ! Surface
        Met_var_names(10)              = "hpbl"
          Met_var_GRIB2_DPcPnSt(10,1:4)= (/0, 3, 196, 1/)
          Met_var_IsAvailable(10)      = .true.
        Met_var_names(11)              = "u"
          Met_var_GRIB2_DPcPnSt(11,1:4)= (/0, 2, 2, 103/)
          Met_var_IsAvailable(11)      = .true.
        Met_var_names(12)              = "v"
          Met_var_GRIB2_DPcPnSt(12,1:4)= (/0, 2, 3, 103/)
          Met_var_IsAvailable(12)      = .true.
        Met_var_names(13)              = "fricv"
          Met_var_GRIB2_DPcPnSt(13,1:4)= (/0, 2, 197, 1/)
          Met_var_IsAvailable(13)      = .true.
        !  14 = Displacement Height
        Met_var_names(15)              = "sd"
          Met_var_GRIB2_DPcPnSt(15,1:4)= (/0, 1, 11, 1/)
          Met_var_IsAvailable(15)      = .true.
        Met_var_names(16)              = "soilw"
          Met_var_GRIB2_DPcPnSt(16,1:4)= (/2, 0, 192, 106/)
          Met_var_IsAvailable(16)      = .true.
        Met_var_names(17)              = "sr"
          Met_var_GRIB2_DPcPnSt(17,1:4)= (/2, 0, 1, 1/)
          Met_var_IsAvailable(17)      = .true.
        Met_var_names(18)              = "gust"
          Met_var_GRIB2_DPcPnSt(18,1:4)= (/0, 2, 22, 1/)
          Met_var_IsAvailable(18)      = .true.
        ! Atmospheric Structure
        Met_var_names(20)              = "pres"
          Met_var_GRIB2_DPcPnSt(20,1:4)= (/0, 3, 0, 2/)
          Met_var_IsAvailable(20)      = .true.
        Met_var_names(21)              = "pres"
          Met_var_GRIB2_DPcPnSt(21,1:4)= (/0, 3, 0, 3/)
          Met_var_IsAvailable(21)      = .true.
        Met_var_names(23)              = "tcc"
          Met_var_GRIB2_DPcPnSt(23,1:4)= (/0, 6, 1, 200/)
          Met_var_IsAvailable(23)      = .true.
        ! Moisture
        Met_var_names(30)              = "r"
          Met_var_GRIB2_DPcPnSt(30,1:4)= (/0, 1, 1, 100/)
          Met_var_IsAvailable(30)      = .true.
        Met_var_names(31)              = "q"
          Met_var_GRIB2_DPcPnSt(31,1:4) = (/0, 1, 0, 100/)
          Met_var_IsAvailable(31)      = .true.
        Met_var_names(32)              = "clwmr"
          Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
          Met_var_IsAvailable(32)      = .true.
        Met_var_names(33)              = "snmr"
          Met_var_GRIB2_DPcPnSt(33,1:4)= (/0, 1, 25, 100/)
          Met_var_IsAvailable(33)      = .true.
        ! Precipitation
        Met_var_names(40)              = "crain"
          Met_var_GRIB2_DPcPnSt(40,1:4)= (/0, 1, 192, 1/)
          Met_var_IsAvailable(40)      = .true.
        Met_var_names(41)              = "csnow"
          Met_var_GRIB2_DPcPnSt(41,1:4)= (/0, 1, 195, 1/)
          Met_var_IsAvailable(41)      = .true.
        Met_var_names(42)              = "cfrzr"
          Met_var_GRIB2_DPcPnSt(42,1:4)= (/0, 1, 193, 1/)
          Met_var_IsAvailable(42)      = .true.
        Met_var_names(43)              = "cicep"
          Met_var_GRIB2_DPcPnSt(43,1:4)= (/0, 1, 194, 1/)
          Met_var_IsAvailable(43)      = .true.
        Met_var_names(44)              = "prate"
          Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 7, 1/)
          Met_var_IsAvailable(44)      = .true.

        fill_value_sp(MR_iwindformat) = -9999.0_sp ! actually NaNf
        MR_ForecastInterval = 1.0_4

      elseif (MR_iwindformat.eq.20.or.MR_iwindformat.eq.22) then
           ! GFS 0.5 (or 0.25) deg from http://www.nco.ncep.noaa.gov/pmb/products/gfs/
           ! or
           ! http://motherlode.ucar.edu/native/conduit/data/nccf/com/gfs/prod/
        Met_dim_IsAvailable=.false.
        Met_var_IsAvailable=.false.

        Met_dim_names(1) = "time"       ! time
          Met_dim_IsAvailable(1)=.true.
        Met_dim_names(2) = "isobaric3"   ! pressure (1.0e5-1.0e3 Pa or 1000 -> 10.0 hPa in 26 levels)
          Met_dim_IsAvailable(2)=.true.
        Met_dim_names(3) = "lat"        ! y        (90.0 -> -90.0)
          Met_dim_IsAvailable(3)=.true.
        Met_dim_names(4) = "lon"        ! x        (0.0 - 359.5)
          Met_dim_IsAvailable(4)=.true.
        Met_dim_names(5) = "isobaric4"  ! pressure coordinate for Vz (to 100 hPa in 21 levels)
          Met_dim_IsAvailable(5)=.true.
        Met_dim_names(6) = "isobaric3"  ! pressure coordinate for RH (to 10 hPa in 25 levels)
          Met_dim_IsAvailable(6)=.true.

        ! Momentum / State variables
        Met_var_names(1) = "Geopotential_height_isobaric"    ! float gpm
          Met_var_IsAvailable(1)=.true.
        Met_var_names(2) = "u-component_of_wind_isobaric"    ! float m/s
          Met_var_IsAvailable(2)=.true.
        Met_var_names(3) = "v-component_of_wind_isobaric"    ! float m/s
          Met_var_IsAvailable(3)=.true.
        Met_var_names(4) = "Vertical_velocity_pressure_isobaric" ! float Pa/s
          Met_var_IsAvailable(4)=.true.
        Met_var_names(5) = "Temperature_isobaric"            ! float deg K
          Met_var_IsAvailable(5)=.true.

        ! Surface
        Met_var_names(10) = "Planetary_Boundary_Layer_Height_surface"
          Met_var_IsAvailable(10)=.true.
        Met_var_names(11) = "U-component_of_wind_height_above_ground" !(time, height_above_ground, lat, lon) 10, 80, 100
          Met_var_IsAvailable(11)=.true.
        Met_var_names(12) = "V-component_of_wind_height_above_ground"
          Met_var_IsAvailable(12)=.true.
        ! Atmospheric Structure
        Met_var_names(22) = "Cloud_Top_Temperature"
          Met_var_IsAvailable(22)=.true.
        Met_var_names(23) = "Total_cloud_cover"
          Met_var_IsAvailable(23)=.true.
        ! Moisture
        Met_var_names(30) = "Relative_humidity_isobaric"      ! float percent
          Met_var_IsAvailable(30)=.true.
        Met_var_names(32) = "Cloud_mixing_ratio"     ! float kg/kg
          Met_var_IsAvailable(32)=.true.
        ! Precipitation
        Met_var_names(40) = "Categorical_Rain"
          Met_var_IsAvailable(40)=.true.
        Met_var_names(41) = "Categorical_Snow"
          Met_var_IsAvailable(41)=.true.
        Met_var_names(42) = "Categorical_Freezing_Rain"
          Met_var_IsAvailable(42)=.true.
        Met_var_names(43) = "Categorical_Ice_Pellets"
          Met_var_IsAvailable(43)=.true.
        Met_var_names(44) = "Precipitation_rate"     ! float liquid surface precip kg/m2/s
          Met_var_IsAvailable(44)=.true.
        Met_var_names(45) = "Convective_Precipitation_Rate"       ! float liquid convective precip kg/m2/s
          Met_var_IsAvailable(45)=.true.

        fill_value_sp(MR_iwindformat) = -9999.0_sp
        MR_ForecastInterval = 1.0_4

      elseif (MR_iwindformat.eq.21) then
           ! Old format GFS 0.5-degree
       elseif (MR_iwindformat.eq.23) then
         ! NCEP / doE reanalysis 2.5 degree files 
       elseif (MR_iwindformat.eq.24) then
         ! NASA-MERRA reanalysis 1.25 degree files 
       elseif (MR_iwindformat.eq.25) then
         ! NCEP/NCAR reanalysis 2.5 degree files 
       elseif (MR_iwindformat.eq.27) then
         ! NOAA-CIRES reanalysis 2.5 degree files 
       elseif (MR_iwindformat.eq.28) then
         ! ECMWF Interim Reanalysis (ERA-Interim)
       elseif (MR_iwindformat.eq.29) then
         ! JRA-55 reanalysis
       elseif (MR_iwindformat.eq.31) then
         ! Catania forecast
       elseif (MR_iwindformat.eq.32) then
         ! Air Force Weather Agency subcenter = 0
        Met_dim_IsAvailable=.false.
        Met_var_IsAvailable=.false.

        Met_dim_names(1) = "time"      ! time
          Met_dim_IsAvailable(1)=.true.
        Met_dim_names(2) = "isobaric3"  ! pressure
          Met_dim_IsAvailable(2)=.true.
        Met_dim_names(3) = "y"         ! y
          Met_dim_IsAvailable(3)=.true.
        Met_dim_names(4) = "x"         ! x
          Met_dim_IsAvailable(4)=.true.
        Met_dim_names(5) = "isobaric1"  ! pressure coordinate for Vz
          Met_dim_IsAvailable(5)=.true.
        Met_dim_names(6) = "isobaric3" ! pressure coordinate for RH
          Met_dim_IsAvailable(6)=.true.

        ! Mechanical / State variables
        Met_var_names(1)               = "gh"
          Met_var_GRIB2_DPcPnSt(1,1:4) = (/0, 3, 5, 100/)
          Met_var_IsAvailable(1)       = .true.
        Met_var_names(2)               = "u"
          Met_var_GRIB2_DPcPnSt(2,1:4) = (/0, 2, 2, 100/)
          Met_var_IsAvailable(2)       = .true.
        Met_var_names(3)               = "v"
          Met_var_GRIB2_DPcPnSt(3,1:4) = (/0, 2, 3, 100/)
          Met_var_IsAvailable(3)       = .true.
        Met_var_names(4)               = "w"
          Met_var_GRIB2_DPcPnSt(4,1:4) = (/0, 2, 8, 100/)
          Met_var_IsAvailable(4)       = .true.
        Met_var_names(5)               = "t"
          Met_var_GRIB2_DPcPnSt(5,1:4) = (/0, 0, 0, 100/)
          Met_var_IsAvailable(5)       =.true.

        ! Surface
        Met_var_names(10)              = "hpbl"
          !Met_var_GRIB2_DPcPnSt(10,1:4)= (/0, 3, 196, 1/)
          Met_var_GRIB2_DPcPnSt(10,1:4)= (/0, 3, 18, 1/)
          Met_var_IsAvailable(10)      = .true.
        Met_var_names(11)              = "u"
          Met_var_GRIB2_DPcPnSt(11,1:4)= (/0, 2, 2, 103/)
          Met_var_IsAvailable(11)      = .true.
        Met_var_names(12)              = "v"
          Met_var_GRIB2_DPcPnSt(12,1:4)= (/0, 2, 3, 103/)
          Met_var_IsAvailable(12)      = .true.
        Met_var_names(13)              = "fricv"
          !Met_var_GRIB2_DPcPnSt(13,1:4)= (/0, 2, 197, 1/)
          Met_var_GRIB2_DPcPnSt(13,1:4)= (/0, 2, 30, 1/)
          Met_var_IsAvailable(13)      = .true.

        !  14 = Displacement Height
        Met_var_names(15)              = "sd"
          !Met_var_GRIB2_DPcPnSt(15,1:4)= (/0, 1, 11, 1/)
          Met_var_GRIB2_DPcPnSt(15,1:4)= (/0, 1, 13, 1/)
          Met_var_IsAvailable(15)      = .true.
        Met_var_names(16)              = "soilw"
          !Met_var_GRIB2_DPcPnSt(16,1:4)= (/2, 0, 192, 106/)
          Met_var_GRIB2_DPcPnSt(16,1:4)= (/2, 3, 20, 106/)
          Met_var_IsAvailable(16)      = .true.
        Met_var_names(17)              = "sr"
          Met_var_GRIB2_DPcPnSt(17,1:4)= (/2, 0, 1, 1/)
          Met_var_IsAvailable(17)      = .true.
        Met_var_names(18)              = "gust"
          !Met_var_GRIB2_DPcPnSt(18,1:4)= (/0, 2, 22, 1/)
          Met_var_GRIB2_DPcPnSt(18,1:4)= (/0, 2, 22, 103/)
          Met_var_IsAvailable(18)      = .true.
        ! Atmospheric Structure
        Met_var_names(20)              = "pres"
          !Met_var_GRIB2_DPcPnSt(20,1:4)= (/0, 3, 0, 2/)
          Met_var_GRIB2_DPcPnSt(20,1:4)= (/0, 6, 11, 2/)
          Met_var_IsAvailable(20)      = .true.
        Met_var_names(21)              = "pres"
          !Met_var_GRIB2_DPcPnSt(21,1:4)= (/0, 3, 0, 3/)
          Met_var_GRIB2_DPcPnSt(21,1:4)= (/0, 6, 12, 3/)
          Met_var_IsAvailable(21)      = .true.
        Met_var_names(23)              = "tcc"
          !Met_var_GRIB2_DPcPnSt(23,1:4)= (/0, 6, 1, 200/)
          Met_var_GRIB2_DPcPnSt(23,1:4)= (/0, 6, 1, 10/)
          Met_var_IsAvailable(23)      = .true.
        ! Moisture
        Met_var_names(30)              = "r"
          Met_var_GRIB2_DPcPnSt(30,1:4)= (/0, 1, 1, 100/)
          Met_var_IsAvailable(30)      = .true.
        Met_var_names(32)              = "clwmr"
          !Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 22, 100/)
          Met_var_GRIB2_DPcPnSt(32,1:4) = (/0, 1, 2, 100/)
          Met_var_IsAvailable(32)      = .true.
        ! Precipitation
        Met_var_names(44)              = "prate"
          !Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 7, 1/)
          Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 54, 1/)
          Met_var_IsAvailable(44)      = .true.
        Met_var_names(45)              = "prate"
          Met_var_GRIB2_DPcPnSt(44,1:4)= (/0, 1, 37, 1/)
          Met_var_IsAvailable(44)      = .true.


        fill_value_sp(MR_iwindformat) = 9999.0_sp ! actually NaNf
        MR_ForecastInterval = 1.0_4


       elseif (MR_iwindformat.eq.40) then
         ! NASA-GEOS Cp
       elseif (MR_iwindformat.eq.41) then
         ! NASA-GEOS Np
      else
        ! Not a recognized MR_iwindformat
        ! call reading of custom windfile variable names
        write(MR_global_error,*)"MR ERROR : windfile format not recognized."
        stop 1
      endif

      if(Met_var_IsAvailable(4)) Have_Vz = .true.

      if(MR_iwindformat.eq.0)then
        ! Custom windfile (example for nam198)
        !  Need to populate
        !call MR_Set_Met_Dims_Custom_GRIB
        write(MR_global_error,*)"MR ERROR : Currently, grib reader only works for known grib2 files."
        write(MR_global_error,*)"           Custom reader forthcoming"
        stop 1
      elseif(MR_iwindformat.eq.2)then
        write(MR_global_error,*)"MR ERROR : MR_iwindformat = 2: should not be here."
        stop 1
      elseif(MR_iwindformat.eq.3.or.MR_iwindformat.eq.4)then
          ! 3 = NARR3D NAM221 32 km North America files
          ! 4 = RAW : assumes full set of variables
        call MR_Set_Met_NCEPGeoGrid(221)
        isGridRelative = .false.

        nt_fullmet = 1
        np_fullmet = 29
        np_fullmet_Vz = 29  ! omega
        np_fullmet_RH = 29  ! rhum
        np_fullmet_P0 = 1   ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, &
             550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, &
             300.0_sp, 275.0_sp, 250.0_sp, 225.0_sp, 200.0_sp, &
             175.0_sp, 150.0_sp, 125.0_sp, 100.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .false.

      elseif(MR_iwindformat.eq.5)then
        ! NAM 45-km Polar Sterographic
        call MR_Set_Met_NCEPGeoGrid(216)

        nt_fullmet = 29
        np_fullmet = 42
        np_fullmet_Vz = 39 ! omega
        np_fullmet_RH = 42  ! rhum

        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, & 
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.6)then
          ! 104 converted automatically from grib2
          ! NAM 90-km Polar Sterographic
        call MR_Set_Met_NCEPGeoGrid(104)

        nt_fullmet = 29
        np_fullmet = 39
        np_fullmet_Vz = 39 ! omega
        np_fullmet_RH = 19  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 950.0_sp, 900.0_sp, 850.0_sp, 800.0_sp, &
             750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
             500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
             250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.7)then
          ! 212 converted automatically from grib2
          ! CONUS 40-km Lambert Conformal
        call MR_Set_Met_NCEPGeoGrid(212)

        nt_fullmet = 29
        np_fullmet = 39
        np_fullmet_Vz = 39 ! omega
        np_fullmet_RH = 39  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.8)then
          !  12 KM CONUS)
        call MR_Set_Met_NCEPGeoGrid(218)

        nt_fullmet = 1
        np_fullmet = 39
        np_fullmet_Vz = 39 ! omega
        np_fullmet_RH = 39  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.9)then
          !  20 KM CONUS)
        call MR_Set_Met_NCEPGeoGrid(215)

        nt_fullmet = 1
        np_fullmet = 39
        np_fullmet_Vz = 39 ! omega
        np_fullmet_RH = 39  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.10)then
          ! NAM242 converted automatically from grib2
        call MR_Set_Met_NCEPGeoGrid(242)

        nt_fullmet = 29
        np_fullmet = 39
        np_fullmet_Vz = 29 ! omega
        np_fullmet_RH = 39  ! rhum
        np_fullmet_P0 = 1   ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp/)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
             250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.

      elseif(MR_iwindformat.eq.11)then
          ! NAM196 converted automatically from grib2
        call MR_Set_Met_NCEPGeoGrid(196)

        nt_fullmet = 1
        np_fullmet = 42
        np_fullmet_Vz = 42 ! omega
        np_fullmet_RH = 42  ! rhum
        np_fullmet_P0 = 1   ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.12)then
          ! NAM198 converted automatically from grib2
        call MR_Set_Met_NCEPGeoGrid(198)

        nt_fullmet = 1
        np_fullmet = 42
        np_fullmet_Vz = 42 ! omega
        np_fullmet_RH = 42  ! rhum
        np_fullmet_P0 = 1   ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.

      elseif(MR_iwindformat.eq.13)then
          ! NAM91 converted automatically from grib2
        call MR_Set_Met_NCEPGeoGrid(91)

        nt_fullmet = 1
        np_fullmet = 42
        np_fullmet_Vz = 42 ! omega
        np_fullmet_RH = 42  ! rhum
        np_fullmet_P0 = 1   ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  75.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.

      elseif(MR_iwindformat.eq.20)then
        ! GFS 0.5
        call MR_Set_Met_NCEPGeoGrid(4)

        nt_fullmet = 1
        np_fullmet = 26   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 21 ! omega
        np_fullmet_RH = 25  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp /)
        p_fullmet_Vz_sp(1:21) =  &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, & 
             100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, & 
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  10.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .true.

      elseif(MR_iwindformat.eq.21)then
        ! GFS 0.5 old style
        call MR_Set_Met_NCEPGeoGrid(4)

        nt_fullmet = 1
        np_fullmet = 26   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 21 ! omega
        np_fullmet_RH = 25  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp /)
        p_fullmet_Vz_sp(1:21) =  &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  10.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.

      elseif(MR_iwindformat.eq.22)then
        ! GFS 0.25
        call MR_Set_Met_NCEPGeoGrid(193)

        nt_fullmet = 1
        np_fullmet = 31   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 21 ! omega
        np_fullmet_RH = 31  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp,   2.0_sp, &
               1.0_sp /)
        p_fullmet_Vz_sp(1:21) =  &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             850.0_sp, 800.0_sp, 750.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp, 550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp,   2.0_sp, & 
               1.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.23)then
        ! NCEP doE reanalysis
        call MR_Set_Met_NCEPGeoGrid(2)

        nt_fullmet = 1
        np_fullmet = 17
        np_fullmet_Vz = 17 ! omega
        np_fullmet_RH = 17  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
             500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
             150.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp = p_fullmet_sp
        p_fullmet_RH_sp = p_fullmet_sp
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.
      elseif(MR_iwindformat.eq.24)then
        ! NASA MERRA reanalysis
        call MR_Set_Met_NCEPGeoGrid(1024)

        nt_fullmet = 1
        np_fullmet = 42
        np_fullmet_Vz = 42 ! omega
        np_fullmet_RH = 42  ! rhum
        np_fullmet_P0 = 42  ! Precip is 3d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, &
             550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, &
             300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp, &
              70.0_sp,  50.0_sp,  40.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   4.0_sp,   3.0_sp, &
               2.0_sp,   1.0_sp,   0.7_sp,   0.5_sp,   0.4_sp, &
               0.3_sp,   0.1_sp /)
        p_fullmet_Vz_sp = p_fullmet_sp
        p_fullmet_RH_sp = p_fullmet_sp
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.
      elseif(MR_iwindformat.eq.25)then
        ! NCEP-1 1948 reanalysis
        call MR_Set_Met_NCEPGeoGrid(2)
        nt_fullmet = 1460 ! might need to add 4 for a leap year
        np_fullmet = 17   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 12 ! omega
        np_fullmet_RH = 8  ! rhum
        np_fullmet_P0 = 1  ! Precip is 2d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
             500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
             150.0_sp, 100.0_sp,  70.0_sp,  50.0_sp, 30.0_sp, &
              20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
             500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
             150.0_sp, 100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
             500.0_sp, 400.0_sp, 300.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        ! These additional grids are needed since surface variables are on a
        ! different spatial grid.
        do i = 1,192
          x_in_iwf25_sp(i)=(i-1)*1.875_sp
        enddo
        y_in_iwf25_sp(1:94) = (/ &
         88.542_sp,  86.6531_sp,  84.7532_sp,  82.8508_sp,  80.9473_sp,   79.0435_sp,  77.1394_sp, 75.2351_sp, &
        73.3307_sp,  71.4262_sp,  69.5217_sp,  67.6171_sp,  65.7125_sp,   63.8079_sp,  61.9033_sp, 59.9986_sp, &
        58.0939_sp,  56.1893_sp,  54.2846_sp,  52.3799_sp,  50.4752_sp,   48.5705_sp,  46.6658_sp, 44.7611_sp, &
        42.8564_sp,  40.9517_sp,   39.047_sp,  37.1422_sp,  35.2375_sp,   33.3328_sp,  31.4281_sp, 29.5234_sp, &
        27.6186_sp,  25.7139_sp,  23.8092_sp,  21.9044_sp,  19.9997_sp,    18.095_sp,  16.1902_sp, 14.2855_sp, &
        12.3808_sp, 10.47604_sp,  8.57131_sp,  6.66657_sp,  4.76184_sp,    2.8571_sp, 0.952368_sp, &
      -0.952368_sp,  -2.8571_sp, -4.76184_sp, -6.66657_sp, -8.57131_sp, -10.47604_sp, -12.3808_sp, &
       -14.2855_sp, -16.1902_sp,  -18.095_sp, -19.9997_sp, -21.9044_sp,  -23.8092_sp, -25.7139_sp, &
       -27.6186_sp, -29.5234_sp, -31.4281_sp, -33.3328_sp, -35.2375_sp,  -37.1422_sp,  -39.047_sp, &
       -40.9517_sp, -42.8564_sp, -44.7611_sp, -46.6658_sp, -48.5705_sp,  -50.4752_sp, -52.3799_sp, &
       -54.2846_sp, -56.1893_sp, -58.0939_sp, -59.9986_sp, -61.9033_sp,  -63.8079_sp, -65.7125_sp, &
       -67.6171_sp, -69.5217_sp, -71.4262_sp, -73.3307_sp, -75.2351_sp,  -77.1394_sp, -79.0435_sp, &
       -80.9473_sp, -82.8508_sp, -84.7532_sp, -86.6531_sp,  -88.542_sp /)
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.
      elseif(MR_iwindformat.eq.26)then
        ! GFS 0.5
        call MR_Set_Met_NCEPGeoGrid(4)

        nt_fullmet = 1
        np_fullmet = 47   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 37 ! omega
        np_fullmet_RH = 41  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
               2.0_sp,   1.0_sp /)

        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 675.0_sp, 650.0_sp, &
             625.0_sp, 600.0_sp, 575.0_sp, 550.0_sp, 525.0_sp, &
             500.0_sp, 475.0_sp, 450.0_sp, 425.0_sp, 400.0_sp, &
             375.0_sp, 350.0_sp, 325.0_sp, 300.0_sp, 275.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
              10.0_sp /)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.
      elseif(MR_iwindformat.eq.27)then
        ! NOAA reanalysis
        call MR_Set_Met_NCEPGeoGrid(1027)

        nt_fullmet = 1460 ! might need to add 4 for a leap year
        np_fullmet = 24   ! This is for HGT, TMP, UGRD, VGRD
        np_fullmet_Vz = 19 ! omega
        np_fullmet_RH = 19 ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
         ! NOTE: This is the order we ultimately want, but what is in the files
         !       is stored top-down (10->1000).
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 950.0_sp, 900.0_sp, 850.0_sp, 800.0_sp, &
             750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
             500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
             250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp,  70.0_sp, &
              50.0_sp,  30.0_sp,  20.0_sp,  10.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1000.0_sp, 950.0_sp, 900.0_sp, 850.0_sp, 800.0_sp, &
             750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
             500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
             250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp /)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_Vz_sp(1:np_fullmet_Vz)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.28)then
        ! ECMWF Global Gaussian Lat/Lon grid 170
          ! Note: grid is not regular
          !       pressure values are from low to high
        call MR_Set_Met_NCEPGeoGrid(170)

        nt_fullmet = 1
        np_fullmet = 37   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 37 ! omega
        np_fullmet_RH = 37  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
         ! NOTE: This is the order we ultimately want, but what is in the files
         !       is stored top-down (1->1000).
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
             500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
             250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
             125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
              20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
               2.0_sp,   1.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .true.
      elseif(MR_iwindformat.eq.29)then
        ! Japanese 25-year reanalysis
        call MR_Set_Met_NCEPGeoGrid(2)

        nt_fullmet = 1
        np_fullmet = 23   ! This is for air, hgt, uwnd, vwnd
        np_fullmet_Vz = 23 ! omega
        np_fullmet_RH = 23  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
         ! NOTE: This is the order we ultimately want, but what is in the files
         !       is stored top-down (1->1000).
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
             500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
             150.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp,   2.0_sp, &
               1.0_sp,   0.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .true.
        z_inverted = .false.


      elseif(MR_iwindformat.eq.31)then
          ! Catania forecasts
        call MR_Set_Met_NCEPGeoGrid(1031)

        nt_fullmet = 1
        np_fullmet = 13
        np_fullmet_Vz = 13 ! omega
        np_fullmet_RH = 13  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 950.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, &
             600.0_sp, 500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, &
             200.0_sp, 150.0_sp, 100.0_sp /)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = p_fullmet_sp(1:np_fullmet)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .false.

      elseif(MR_iwindformat.eq.32)then
          ! Air Force Weather Agency
        call MR_Set_Met_NCEPGeoGrid(1032)

        nt_fullmet = 1
        np_fullmet = 39
        np_fullmet_Vz = 31 ! omega
        np_fullmet_RH = 39  ! rhum
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1013.0_sp, 1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, &
             900.0_sp,  875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, &
             775.0_sp,  750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp,  550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp,  300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,   70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,    3.0_sp,   2.0_sp,   1.0_sp,   0.5_sp, &
               0.3_sp,    0.2_sp,   0.1_sp,   0.05_sp/)
        p_fullmet_Vz_sp(1:np_fullmet_Vz) = &
          (/1013.0_sp, 1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, &
             900.0_sp,  875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, &
             775.0_sp,  750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, &
             600.0_sp,  550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, &
             350.0_sp,  300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, &
             100.0_sp,   70.0_sp,  50.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp/)
        p_fullmet_RH_sp(1:np_fullmet_RH) = p_fullmet_sp(1:np_fullmet)

        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)) 
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .true.

      elseif(MR_iwindformat.eq.40)then
        ! NASA GEOS Cp
        call MR_Set_Met_NCEPGeoGrid(1040)

        nt_fullmet = 1
        np_fullmet = 37
        np_fullmet_Vz = 37 ! omega
        np_fullmet_RH = 37  ! rhum
        np_fullmet_P0 = 37  ! Precip is 3d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, &
             550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, &
             300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp, &
              70.0_sp,  50.0_sp,  40.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   4.0_sp,   3.0_sp, &
               2.0_sp,   1.0_sp /)
        p_fullmet_Vz_sp = p_fullmet_sp
        p_fullmet_RH_sp = p_fullmet_sp
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .false.

      elseif(MR_iwindformat.eq.41)then
        ! NASA GEOS Np
        call MR_Set_Met_NCEPGeoGrid(1041)

        nt_fullmet = 1
        np_fullmet = 42
        np_fullmet_Vz = 42 ! omega
        np_fullmet_RH = 42  ! rhum
        np_fullmet_P0 = 42  ! Precip is 3d
        allocate(p_fullmet_sp(np_fullmet))
        allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
        allocate(p_fullmet_RH_sp(np_fullmet_RH))
        p_fullmet_sp(1:np_fullmet) = &
          (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
             875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
             750.0_sp, 725.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, &
             550.0_sp, 500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, &
             300.0_sp, 250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp, &
              70.0_sp,  50.0_sp,  40.0_sp,  30.0_sp,  20.0_sp, &
              10.0_sp,   7.0_sp,   5.0_sp,   4.0_sp,   3.0_sp, &
               2.0_sp,   1.0_sp,   0.7_sp,   0.5_sp,   0.4_sp, &
               0.3_sp,   0.1_sp /)
        p_fullmet_Vz_sp = p_fullmet_sp
        p_fullmet_RH_sp = p_fullmet_sp
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))
        p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
        p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa
        x_inverted = .false.
        y_inverted = .false.
        z_inverted = .false.

      else
        ! Not a recognized MR_iwindformat
        ! call reading of custom windfile pressure,grid values
        write(MR_global_error,*)"MR ERROR : windfile format not recognized."
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
        if(MR_runAsForecast)then
          write(MR_global_error,*)"MR ERROR: iwf=27 cannot be used for forecast runs."
          write(MR_global_error,*)"          These are reanalysis files."
          stop 1
        endif
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
        MR_ForecastInterval = 6.0_4
        do iwstep = 1,nt_fullmet
          MR_windfile_stephour(:,iwstep) = (iwstep-1)*MR_ForecastInterval
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
            write(MR_global_info,*)"   For grib2, assume one time step per file."
            write(MR_global_info,*)"   Allocating time arrays for ",MR_iwindfiles,"files"
            write(MR_global_info,*)"                              ",nt_fullmet,"step(s) each"
            allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
          endif

          call grib_open_file(ifile,trim(ADJUSTL(MR_windfiles(iw))),'R')
          call grib_new_from_file(ifile,igrib,iret)
          call grib_get(igrib,'dataDate',dataDate)
          call grib_get(igrib,'dataTime',dataTime)
          call grib_get(igrib,'forecastTime',forecastTime)

          itstart_year  = int(dataDate/10000)
          itstart_month = int((dataDate-10000*itstart_year)/100)
          itstart_day   = mod(dataDate,100)
          itstart_hour  = dataTime
          itstart_min   = 0
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
2100  FORMAT(20x,a11,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)

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
      character(len=130) :: grib2_file

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
        ! Used for keys
      integer(kind=4)  :: discipline
      integer(kind=4)  :: parameterCategory
      integer(kind=4)  :: parameterNumber
      integer(kind=4)  :: level
      integer(kind=4)  :: scaledValueOfFirstFixedSurface
      !character(len=9)   :: sName
      !integer(kind=4)  :: forecastTime
      integer(kind=4),dimension(:),allocatable :: discipline_idx
      integer(kind=4),dimension(:),allocatable :: parameterCategory_idx
      integer(kind=4),dimension(:),allocatable :: parameterNumber_idx
      integer(kind=4),dimension(:),allocatable :: level_idx
      integer(kind=4),dimension(:),allocatable :: forecastTime_idx
      integer(kind=4)  :: disciplineSize
      integer(kind=4)  :: parameterCategorySize
      integer(kind=4)  :: parameterNumberSize
      integer(kind=4)  :: levelSize
      integer(kind=4)  :: forecastTimeSize

      integer(kind=4)  :: iv_discpl
      integer(kind=4)  :: iv_paramC
      integer(kind=4)  :: iv_paramN
      integer(kind=4)  :: iv_typeSf

      real(kind=sp) :: Z_top, T_top
      real(kind=sp),dimension(:,:,:),allocatable :: full_values

      logical :: Use_GRIB2_Index = .false.

      if(.not.Met_var_IsAvailable(ivar))then
        write(MR_global_error,*)"MR ERROR:  Variable not available for this windfile"
        write(MR_global_error,*)"             ivar = ",ivar
        write(MR_global_error,*)"            vname = ",Met_var_names(ivar)
        write(MR_global_error,*)"             iwf  = ",MR_iwindformat
        stop 1
      endif

      !if(ivar.eq.3 .or. ivar.eq.12)then
      !  Use_GRIB2_Index = .false.
      !else
      !  Use_GRIB2_Index = .true.
      !endif

      ! Get the variable discipline, Parameter Catagory, Parameter Number, and
      ! level type for this variable
      iv_discpl = Met_var_GRIB2_DPcPnSt(ivar,1)
      iv_paramC = Met_var_GRIB2_DPcPnSt(ivar,2)
      iv_paramN = Met_var_GRIB2_DPcPnSt(ivar,3)
      iv_typeSf = Met_var_GRIB2_DPcPnSt(ivar,4)

      iw     = MR_MetStep_findex(istep)
      iwstep = MR_MetStep_tindex(istep)

      if(Met_var_names(ivar).eq."")then
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
          np_met_loc = np_fullmet_Vz
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
          np_met_loc = np_fullmet_RH
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
        if(ivar.eq.4)then      ! Vertical_velocity_pressure_isobaric
          np_met_loc = np_fullmet_Vz
          !p_met_loc(1:np_met_loc)  = int(p_fullmet_Vz_sp(1:np_met_loc)/100.0)
          p_met_loc(1:np_met_loc)  = int(p_fullmet_Vz_sp(1:np_met_loc))
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
          np_met_loc = np_fullmet_RH
          !p_met_loc(1:np_met_loc)  = int(p_fullmet_RH_sp(1:np_met_loc)/100.0)
          p_met_loc(1:np_met_loc)  = int(p_fullmet_RH_sp(1:np_met_loc))
        elseif(ivar.eq.40.or.ivar.eq.41.or.ivar.eq.42.or.ivar.eq.43)then ! categorical precip
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        elseif(ivar.eq.44)then ! Precipitation_rate_surface
          np_met_loc = 1
          p_met_loc(1:np_met_loc)  = 0
        else
          np_met_loc = np_fullmet
          !p_met_loc(1:np_met_loc)  = int(p_fullmet_sp(1:np_met_loc)/100.0)
          p_met_loc(1:np_met_loc)  = int(p_fullmet_sp(1:np_met_loc))
        endif
        !write(MR_global_info,*)p_met_loc
        !write(MR_global_info,*)"Allocating full_values ",nx_fullmet,ny_fullmet,np_met_loc
        allocate(full_values(nx_fullmet,ny_fullmet,np_met_loc))
          ! Files are listed directly, not through directories (as in MR_iwindformat=25,27)
        grib2_file = trim(adjustl(MR_MetStep_File(istep)))
        index_file = trim(adjustl(MR_MetStep_File(istep))) // ".index"
      endif
      invar = Met_var_names(ivar)

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

      if(Use_GRIB2_Index)then
        write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                  trim(adjustl(index_file))!,nx_submet,ny_submet,np_met_loc
  
        call grib_index_read(idx,index_file)
        call grib_multi_support_on()
  
          ! get the number of distinct values of all the keys in the index
        call grib_index_get_size(idx,'discipline',disciplineSize)
        call grib_index_get_size(idx,'parameterCategory',parameterCategorySize)
        call grib_index_get_size(idx,'parameterNumber',parameterNumberSize)
        !call grib_index_get_size(idx,'level',levelSize)
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
        !call grib_index_get(idx,'level',level_idx)
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
          !write(MR_global_info,*)discipline_idx(l),parameterCategory_idx(j),&
          !          parameterNumber_idx(k),level_idx(i),forecastTime_idx(t),&
          !          typeOfFirstFixedSurface
          if ( discipline_idx(l)       .eq. iv_discpl .and. &
               parameterCategory_idx(j).eq. iv_paramC .and. &
               parameterNumber_idx(k)  .eq. iv_paramN .and. &
               typeOfFirstFixedSurface .eq. iv_typeSf) then

            !call grib_get(igrib,'shortName',sName)
            !write(MR_global_info,*)"       FOUND :: ",sName
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
               !write(MR_global_info,*)"Checking levels: ",kk,p_met_loc(kk),level_idx(i)
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
                  trim(adjustl(grib2_file))!,nx_submet,ny_submet,np_met_loc
        ifile=5
        call grib_open_file(ifile,grib2_file,'R')
      
        !     turn on support for multi fields messages */
        call grib_multi_support_on()
      
        ! Loop on all the messages in a file.
        call grib_new_from_file(ifile,igrib,iret)
        !write(MR_global_info,*)"  ifile,igrib,iret: ",ifile,igrib,iret
        count1=0
        do while (iret/=GRIB_END_OF_FILE)
          count1=count1+1
          call grib_get(igrib,'discipline',              discipline)
          call grib_get(igrib,'parameterCategory',       parameterCategory)
          call grib_get(igrib,'parameterNumber',         parameterNumber)
          call grib_get(igrib,'typeOfFirstFixedSurface', typeOfFirstFixedSurface)
          !call grib_get(igrib,'level',level)
          call grib_get(igrib,'scaledValueOfFirstFixedSurface',level)

          call grib_get(igrib,'scaledValueOfFirstFixedSurface',scaledValueOfFirstFixedSurface)
          !write(MR_global_info,*)count1,discipline,parameterCategory,parameterNumber,typeOfFirstFixedSurface
          if ( discipline              .eq. iv_discpl .and. &
               parameterCategory       .eq. iv_paramC .and. &
               parameterNumber         .eq. iv_paramN .and. &
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
             ! Now loop through the pressure values for this variable and put
             ! this
             ! slice at the correct level.
             if(ivar.eq.16)level = scaledValueOfFirstFixedSurface
             do kk = 1,np_met_loc
               if(p_met_loc(kk).eq.level)then
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
!          allocate(temp2d_int(nx_submet,ny_submet,1))
!          do i=1,ict        !read subgrid at current time step
!            if(MR_iwindformat.eq.25)then
!              ! No catagorical variables for MR_iwindformat = 25
!            else
!              nSTAT = nf90_get_var(ncid,in_var_id,temp2d_int(ileft(i):iright(i),:,:), &
!                         start = (/iistart(i),jstart,iwstep/),       &
!                         count = (/iicount(i),ny_submet,1/))
!              do j=1,ny_submet
!                itmp = ny_submet-j+1
!                if(y_inverted)then
!                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,itmp,1)
!                else
!                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,j,1)
!                endif
!              enddo
!            endif
!            if(nSTAT.ne.0)then
!               write(MR_global_error,*)'MR ERROR: get_var:Vx ',invar,nf90_strerror(nSTAT)
!               write(MR_global_log  ,*)'MR ERROR: get_var:Vx ',invar,nf90_strerror(nSTAT)
!               stop 1
!             endif
!          enddo
!          deallocate(temp2d_int)
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
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=0.0_sp)
          else
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=0.0_sp)
          endif
        elseif(ivar.eq.5)then
          ! set ground and top-level conditions for temperature
          Z_top = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/real(100.0,kind=sp))
          T_top = MR_Temp_US_StdAtm(Z_top)
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          else
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat), &
                          bc_low_sp=293.0_sp, bc_high_sp=T_top)
          endif
        else
          ! For other variables, use the top and bottom non-fill values
          if(istep.eq.MR_iMetStep_Now)then
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_last,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat))
          else
            call MR_QC_3dvar(nx_submet,ny_submet,np_fullmet,MR_geoH_metP_next,       &
                          np_met_loc,MR_dum3d_metP,fill_value_sp(MR_iwindformat))
          endif
        endif
      endif

      if(ivar.eq.4)then
          ! For pressure vertical velocity, convert from Pa s to m/s by dividing
          ! by pressure gradient
        do k=1,np_met_loc
          do i=1,nx_submet
            do j=1,ny_submet
              if(k.eq.1)then
                ! Use one-sided gradients for bottom
                del_P = p_fullmet_Vz_sp(2)-p_fullmet_Vz_sp(1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,2) - MR_geoH_metP_last(i,j,1)
                else
                  del_H = MR_geoH_metP_next(i,j,2) - MR_geoH_metP_next(i,j,1)
                endif
              elseif(k.eq.np_met_loc)then
                ! Use one-sided gradients for top
                del_P = p_fullmet_Vz_sp(np_met_loc) - &
                         p_fullmet_Vz_sp(np_met_loc-1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,np_met_loc) - &
                           MR_geoH_metP_last(i,j,np_met_loc-1)
                else
                  del_H = MR_geoH_metP_next(i,j,np_met_loc) - &
                           MR_geoH_metP_next(i,j,np_met_loc-1)
                endif
              else
                ! otherwise, two-sided calculation
                del_P = p_fullmet_Vz_sp(k+1)-p_fullmet_Vz_sp(k-1)
                if(istep.eq.MR_iMetStep_Now)then
                  del_H = MR_geoH_metP_last(i,j,k+1) - MR_geoH_metP_last(i,j,k-1)
                else
                  del_H = MR_geoH_metP_next(i,j,k+1) - MR_geoH_metP_next(i,j,k-1)
                endif
              endif
              del_h = del_H * 1000.0_sp ! convert to m
              dpdz  = del_P/del_H
              MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,k) / dpdz
            enddo
          enddo
        enddo
      endif
      MR_dum3d_metP = MR_dum3d_metP * Met_var_conversion_factor(ivar)

      end subroutine MR_Read_MetP_Variable_GRIB

