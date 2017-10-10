![hschwaiger@ramp RAW]$ more ../get_RadioSond.sh 
!#!/bin/bash
!
!StationNum=("70414" "70316" "70326" "70350" "70273")
!StaTionCode=("PASY" "PACD" "PAKN" "PADQ" "PANC")
!StationName=("Shemya Afb" "Cold Bay" "King Salmon" "Kodiak" "Anchorage")
!#yearmonthday=`date -u +%Y%m%d`                        #current year, month &
!day (e.g. 20110119)
!Y=`date -u +%Y`
!M=`date -u +%m`
!DD=`date -u +%d`
!Y=2015
!M=02
!DD=24
!HH=(00 12)
!
!cd /data/windfiles/MetProfiles
!
!for (( si=0;si<=4;si++))
!do
!  #TYPE=TEXT%3AUNMERGED
! echo
!"http://weather.uwyo.edu/cgi-bin/sounding?region=naconf&TYPE=TEXT%3ALIST&YEAR=${Y}&MONTH=${M}&FROM=${DD}${HH[$1]}&TO=${DD}${HH[$1]}&STNM=${StationNum[si]}"
!> tmp.lnk
! wget -i tmp.lnk -O out.html
! cp out.html RAW/${StaTionCode[si]}_${Y}${M}${DD}${HH[$1]}_raw.dat
!
! sed '1,10d' out.html > headless.dat
! sed '/PRE/,$d' headless.dat > ${StaTionCode[si]}_${Y}${M}${DD}${HH[$1]}.dat
! rm tmp.lnk headless.dat out.html
!
!done
!
!PRES:   Atmospheric Pressure    [hPa]
!HGHT:   Geopotential Height     [meter]
!TEMP:   Temperature     [celsius]
!DWPT:   Dewpoint Temperature    [celsius]
!RELH:   Relative Humidity       [%]
!MIXR:   Mixing Ratio    [gram/kilogram]
!DRCT:   Wind Direction  [degrees true]
!SKNT:   Wind Speed      [knot]
!THTA:   Potential Temperature   [kelvin]
!THTE:   Equivalent Potential Temperature        [kelvin]
!THTV:   Virtual Potential Temperature   [kelvin]


!##############################################################################
!
!     Set_MetComp_Grids_1dascii
!
!     Called once from MR_Read_Met_DimVars
!
!     This subroutine needs to determine the full Met grid (i.e. determine the
!     pressure grid and the spatial grid).  It also needs to determine which
!     variables are available.  Generally, the actual region used (e.g. the submet)
!     is determinined through MR_Initialize_Met_Grids after the comp grid is
!     defined.  Values are normally only read on this submet grid.  For 1d
!     ASCII cases, all sonde data are loaded into memory in this subroutine.
!     MR_Initialize_Met_Grids will be used to define the mapping of sonde data
!     to the comp grid.
!
!     This subroutine also allocates and populates the time arrays for the sonde
!     data, something usually done by MR_Read_Met_Times_[netcdf,grib,etc]

!     Note: Every sonde profile will have its own pressure levels, however,
!           the 1d array, p_fullmet_sp, is allocated as the intersection of all
!           the individual pressure profiles so that programs can be data on the
!           metP grid
!     
!     Allocated the dummy arrays for storing met data on met and computational
!     grids:  MR_dum2d_met(nx_submet,ny_submet)
!             MR_dum3d_metP(nx_submet,ny_submet,np_fullmet)
!             MR_dum3d_metH(nx_submet,ny_submet,nz_comp)
!             MR_dum2d_comp(nx_comp,ny_comp)
!             MR_dum3d_compH(nx_comp,ny_comp,nz_comp)
!             MR_windfile_starthour(MR_Snd_nt_fullmet)
!             MR_windfile_stephour(MR_Snd_nt_fullmet,1)
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_ASCII_1d

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=dp), parameter :: PI        = 3.141592653589793_dp
      real(kind=dp), parameter :: DEG2RAD   = 1.7453292519943295e-2_dp

      integer :: ioerr
      integer :: i
      integer :: fid

      integer :: nlev
      integer :: iw_idx
      integer :: iloc, itime
      integer :: p_lidx, p_tidx
      real(kind=sp) :: p_maxtop, p_top

      real(kind=sp),dimension(:),allocatable :: WindVelocity
      real(kind=sp),dimension(:),allocatable :: WindDirection
      real(kind=sp) :: WindTime
      integer :: iw,iws
      real(kind=sp) :: rvalue1,rvalue2,rvalue3,rvalue4,rvalue5
      integer       :: ivalue1,ivalue2,ivalue3,ivalue4,ivalue5

      character(len=80)  :: linebuffer
      logical :: In_hPa = .true.

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------                MR_Read_Met_Times_ASCII_1d                  ----------"
      write(*,*)"--------------------------------------------------------------------------------"

!------------------------------------------------------------------------------
!    MR_iwind.eq.1
!--------------------------------------------------------------------
!      MR_iwindformat.eq.1 (custom 1-d profile)
! L1 string header line
! L2 time(hr) nlev [ncol] [ivar(ncol)]
! L3 x/lon y/lat [LLflag] []
!  data
!    if no optional params, then
!      3 col -> alt(m) windspeed(m/s) winddir(E of N wind FROM)
!      5 col -> alt(m) windspeed(m/s) winddir(E of N wind FROM) pres(hPa) T(degC)
!      ncol  -> variables can be in any order but must include:
!                  1) either pres or alt or both (if one is given, the other is filled from stdAtm)
!                  2) either vx/vy (vector componets) or speed/dir (with 'dir' as wind from)
!                If T is not provided, it is auto-filled from StdAtm
!-------------
!Input wind file for Ashfall tests
!2    41                             #wind time, number of levels
!-122.18   46.20                     #Longitude, latitude of wind sounding (fictitious)
!0       10.000   90.00    101300.    15.0          z, windspeed, direction, p, T
!1000    10.000   90.00     89846.     8.5
!2000    10.000   90.00     79464.     2.0
!--------------------------------------------------------------------
!      MR_iwindformat.eq.2 (radiosonde format as downloaded from "http://weather.uwyo.edu)
!           This format expects 10 lines for header (html lines + column headers)
!           Followed by a block of data in 11 columns
!           Followed by an html tag </PRE><H3>
!           Followed by station information
!
!<HTML>
!<TITLE>University of Wyoming - Radiosonde Data</TITLE>
!<LINK REL="StyleSheet" HREF="/resources/select.css" TYPE="text/css">
!<BODY BGCOLOR="white">
!<H2>70273 PANC Anchorage Observations at 00Z 30 May 2017</H2>
!<PRE>
!-----------------------------------------------------------------------------
!   PRES   HGHT   TEMP   DWPT   RELH   MIXR   DRCT   SKNT   THTA   THTE   THTV
!    hPa     m      C      C      %    g/kg    deg   knot     K      K      K
!-----------------------------------------------------------------------------
! 1008.0     40    7.6    4.9     83   5.41    240      3  280.1  295.2  281.0
! 1003.0     91    5.8    3.8     87   5.03    243      4  278.7  292.7  279.6
! 1000.0    121    5.8    3.7     86   5.01    245      4  278.9  292.9  279.8
!:
!:
!:
!   13.3  29566  -38.5  -66.3      4   0.38     65     16  805.6  809.9  805.8
!   12.2  30175  -37.1  -66.2      3   0.42    110     23  831.2  836.1  831.4
!   11.7  30480  -36.4  -66.2      3   0.45    115     12  844.4  849.5  844.6
!   11.2  30785  -35.7  -66.1      3   0.47     90     20  857.7  863.2  857.9
!   10.7  31072  -35.1  -66.1      3   0.49                870.4  876.3  870.6
!</PRE><H3>Station information and sounding indices</H3><PRE>
!                         Station identifier: PANC
!                             Station number: 70273
!                           Observation time: 170530/0000
!                           Station latitude: 61.16
!                          Station longitude: -150.01
!                          Station elevation: 40.0
!                            Showalter index: 1.68
!                               Lifted index: 2.33
!    LIFT computed using virtual temperature: 2.32
!                                SWEAT index: 228.00
!                                    K index: 21.60
!                         Cross totals index: 28.10
!                      Vertical totals index: 29.80
!                        Totals totals index: 57.90
!      Convective Available Potential Energy: 0.00
!             CAPE using virtual temperature: 0.00
!                      Convective Inhibition: 0.00
!             CINS using virtual temperature: 0.00
!                     Bulk Richardson Number: 0.00
!          Bulk Richardson Number using CAPV: 0.00
!  Temp [K] of the Lifted Condensation Level: 274.14
!Pres [hPa] of the Lifted Condensation Level: 934.80
!     Mean mixed layer potential temperature: 279.49
!              Mean mixed layer mixing ratio: 4.44
!              1000 hPa to 500 hPa thickness: 5289.00
!Precipitable water [mm] for entire sounding: 11.68
!</PRE>

      Met_dim_IsAvailable = .false.
      Met_var_IsAvailable = .false.
      If(MR_iwind.eq.1.and.MR_iwindformat.eq.1)THEN
        ! We are reading just one windfile with the following format
        !  If three values per line:
        !    Elevation(m) , Velocity(m/s) , Direction(degree E of N)
        !  If five values per line:
        !    Elevation(m) , Velocity(m/s) , Direction(degree E of N) , Pressure(hPa) ,
        !     Temperature(C)

        ! x and y fullmet array will just be the coordinates in the order listed
        allocate(x_fullmet_sp(MR_nSnd_Locs))
        allocate(y_fullmet_sp(MR_nSnd_Locs))
        ! The number of time steps/file is fixed to 1
        nt_fullmet = 1
        ! Since there are potentially multiple files for each time (one for each
        ! sonde location), we allocate the following based on MR_Snd_nt_fullmet instead
        ! of MR_iwindfiles
        allocate(MR_windfile_starthour(MR_Snd_nt_fullmet))
        allocate(MR_windfile_stephour(MR_Snd_nt_fullmet,nt_fullmet))
        MR_windfile_starthour = 0.0 ! This is the initialization; will be set below
        MR_windfile_stephour  = 0.0 ! This will be the final value since all files
                                    ! have one step and so no offset.
        MR_windfiles_nt_fullmet(:) = nt_fullmet

        Have_Vz = .false.
        DO itime = 1,MR_Snd_nt_fullmet
          DO iloc = 1,MR_nSnd_Locs
        ! In general, we will want to have multiple locations and multiple times
        ! For now, just set the index of each to 1

        iw_idx = (itime-1)*MR_nSnd_Locs + iloc
        write(*,*)"Opening sonde file ",iw_idx,MR_windfiles(iw_idx)
        fid = 127
        open(unit=fid,file=trim(adjustl(MR_windfiles(1))), status='unknown',err=1971)
        read(fid,*)!skip over first line
        read(fid,'(a80)')linebuffer
        ! Assume we can read at least two values (a real and an interger)
        read(linebuffer,*) rvalue1, ivalue1
        WindTime = rvalue1
        MR_windfile_starthour(iw_idx) = WindTime
        nlev     = ivalue1
        ! Try for three values
        read(linebuffer,*,iostat=ioerr) rvalue1,ivalue1, ivalue2
        IF(ioerr.eq.0)THEN
          ! Success: third value is the number of variables
          !    Note: We need at least 5 variables for P,H,U,V,T
          ! First check if this is the first file read, otherwise do not allocate
          IF(iw_idx.eq.1)THEN
            MR_Snd_nvars = max(5,ivalue2)
            allocate(MR_SndVarsID(MR_Snd_nvars))
          ENDIF
          read(linebuffer,*,iostat=ioerr) rvalue1,ivalue1, ivalue2, MR_SndVarsID(1:MR_Snd_nvars)
        ELSE
          ! If no list of variables is provided, we will still need a list up to 5
          IF(iw_idx.eq.1)THEN
            MR_Snd_nvars = 5
            allocate(MR_SndVarsID(MR_Snd_nvars))
            MR_SndVarsID(1) = 0 ! P
            MR_SndVarsID(2) = 1 ! H
            MR_SndVarsID(3) = 2 ! U
            MR_SndVarsID(4) = 3 ! V
            MR_SndVarsID(5) = 5 ! T
          ENDIF
        ENDIF
          ! This variable will hold all the sonde data
        IF(iw_idx.eq.1)THEN
          allocate(MR_SndVars_metP(MR_nSnd_Locs,MR_Snd_nt_fullmet,MR_Snd_nvars,300))
          allocate(MR_Snd_np_fullmet(MR_nSnd_Locs,MR_Snd_nt_fullmet))
          MR_SndVars_metP   = 0.0
          MR_Snd_np_fullmet = 0
        ENDIF
        MR_Snd_np_fullmet(iloc,itime) = nlev
        read(fid,*) rvalue1,rvalue2
        x_fullmet_sp(iloc) = rvalue1
        y_fullmet_sp(iloc) = rvalue2

        allocate( WindVelocity(nlev));  WindVelocity = 0.0_sp
        allocate(WindDirection(nlev)); WindDirection = 0.0_sp
         !Read elevation (m), speed, direction at each level
         ! Speed is given in m/s (multiply by 0.514444444 to convert
         ! from knots to m/s)
         ! Direction is given as degrees east of north and specifies the
         ! direction FROM which the wind is blowing.
        do i=1,nlev
          rvalue1 = -1.99_4
          rvalue2 = -1.99_4
          rvalue3 = -1.99_4
          rvalue4 = -1.99_4
          rvalue5 = -1.99_4
          read(fid,'(a80)')linebuffer
          read(linebuffer,*,iostat=ioerr) rvalue1, rvalue2, rvalue3
            ! Recall the first five columns are P,H,U,V,T
          MR_SndVars_metP(iloc,itime,2,i) = rvalue1*1.0e-3_sp
          WindVelocity(i)  = rvalue2
          WindDirection(i) = rvalue3
          ! Assume we can at least read three values, try for five
          if (ioerr.eq.0)then
            read(linebuffer,*,iostat=ioerr) rvalue1, rvalue2, rvalue3, &
                                            rvalue4, rvalue5
            if (ioerr.eq.0)then
              Snd_Have_PT = .true.
              ! Five values were successfully read, interpret as:
              ! WindElevation,,WindVelocity,WindDirection,Pressure,Temperature
              IF(iw_idx.eq.1.and.i.eq.1.and.rvalue4.gt.1500.0)THEN
                ! For the first file and first level (ground), check if pressure is greater
                ! than the extected 1013 hPa.  If so, assume pressure is in Pa
                In_hPa = .false.
              ENDIF
              IF(In_hPa)THEN
                MR_SndVars_metP(iloc,itime,1,i) = rvalue4            ! already in Pa
              ELSE
                MR_SndVars_metP(iloc,itime,1,i) = rvalue4 * 100.0_sp ! convert to Pa
              ENDIF
              MR_SndVars_metP(iloc,itime,5,i) = rvalue5 + 273.0_sp
            else
              MR_SndVars_metP(iloc,itime,1,i) = MR_Pres_US_StdAtm(MR_SndVars_metP(iloc,itime,2,i))&
                                                 * 100.0_sp
              MR_SndVars_metP(iloc,itime,5,i) = MR_Temp_US_StdAtm(MR_SndVars_metP(iloc,itime,2,i))
            endif
          endif
          MR_SndVars_metP(iloc,itime,3,i) = real(WindVelocity(i)*sin(pi + DEG2RAD*WindDirection(i)),kind=sp)
          MR_SndVars_metP(iloc,itime,4,i) = real(WindVelocity(i)*cos(pi + DEG2RAD*WindDirection(i)),kind=sp)

          !Met_dim_IsAvailable(1) = .true.  ! Time
          Met_dim_IsAvailable(2) = .true.  ! P
          !Met_dim_IsAvailable(3) = .true.  ! y
          !Met_dim_IsAvailable(4) = .true.  ! x

          Met_var_IsAvailable(1) = .true.  ! GPH
          Met_var_IsAvailable(2) = .true.  ! U
          Met_var_IsAvailable(3) = .true.  ! V
          Met_var_IsAvailable(5) = .true.  ! T

        end do
        close(fid)
        ENDDO
      ENDDO

      ! Here we look for the highest pressure value (lowest altitude) of all the pressure
      ! values to be used for setting the master pressure array (p_fullmet_sp)
      p_maxtop = 0.0
      p_tidx = 0
      p_lidx = 0
      DO itime = 1,MR_Snd_nt_fullmet
        DO iloc = 1,MR_nSnd_Locs
          p_top =  MR_SndVars_metP(iloc,itime,1,MR_Snd_np_fullmet(iloc,itime))
          IF(p_top.gt.p_maxtop)THEN
            p_maxtop = MR_Snd_np_fullmet(iloc,itime)
            p_tidx = itime
            p_lidx = iloc
          ENDIF
        ENDDO
      ENDDO
      np_fullmet = MR_Snd_np_fullmet(p_lidx,p_tidx)
      np_fullmet_Vz = np_fullmet
      np_fullmet_RH = np_fullmet
      allocate(p_fullmet_sp(np_fullmet))
      allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
      allocate(p_fullmet_RH_sp(np_fullmet_RH))
      p_fullmet_sp(1:np_fullmet)=MR_SndVars_metP(p_lidx,p_tidx,1,1:np_fullmet)
      p_fullmet_Vz_sp = p_fullmet_sp
      p_fullmet_RH_sp = p_fullmet_sp
      MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet))

      IF(nt_fullmet.gt.1)THEN
        MR_ForecastInterval = MR_windfile_starthour(MR_nSnd_Locs+1) - &
                               MR_windfile_starthour(1)
      ELSE
        MR_ForecastInterval = 2400.0_dp
      ENDIF

      ELSEIF(MR_iwind.eq.1.and.MR_iwindformat.eq.2)THEN
        ! We are reading radiosonde data from http://weather.uwyo.edu/
        ! Allocate arrays for iGridCode points at iwindfiles times, and 300 levels
        !  
        Have_Vz = .false.


      ENDIF

      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)

      write(*,*)"File, step, Ref, Offset, HoursSince"
      DO iw = 1,MR_iwindfiles
        DO iws = 1,nt_fullmet
          write(*,*)iw,iws,real(MR_windfile_starthour(iw),kind=4),&
                           real(MR_windfile_stephour(iw,iws),kind=4),&
                           real(MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws),kind=4)
        ENDDO
      ENDDO

      write(*,*)"--------------------------------------------------------------------------------"

      return

1971  write(6,*)  'error: cannot find file input wind file.',&
                  '  Program stopped.'
      write(9,*)  'error: cannot find file input wind file.',&
                  '  Program stopped.'
      stop 1

      end subroutine MR_Read_Met_DimVars_ASCII_1d

!##############################################################################
!
!     MR_Set_MetComp_Grids_ASCII_1d
!
!     Called once from MR_Initialize_Met_Grids
!
!     This subroutine sets the mapping between the network of sonde points and
!     the computational grid
!
!     Allocated the dummy arrays for storing met data on met and computational
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

      subroutine MR_Set_MetComp_Grids_ASCII_1d

      use MetReader
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      write(*,*)"--------------------------------------------------------------------------------"
      write(*,*)"----------                          MR_Set_MetComp_Grids_ASCII_1d     ----------"
      write(*,*)"--------------------------------------------------------------------------------"

      If(MR_iwind.eq.1.and.MR_iwindformat.eq.1.or.&
         MR_iwind.eq.1.and.MR_iwindformat.eq.2)THEN
          ! 1d ascii file or Radiosonde case
        nx_submet     = 1
        ny_submet     = 1
        allocate(x_submet_sp(nx_submet))
        allocate(y_submet_sp(ny_submet))
        allocate(MR_dum3d_metP(nx_submet,ny_submet,np_fullmet))
        allocate(MR_dum3d_metH(nx_submet,ny_submet,nz_comp))
        allocate(MR_dum2d_comp_int(nx_comp,ny_comp))
        allocate(MR_dum2d_comp(nx_comp,ny_comp))
        allocate(MR_dum3d_compH(nx_comp,ny_comp,nz_comp))
        allocate(MR_geoH_metP_last(nx_submet,ny_submet,np_fullmet)) ! This will hold elevation
        allocate(MR_geoH_metP_next(nx_submet,ny_submet,np_fullmet)) ! Copy of geoH_metP_last

        allocate(MR_Snd2Comp_tri_map_wgt(nx_comp,ny_comp,3))
        allocate(MR_Snd2Comp_tri_map_idx(nx_comp,ny_comp,3))
        IF(MR_nSnd_Locs.eq.1)THEN
          ! Easy case: all weights are on the one sonde point
          MR_Snd2Comp_tri_map_wgt = 0.0_sp
          MR_Snd2Comp_tri_map_wgt(1:nx_comp,1:ny_comp,1)   = 1.0_sp
          MR_Snd2Comp_tri_map_idx(1:nx_comp,1:nx_comp,1:3) = 1
        ELSEIf(MR_nSnd_Locs.eq.2)THEN
          ! Also somewhat easy, but not yet implemented
          ! Find each comp point wrt the the two sonde locations.
          !  Use bilinear for points between and sonde points for points beyond
          write(*,*)"MR ERROR:  Multiple sonde locations not yet implemented"
          stop 1
        ELSE
          ! Difficult case:  Here we need to triangulate each comp point with
          ! nearby sonde locations.
          write(*,*)"MR ERROR:  Multiple sonde locations not yet implemented"
          stop 1
        ENDIF

      ELSE
        write(*,*)"Unknown ASCII wind file format."
        write(*,*)"  MR_iwind = ",MR_iwind
        write(*,*)"MR_iwindformat = ",MR_iwindformat
        stop 1
      ENDIF

      write(*,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Set_MetComp_Grids_ASCII_1d


!##############################################################################
!
!     MR_Set_Met_Times_ASCII
!
!     Called once from 
!
!     This subroutine 
!
!     Allocated the dummy arrays for storing met data on met and computational
!
!##############################################################################
!
!      subroutine MR_Set_Met_Times_ASCII
!
!      use MetReader
!      use projection
!
!      implicit none
!
!      integer, parameter :: sp        = 4 ! single precision
!      integer, parameter :: dp        = 8 ! double precision
!
!      If(MR_iwind.eq.1.and.MR_iwindformat.eq.1)THEN
!          ! 1d ascii file
!        allocate(MR_MetStep_File(2))
!        allocate(MR_MetStep_tindex(2))
!        allocate(MR_MetStep_Hour_since_baseyear(2))
!        allocate(MR_MetStep_Interval(1))
!        MR_MetSteps_Total    = 2
!        MR_ForecastInterval = 10000.0_dp       ! give it a long forecast interval
!        MR_MetStep_Interval(1) = MR_ForecastInterval
!        MR_MetStep_File(1)   = MR_windfiles(1) ! Currently, only one sounding
!        MR_MetStep_tindex(1) = 1            ! Only one time with sounding file
!        MR_MetStep_Hour_since_baseyear(1) = MR_Comp_StartHour - 1.0_sp
!
!        MR_MetStep_File(2)   = MR_windfiles(1) ! Currently, only one sounding
!        MR_MetStep_tindex(2) = 2            ! Only one time with sounding file
!        MR_MetStep_Hour_since_baseyear(2) = MR_Comp_StartHour + MR_Comp_Time_in_hours + 1.0_sp
!      ELSEIF(MR_iwind.eq.1.and.MR_iwindformat.eq.2)THEN
!        ! Radiosonde case
!      ELSEIF(MR_iwind.eq.2.and.MR_iwindformat.eq.1)THEN
!        ! 3d ascii windfiles
!      ELSE
!        write(*,*)"Unknown ASCII wind file format."
!        write(*,*)"  MR_iwind = ",MR_iwind
!        write(*,*)"MR_iwindformat = ",MR_iwindformat
!        stop 1
!      ENDIF
!
!      end subroutine MR_Set_Met_Times_ASCII


!##############################################################################
!
!     MR_Read_MetP_Variable_ASCII
!
!     Called from MR_Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_ASCII_1d(ivar,istep)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: ivar,istep

      integer :: icol, i, j
      integer :: itime
      integer :: iloc1, iloc2, iloc3

      itime = MR_MetStep_findex(istep)
            ! these variables are set in MR_Read_Met_DimVars_ASCII_1d
      !  P,H,U,V,T
      IF(Met_var_IsAvailable(ivar))THEN
        DO i=1,MR_Snd_nvars
          IF(ivar.eq.MR_SndVarsID(i))THEN
            exit
          ENDIF
        ENDDO
        icol  = i
        ! Now loop over all the submet points and build the value needed by multiplying by the
        ! weights determined in MR_Set_MetComp_Grids_ASCII_1d
        ! Note: This is a placeholder for multi-sonde triangulation, but only 1 sonde location is
        !       currently supported.
        DO i = 1,nx_submet
          DO j = 1,ny_submet
            iloc1 = MR_Snd2Comp_tri_map_idx(i,j,1)  ! Right now, these should all be 1, with weights of (1,0,0)
            iloc2 = MR_Snd2Comp_tri_map_idx(i,j,2)
            iloc3 = MR_Snd2Comp_tri_map_idx(i,j,3)
            MR_dum3d_metP(i,j,1:np_fullmet) = &
              MR_SndVars_metP(iloc1,itime,icol,1:np_fullmet) * MR_Snd2Comp_tri_map_wgt(i,j,1) + &
              MR_SndVars_metP(iloc2,itime,icol,1:np_fullmet) * MR_Snd2Comp_tri_map_wgt(i,j,2) + &
              MR_SndVars_metP(iloc3,itime,icol,1:np_fullmet) * MR_Snd2Comp_tri_map_wgt(i,j,3)
          ENDDO
        ENDDO
      ELSE
        ! W is typically not provided
        MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_fullmet) = 0.0_sp
      ENDIF

      return

      end subroutine MR_Read_MetP_Variable_ASCII_1d

