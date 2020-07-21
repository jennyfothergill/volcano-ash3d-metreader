!##############################################################################
!##############################################################################
      program ncMetTraj

      use MetReader

      implicit none

      integer             :: nargs
      integer             :: status
      character (len=100) :: arg

      real(kind=8)        :: inlon,inlat
      integer             :: inyear,inmonth,inday
      real(kind=8)        :: inhour
      integer             :: TrajFlag  ! .ge. 0 for forward, .lt.0 for backward

      integer             :: FC_freq = 12
      integer             :: nxmax,nymax,nzmax,nsize
      real(kind=8)        :: dx,dy
      real(kind=4),dimension(:)    ,allocatable :: lon_grid
      real(kind=4),dimension(:)    ,allocatable :: lat_grid
      real(kind=4),dimension(:)    ,allocatable :: z_cc
      real(kind=4),dimension(:)    ,allocatable :: gsdiam
      logical             :: IsPeriodic
      integer             :: i

      integer             :: ntraj             = 0
      real(kind=8)        :: Simtime_in_hours
      real(kind=8)        :: starty

      real(kind=8), dimension(:),allocatable :: OutputLevels

      integer      :: iprojflag
      real(kind=8) :: lambda0,phi0,phi1,phi2,k0,radius_earth
      logical      :: IsLatLon
      real(kind=4) :: tmp_4

      integer :: BaseYear = 1900
      logical :: useLeap  = .true.

      INTERFACE
        subroutine GetWindFile_FC(inyear,inmonth,inday,inhour,&
                                FC_freq,Simtime_in_hours,TrajFlag)
          integer,parameter   :: dp        = 8 ! double precision
          integer             :: inyear,inmonth,inday
          real(kind=dp)       :: inhour
          integer             :: FC_freq
          real(kind=dp)       :: Simtime_in_hours
          integer             :: TrajFlag
        end subroutine
        subroutine Integrate_Trajectory(inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,OutputLevels)

          integer,parameter   :: dp        = 8 ! double precision
          real(kind=dp)       :: inlon
          real(kind=dp)       :: inlat
          integer             :: inyear,inmonth,inday
          real(kind=dp)       :: inhour
          real(kind=dp)       :: Simtime_in_hours
          integer             :: TrajFlag
          integer             :: ntraj
          real(kind=dp), dimension(ntraj) :: OutputLevels
        end subroutine
      END INTERFACE

      ! Make user MetReader is using the same calendar
      MR_BaseYear = BaseYear
      MR_useLeap  = useLeap

!     TEST READ COMMAND LINE ARGUMENTS
      nargs = command_argument_count()
      if (nargs.lt.6) then
        write(MR_global_info,*)"Enter lon,lat,YYYY MM DD HH (FC_hours nlev lev1 lev2 ...)"
        stop 1
      else
        call get_command_argument(1, arg, status)
        read(arg,*)inlon
        if(inlon.lt.-360.0)then
          write(MR_global_info,*)"ERROR: Longitude must be gt -360"
          stop 1
        endif
        if(inlon.lt.0.0_8.or.inlon.gt.360.0_8)inlon=mod(inlon+360.0_8,360.0_8)
        call get_command_argument(2, arg, status)
        read(arg,*)inlat
        call get_command_argument(3, arg, status)
        read(arg,*)inyear
        call get_command_argument(4, arg, status)
        read(arg,*)inmonth
        call get_command_argument(5, arg, status)
        read(arg,*)inday
        call get_command_argument(6, arg, status)
        read(arg,*)inhour

        if(nargs.gt.6)then
          call get_command_argument(7, arg, status)
          read(arg,*)Simtime_in_hours
          write(MR_global_info,*)"Calculating trajectories for ",&
                     Simtime_in_hours," hours."
          if(nargs.gt.7)then
            call get_command_argument(8, arg, status)
            read(arg,*)ntraj
            if(ntraj.le.0)then
              write(MR_global_info,*)"Error reading ntraj."
              write(MR_global_info,*)"ntraj = ",ntraj
              write(MR_global_info,*)"ntraj must be positive."
              stop 1
            elseif(ntraj.gt.9)then
              write(MR_global_info,*)"ERROR: ntraj is currently limited to 9"
              stop 1
            endif
            allocate(OutputLevels(ntraj))
            if(nargs-8.lt.ntraj)then
              write(MR_global_info,*)"ERROR:  There are not enough arguments for ",&
                        ntraj," levels"
            elseif(nargs-8.gt.ntraj)then
              write(MR_global_info,*)"WARNING:  There are more trajectory levels given than needed"
              write(MR_global_info,*)"  Expected ntraj = ",ntraj
              write(MR_global_info,*)"  Extra sommand line arguments = ",nargs-8
            endif
            do i=1,ntraj
              call get_command_argument(8+i, arg, status)
              read(arg,*)OutputLevels(i)
              if(OutputLevels(i).le.0.0.or.OutputLevels(i).gt.30.0)then
                write(MR_global_info,*)"ERROR: trajectory levels must be in range 0-30 km"
                write(MR_global_info,*)"Failing on trajectory ",i,OutputLevels(i)
                stop 1
              endif
            enddo
          endif
        else
          Simtime_in_hours = 24.0_8
        endif
        if(ntraj.eq.0)then
          ntraj = 6
          allocate(OutputLevels(ntraj))
          OutputLevels(1) =  1.524_4 !  5000 ft
          OutputLevels(2) =  3.048_4 ! 10000 ft
          OutputLevels(3) =  6.096_4 ! 20000 ft
          OutputLevels(4) =  9.144_4 ! 30000 ft
          OutputLevels(5) = 12.192_4 ! 40000 ft
          OutputLevels(6) = 15.240_4 ! 50000 ft
        endif

        write(MR_global_info,*)"Calculating ",ntraj," trajectories:"
        do i=1,ntraj
          tmp_4 = real(OutputLevels(i),kind=4)
          write(MR_global_info,*)i," at ",tmp_4,"km (",tmp_4*3280.8_4," ft)."
        enddo

      endif

      TrajFlag = 0
#ifdef FORWARD
      TrajFlag = 1
#endif
#ifdef BACKWARD
      TrajFlag = -1
#endif

      if(TrajFlag.eq.0)then
        write(MR_global_info,*)"ERROR: Forward/Backward not specified."
        write(MR_global_info,*)"       Recompile with preprocessor flags"
        write(MR_global_info,*)"         FORWARD for forward trajectories"
        write(MR_global_info,*)"         BACKWARD for backward trajectories"
        stop 1
      elseif(TrajFlag.lt.0)then
        write(MR_global_info,*)"Calculating Backward trajectories from "
      elseif(TrajFlag.gt.0)then
        write(MR_global_info,*)"Calculating Forward trajectories from "
      endif
      write(MR_global_info,*)"Longitude = ",inlon
      write(MR_global_info,*)"Latitude  = ",inlat
      write(MR_global_info,*)"Year      = ",inyear
      write(MR_global_info,*)"Month     = ",inmonth
      write(MR_global_info,*)"Day       = ",inday
      write(MR_global_info,*)"Hour      = ",inhour

      nzmax = ntraj
      ! Define grid padding based on the integration time
      if(Simtime_in_hours.le.8.0_8)then
        ! +-15 in lon ; +-10 in lat
        nxmax =  61; dx = 0.5_4
        nymax =  41; dy = 0.5_4
      elseif(Simtime_in_hours.le.16.0_8)then
        ! +-25 in lon ; +-15 in lat
        nxmax = 101; dx = 0.5_4
        nymax =  61; dy = 0.5_4
      elseif(Simtime_in_hours.le.24.0_8)then
        ! +-35 in lon ; +-20 in lat
        nxmax = 141; dx = 0.5_4
        nymax =  81; dy = 0.5_4
      else
        write(MR_global_info,*)"ERROR:  Simulation times greater that 24 hours are"
        write(MR_global_info,*)"        currently disallowed."
        stop 1
      endif
      nsize = 1 ! This is also not really used in this utility
      allocate(lon_grid(nxmax))
      allocate(lat_grid(nymax))
      allocate(z_cc(nzmax))
      allocate(gsdiam(nsize))  ; gsdiam(1) = 1.0_4
      do i=1,nxmax
        lon_grid(i) = real(inlon - 0.5*(nxmax-1) * dx + (i-1) * dx,kind=4)
      enddo
      ! We can specify the longitude grid with no problems, but the latitude
      ! grid might be padded up across the pole.  We need to set the cap at
      ! an extreme point (89 degrees N or S)
      if((inlat + (nymax-1)/2 * dy).gt.89.0)then
        ! Start from 89.0 N and count down nymax
        starty = 89.0 - (nymax-1) * dy
        do i=1,nymax
          lat_grid(i) = real(starty + (i-1) * dy,kind=4)
        enddo
      elseif((inlat - (nymax-1)/2 * dy).lt.-89.0)then
        ! Start from 89.0 N and count down nymax
        starty = -89.0
        do i=1,nymax
          lat_grid(i) = real(starty + (i-1) * dy,kind=4)
        enddo
      else
        ! lat grid doesn't involve poles; center grid over inlat
        do i=1,nymax
          lat_grid(i) = real(inlat - 0.5*(nymax-1) * dy + (i-1) * dy,kind=4)
        enddo
      endif
      IsPeriodic = .false.
      do i = 1,ntraj
        z_cc(i) = real(OutputLevels(i),4)
      enddo

      ! Need to modify this to start correctly for backtrajectories
      call GetWindFile_FC(inyear,inmonth,inday,inhour,FC_freq, &
                          Simtime_in_hours,TrajFlag)

      IsLatLon  = .true.
      iprojflag = 1
      lambda0   = -105.0
      phi0      = 90.0
      phi1      = 90.0
      phi2      = 90.0
      k0        = 0.933
      radius_earth = 6371.229
      call MR_Set_CompProjection(IsLatLon,iprojflag,lambda0,phi0,phi1,phi2,&
                                 k0,radius_earth)
      write(MR_global_info,*)"Setting up wind grids"
      call MR_Initialize_Met_Grids(nxmax,nymax,nzmax,&
                              lon_grid(1:nxmax), &
                              lat_grid(1:nymax), &
                              z_cc(1:nzmax)    , &
                              IsPeriodic)

      call Integrate_Trajectory(inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,OutputLevels)

      write(MR_global_info,*)"Program ended normally."

      end program ncMetTraj

!##############################################################################
!
!     Subroutines
!
!##############################################################################

!##############################################################################
!##############################################################################

      subroutine GetWindFile_FC(inyear,inmonth,inday,inhour,&
                                FC_freq,Simtime_in_hours,TrajFlag)

      use MetReader

      implicit none

      integer            :: inyear,inmonth,inday
      real(kind=8)       :: inhour
      integer            :: FC_freq
      real(kind=8)       :: Simtime_in_hours
      integer            :: TrajFlag

      real(kind=8)       :: HS_hours_since_baseyear
      character(len=13)  :: HS_yyyymmddhhmm_since   ! function that calculates date
                                                      !  string given
                                                      !  hours since
                                                      !  MR_BaseYear
      integer            :: HS_YearOfEvent
      integer            :: HS_MonthOfEvent
      integer            :: HS_DayOfEvent
      real(kind=8)       :: HS_HourOfDay


      character(len=8)   :: date
      character(LEN=10)  :: time2                 ! time argument used to get current
                                                 !  date and time.
      character(len=5)   :: zone  !variables used by the date_and_time subroutine
      integer            :: values(8)                !return values from date_and_time
      integer            :: timezone  ! timezone of grid relative to UTC

      real(kind=8)      :: StartHour
      real(kind=8)      :: RunStartHour    ! Current UTC time, in hours since MR_BaseYear
      character(len=17) :: RunStartHour_ch
      real(kind=8)      :: Probe_StartHour
      real(kind=8)      :: FC_Package_StartHour
      real(kind=8)      :: Met_needed_StartHour

      real(kind=8) :: GFS_Avail_Delay = 10.0
      integer      :: GFS_Archive_Days = 14

      integer      :: RunStartYear
      integer      :: RunStartMonth
      integer      :: RunStartDay
      integer      :: RunStartHr

      integer      :: FC_Package_year
      integer      :: FC_Package_month
      integer      :: FC_Package_day
      integer      :: FC_Package_hour
      integer      :: FC_hour_int

      integer      :: iw,iwf,igrid,iwfiles
      character(len=100) :: WINDROOT = '/data/WindFiles'
      character(len=47) :: string1,string2

      integer :: i
      integer :: FC_year,FC_mon,FC_day
      real(kind=8) :: FC_hour,FC_Package_hour_dp,FC_intvl

       ! Get the UTC time that the program is called
       !   This will be used to determine if gfs or NCEP winds are to be used
      call date_and_time(date,time2,zone,values)
      read(zone,'(i3)') timezone
        ! FIND TIME IN UTC
      StartHour = real(values(5)-timezone,kind=8) + &
                  real(values(6)/60.0,kind=8)
        ! find time in hours since BaseYear
      RunStartHour = HS_hours_since_baseyear(values(1),values(2),values(3),&
                                             StartHour,MR_BaseYear,MR_useLeap)
        ! get character string
      RunStartHour_ch = HS_yyyymmddhhmm_since(RunStartHour,MR_BaseYear,MR_useLeap)
      read(RunStartHour_ch,'(i4)') RunStartYear
      read(RunStartHour_ch,'(4x,i2)') RunStartMonth
      read(RunStartHour_ch,'(6x,i2)') RunStartDay
      read(RunStartHour_ch,'(8x,i2)') RunStartHr

      ! Find the start time given on the command line
      Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      MR_Comp_StartHour = Probe_StartHour
      MR_Comp_Time_in_hours = Simtime_in_hours

      ! Calculate the earliest Met data needed
      if(TrajFlag.gt.0)then
        Met_needed_StartHour = Probe_StartHour
      else
        Met_needed_StartHour = Probe_StartHour - Simtime_in_hours
      endif
      ! Now find the Forecast block that starts immediately before the
      ! called time
      FC_Package_year    = HS_YearOfEvent(Met_needed_StartHour,MR_BaseYear,MR_useLeap)
      FC_Package_month   = HS_MonthOfEvent(Met_needed_StartHour,MR_BaseYear,MR_useLeap)
      FC_Package_day     = HS_DayOfEvent(Met_needed_StartHour,MR_BaseYear,MR_useLeap)
      FC_Package_hour_dp = HS_HourOfDay(Met_needed_StartHour,MR_BaseYear,MR_useLeap)
      FC_Package_hour  = floor(FC_Package_hour_dp/FC_freq) * FC_freq
      FC_Package_StartHour = HS_hours_since_baseyear(FC_Package_year,&
                                FC_Package_month,FC_Package_day,real(FC_Package_hour,kind=8),&
                                MR_BaseYear,MR_useLeap)

      write(MR_global_info,*)"Probe_StartHour      = ",Probe_StartHour
      write(MR_global_info,*)"Met_needed_StartHour = ",Met_needed_StartHour
      write(MR_global_info,*)"FC_Package_year    = ",FC_Package_year
      write(MR_global_info,*)"FC_Package_month   = ",FC_Package_month
      write(MR_global_info,*)"FC_Package_day     = ",FC_Package_day
      write(MR_global_info,*)"FC_Package_hour_dp = ",FC_Package_hour_dp
      write(MR_global_info,*)"FC_Package_hour    = ",FC_Package_hour
      write(MR_global_info,*)"FC_Package_StartHour = ",FC_Package_StartHour

      if(RunStartHour-FC_Package_StartHour.lt.GFS_Avail_Delay)then
        ! The closest forecast package to the probe time is too close to
        ! the current (real) time.  The GFS files are probably not yet
        ! available.  Decrement the forecast package.
        FC_Package_StartHour = FC_Package_StartHour - real(FC_freq,kind=8)
      endif
      FC_year = HS_YearOfEvent(FC_Package_StartHour,MR_BaseYear,MR_useLeap)
      FC_mon  = HS_MonthOfEvent(FC_Package_StartHour,MR_BaseYear,MR_useLeap)
      FC_day  = HS_DayOfEvent(FC_Package_StartHour,MR_BaseYear,MR_useLeap)
      FC_hour = HS_HourOfDay(FC_Package_StartHour,MR_BaseYear,MR_useLeap)
      FC_Package_hour = floor(FC_hour/FC_freq) * FC_freq

      if(RunStartHour-Probe_StartHour.gt.(GFS_Archive_Days*24))then
        ! Run is older than 2 weeks, use NCEP winds
        write(MR_global_info,*)"Start time is older than the hardwired threshold of ",&
                  GFS_Archive_Days," days"
        write(MR_global_info,*)"Using the NCEP reanalysis"
        iw  = 5
        iwf = 25
        igrid   = 0
        iwfiles = 1

         !  Gather all the time slices needed in chronological order
        call MR_Allocate_FullMetFileList(iw,iwf,igrid,2,iwfiles) !,FC_year
                  !Met_needed_StartHour,Simtime_in_hours)
        do i=1,MR_iwindfiles
          write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)), &
                               '/NCEP'
        enddo

      elseif(RunStartHour-Probe_StartHour.lt.-90)then
        ! Run is too far in the future
        write(MR_global_info,*)"ERROR: run is too far in future"
        stop 1
      else
        ! Run is newer than 2 weeks, use GFS winds
        iw      = 4
        iwf     = 20
        igrid   = 0
        iwfiles = 34
        FC_intvl = 3.0

          !  Gather all the time slices needed in chronological order
        call MR_Allocate_FullMetFileList(iw,iwf,igrid,2,iwfiles) !,FC_year,   &
                  !Met_needed_StartHour,Simtime_in_hours)
        do i=1,MR_iwindfiles
          FC_hour_int = nint((i-1)*FC_intvl)
          write(string1,'(a9,I4.4,I2.2,I2.2,I2.2,a1)')'/gfs/gfs.', &
                        FC_year,FC_mon,FC_day,FC_Package_hour,'/'
          write(string2,'(I4.4,I2.2,I2.2,I2.2,a2,I3.3,a3)')&
                        FC_year,FC_mon,FC_day,FC_Package_hour, &
                        '.f',FC_hour_int,'.nc'

          write(MR_windfiles(i),*)trim(ADJUSTL(WINDROOT)), &
                               trim(ADJUSTL(string1)), &
                               trim(ADJUSTL(string2))
        enddo

      endif
        ! Check for existance and compatibility with simulation time requirements
      call MR_Read_Met_DimVars(FC_year)

      call MR_Set_Met_Times(Met_needed_StartHour, Simtime_in_hours)


      write(MR_global_info,*)"Traj time: ",inyear,inmonth,inday,inhour
      write(MR_global_info,*)"Now      : ",RunStartYear,RunStartMonth,RunStartDay,RunStartHr
      write(MR_global_info,*)"FC  time : ",FC_year,FC_mon,FC_day,FC_Package_hour

      end subroutine GetWindFile_FC

!##############################################################################
!##############################################################################

      subroutine Integrate_Trajectory(inlon,inlat,inyear,inmonth,inday,inhour,&
                                Simtime_in_hours,TrajFlag,ntraj,OutputLevels)

      use MetReader

      implicit none

      real(kind=8)       :: inlon
      real(kind=8)       :: inlat
      integer            :: inyear,inmonth,inday
      real(kind=8)       :: inhour
      real(kind=8)       :: Simtime_in_hours
      integer            :: TrajFlag
      integer            :: ntraj
      real(kind=8), dimension(ntraj) :: OutputLevels

      real(kind=8), parameter :: PI        = 3.141592653589793
      real(kind=8), parameter :: DEG2RAD   = 1.7453292519943295e-2
      real(kind=8), parameter :: KM_2_M       = 1.0e3
      real(kind=8), parameter :: RAD_EARTH   = 6371.229 ! Radius of Earth in km


      real(kind=8)       :: HS_hours_since_baseyear
      real(kind=8)       :: Probe_StartHour
      integer            :: ivar
      integer            :: i,k,kk
      integer            :: ix,iy

      real(kind=4),dimension(:,:,:),allocatable :: Vx_meso_last_step_MetP_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vx_meso_next_step_MetP_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vy_meso_last_step_MetP_sp
      real(kind=4),dimension(:,:,:),allocatable :: Vy_meso_next_step_MetP_sp
      real(kind=4),dimension(:),allocatable :: GPHprof1,GPHprof2

      real(kind=8) :: tfrac,tc
      real(kind=8) :: xfrac,xc,yfrac,yc
      real(kind=4) :: a1,a2,a3,a4
      real(kind=8),dimension(:),allocatable :: Hinit
      integer,       dimension(ntraj) :: OutLev_indx
      real(kind=8), dimension(ntraj) :: zfrac
      real(kind=8) :: zc
      real(kind=8), dimension(ntraj) :: x1,y1

      real(kind=8),dimension(:,:,:,:),allocatable :: Vx_full
      real(kind=8),dimension(:,:,:,:),allocatable :: Vy_full
      real(kind=8),dimension(:)      ,allocatable :: Step_Time_since1900
      real(kind=8),dimension(:,:,:)  ,allocatable :: dvxdt
      real(kind=8),dimension(:,:,:)  ,allocatable :: dvydt

      integer      :: istep,stepindx
      integer      :: ti,iit,it
      real(kind=8) :: vx1,vx2,vx3,vx4
      real(kind=8) :: vy1,vy2,vy3,vy4
      real(kind=8),dimension(2)  :: vel_1
      real(kind=8) :: dt
      real(kind=8) :: mstodeghr
      real(kind=8) :: t1
      real(kind=8) :: x_fin,y_fin
      real(kind=8) :: xstep,ystep
      real(kind=8) :: lonmin,lonmax,latmin,latmax

      allocate(Hinit(np_fullmet))
      allocate(GPHprof1(np_fullmet))
      allocate(GPHprof2(np_fullmet))
      ! These store the full subgrid from the Met files
      allocate(Vx_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Vx_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Vy_meso_last_step_MetP_sp(nx_submet,ny_submet,np_fullmet))
      allocate(Vy_meso_next_step_MetP_sp(nx_submet,ny_submet,np_fullmet))

      ! These store just the layers relevant, but for all time
      allocate(Vx_full(nx_submet,ny_submet,ntraj,MR_MetSteps_Total))
      allocate(Vy_full(nx_submet,ny_submet,ntraj,MR_MetSteps_Total))
      allocate(Step_Time_since1900(MR_MetSteps_Total))
      ! These are needed for each integration point
      allocate(dvxdt(nx_submet,ny_submet,ntraj))
      allocate(dvydt(nx_submet,ny_submet,ntraj))

      lonmin = 360.0
      lonmax =   0.0
      latmin =  90.0
      latmax = -90.0

       ! Load the full sub-grid for all times
        ! First load the Met grids for Geopotential
      if(TrajFlag.gt.0)then
        ! Foreward trajectory
        MR_iMetStep_Now = 1
      else
        ! Backward trajectory
        MR_iMetStep_Now = MR_MetSteps_Total-1 
      endif
      call MR_Read_HGT_arrays(MR_iMetStep_Now)

      ! Get the fractional time between forecast steps
      Probe_StartHour = HS_hours_since_baseyear(inyear,inmonth,inday,inhour,&
                                                MR_BaseYear,MR_useLeap)
      tfrac = (Probe_StartHour-MR_MetStep_Hour_since_baseyear(MR_iMetStep_Now))/ &
               !MR_ForecastInterval
               MR_MetStep_Interval(MR_iMetStep_Now)
      tc    = 1.0_8-tfrac
      ! Get the fractional position in cell and corner weights
        ! First get indicies of LL of strat point
      do i = 1,nx_submet-1
        if(inlon.ge.x_submet_sp(i).and.inlon.lt.x_submet_sp(i+1))ix=i
      enddo
      do i = 1,ny_submet-1
        if(inlat.ge.y_submet_sp(i).and.inlat.lt.y_submet_sp(i+1))iy=i
      enddo
      xfrac=(inlon - x_submet_sp(ix))/abs(dx_met_const)
      yfrac=(inlat - y_submet_sp(iy))/abs(dy_met_const)
      xc = 1.0_4-xfrac
      yc = 1.0_4-yfrac
      a1 = real(xc*yc,kind=4)
      a2 = real(xfrac*yc,kind=4)
      a3 = real(xfrac*yfrac,kind=4)
      a4 = real(yfrac*xc,kind=4)
      GPHprof1  = a1*MR_geoH_metP_last(ix  ,iy  ,:) + &
                  a2*MR_geoH_metP_last(ix+1,iy  ,:) + &
                  a3*MR_geoH_metP_last(ix+1,iy+1,:) + &
                  a4*MR_geoH_metP_last(ix  ,iy+1,:)
      GPHprof2  = a1*MR_geoH_metP_next(ix  ,iy  ,:) + &
                  a2*MR_geoH_metP_next(ix+1,iy  ,:) + &
                  a3*MR_geoH_metP_next(ix+1,iy+1,:) + &
                  a4*MR_geoH_metP_next(ix  ,iy+1,:)
      do k = 1,np_fullmet
        Hinit(k) = (GPHprof1(k)*tc + GPHprof2(k)*tfrac)
      enddo

      ! The mapping should vary in space and time, but we'll just use
      ! the mapping at the start point/time for now
      !do i = 1,nx_met_per
      !  do j = 1,ny_met
          do kk = 1,ntraj
            do k = 1,np_fullmet-1
              if(kk.eq.1.and.k.eq.1.and.Hinit(k).ge.OutputLevels(kk))then
                OutLev_indx(kk) = k
                zfrac(kk)    = 1.0
                write(MR_global_info,*)kk,OutputLevels(kk),OutLev_indx(kk),zfrac(kk)
                exit
              elseif(Hinit(k).lt.OutputLevels(kk).and.Hinit(k+1).gt.OutputLevels(kk))then
                OutLev_indx(kk) = k
                zfrac(kk)    = (OutputLevels(kk)-Hinit(k))/(Hinit(k+1)-Hinit(k))
                exit
              endif
            enddo
          enddo
      !  enddo
      !enddo

      ! We need to reset the iMetStep_Now since we are now loading all
      ! the wind values and back trajectories will have iMetStep_Now at
      ! MetSteps_Total-1 for probe point calculations
      MR_iMetStep_Now = 1

      ivar = 2 ! Vx
      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
      Vx_meso_last_step_MetP_sp = MR_dum3d_MetP
      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
      Vx_meso_next_step_MetP_sp = MR_dum3d_MetP
      ivar = 3 ! Vy
      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now)
      Vy_meso_last_step_MetP_sp = MR_dum3d_MetP
      call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
      Vy_meso_next_step_MetP_sp = MR_dum3d_MetP

      ! Now interpolate onto the first few steps of the trajectory level
      ! array
      if(TrajFlag.gt.0)then
        istep = 1
      else
        istep = MR_MetSteps_Total
      endif
      Step_Time_since1900(1) = MR_MetStep_Hour_since_baseyear(istep)
      do kk = 1,ntraj
        k  = OutLev_indx(kk)
        zc = 1.0_4-zfrac(kk)
        Vx_full(:,:,kk,istep) = Vx_meso_last_step_MetP_sp(:,:,k  )*zc + &
                                Vx_meso_last_step_MetP_sp(:,:,k+1)*zfrac(kk)
        Vy_full(:,:,kk,istep) = Vy_meso_last_step_MetP_sp(:,:,k  )*zc + &
                                Vy_meso_last_step_MetP_sp(:,:,k+1)*zfrac(kk)
      enddo
      if(TrajFlag.gt.0)then
        istep = 2
      else
        istep = MR_MetSteps_Total-1
      endif
      Step_Time_since1900(2) = MR_MetStep_Hour_since_baseyear(istep)
      do kk = 1,ntraj
        k  = OutLev_indx(kk)
        zc = 1.0_4-zfrac(kk)
        Vx_full(:,:,kk,istep) = Vx_meso_next_step_MetP_sp(:,:,k  )*zc + &
                                Vx_meso_next_step_MetP_sp(:,:,k+1)*zfrac(kk)
        Vy_full(:,:,kk,istep) = Vy_meso_next_step_MetP_sp(:,:,k  )*zc + &
                                Vy_meso_next_step_MetP_sp(:,:,k+1)*zfrac(kk)
      enddo

      ! And now fill in all the remaining steps
      if(MR_MetSteps_Total.ge.3)then
        do istep = 3,MR_MetSteps_Total
          if(TrajFlag.gt.0)then
            stepindx=istep
          else
            stepindx=MR_MetSteps_Total-istep+1
          endif
          Step_Time_since1900(istep) = MR_MetStep_Hour_since_baseyear(stepindx)
          MR_iMetStep_Now = MR_iMetStep_Now + 1
          ivar = 2 ! Vx
          call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
          Vx_meso_next_step_MetP_sp = MR_dum3d_MetP
          ivar = 3 ! Vy
          call MR_Read_3d_MetP_Variable(ivar,MR_iMetStep_Now+1)
          Vy_meso_next_step_MetP_sp = MR_dum3d_MetP
          do kk = 1,ntraj
            k  = OutLev_indx(kk)
            zc = 1.0_4-zfrac(kk)
            Vx_full(:,:,kk,stepindx) = Vx_meso_next_step_MetP_sp(:,:,k  )*zc + &
                                       Vx_meso_next_step_MetP_sp(:,:,k+1)*zfrac(kk)
            Vy_full(:,:,kk,stepindx) = Vy_meso_next_step_MetP_sp(:,:,k  )*zc + &
                                       Vy_meso_next_step_MetP_sp(:,:,k+1)*zfrac(kk)
          enddo
        enddo
      endif

      ! We now have the full x,y,z,vx,vy data needed from the Met file
      ! for the full forward/backward simulation

      ! Initialize the start coordinates
      x1(:) = inlon
      y1(:) = inlat
      t1    = Probe_StartHour
      it    = 1

      ! Assume an integration step of 1min and a max v of around 100m/s
      if(TrajFlag.ge.0)then
        dt = 1.0_8/60.0_8
      else
        dt = -1.0_8/60.0_8
      endif

      mstodeghr = 3600.0_8*360.0_8/(2.0_8*PI*RAD_EARTH*KM_2_M)

      if(TrajFlag.ge.0)then
        if(ntraj.ge.1)open(unit=21,file='ftraj1.dat')
        if(ntraj.ge.2)open(unit=22,file='ftraj2.dat')
        if(ntraj.ge.3)open(unit=23,file='ftraj3.dat')
        if(ntraj.ge.4)open(unit=24,file='ftraj4.dat')
        if(ntraj.ge.5)open(unit=25,file='ftraj5.dat')
        if(ntraj.ge.6)open(unit=26,file='ftraj6.dat')
        if(ntraj.ge.7)open(unit=27,file='ftraj7.dat')
        if(ntraj.ge.8)open(unit=28,file='ftraj8.dat')
        if(ntraj.ge.9)open(unit=29,file='ftraj9.dat')
      else
        if(ntraj.ge.1)open(unit=21,file='btraj1.dat')
        if(ntraj.ge.2)open(unit=22,file='btraj2.dat')
        if(ntraj.ge.3)open(unit=23,file='btraj3.dat')
        if(ntraj.ge.4)open(unit=24,file='btraj4.dat')
        if(ntraj.ge.5)open(unit=25,file='btraj5.dat')
        if(ntraj.ge.6)open(unit=26,file='btraj6.dat')
        if(ntraj.ge.7)open(unit=27,file='btraj7.dat')
        if(ntraj.ge.8)open(unit=28,file='btraj8.dat')
        if(ntraj.ge.9)open(unit=29,file='btraj9.dat')
      endif

      if(ntraj.ge.1)write(21,*)real(x1(1),kind=4),real(y1(1),kind=4)
      if(ntraj.ge.2)write(22,*)real(x1(2),kind=4),real(y1(2),kind=4)
      if(ntraj.ge.3)write(23,*)real(x1(3),kind=4),real(y1(3),kind=4)
      if(ntraj.ge.4)write(24,*)real(x1(4),kind=4),real(y1(4),kind=4)
      if(ntraj.ge.5)write(25,*)real(x1(5),kind=4),real(y1(5),kind=4)
      if(ntraj.ge.6)write(26,*)real(x1(6),kind=4),real(y1(6),kind=4)
      if(ntraj.ge.7)write(27,*)real(x1(7),kind=4),real(y1(7),kind=4)
      if(ntraj.ge.8)write(28,*)real(x1(8),kind=4),real(y1(8),kind=4)
      if(ntraj.ge.9)write(29,*)real(x1(9),kind=4),real(y1(9),kind=4)

      ! Find location of initial trajectory heights
      ! Get interpolation coefficients
      it = 1
      if(TrajFlag.gt.0)then
        tfrac = (t1-Step_Time_since1900(it))/MR_MetStep_Interval(it)
      else
        tfrac = -(t1-Step_Time_since1900(it))/MR_MetStep_Interval(it)
      endif

      ! integrate out Simtime_in_hours hours
      write(*,*)"Now integrating out ",abs(int(Simtime_in_hours/dt))," steps"
      it = 0
      do ti = 1,abs(int(Simtime_in_hours/dt))
        if(TrajFlag.gt.0)then
          ! Get the interval by assuming all MetStep_Intervals are the same
          iit = floor((t1-Step_Time_since1900(1))/MR_MetStep_Interval(1)) + 1
          tfrac = (t1-Step_Time_since1900(iit))/MR_MetStep_Interval(1)
        else
          ! Get the interval by assuming all MetStep_Intervals are the same
          iit = floor(-(t1-Step_Time_since1900(1))/MR_MetStep_Interval(1)) + 1
          tfrac = -(t1-Step_Time_since1900(iit))/MR_MetStep_Interval(1)
        endif

        if(iit.ne.it)then
          it = iit
           ! Get the change in velocity in m/s/hr
          dvxdt(:,:,:) = (Vx_full(:,:,:,it+1)-Vx_full(:,:,:,it))/MR_MetStep_Interval(it)
          dvydt(:,:,:) = (Vy_full(:,:,:,it+1)-Vy_full(:,:,:,it))/MR_MetStep_Interval(it)
        endif

        do kk = 1,ntraj
          ! Get current time and position indecies
          ix = floor((x1(kk)-x_submet_sp(1))/abs(dx_met_const)) + 1
          iy = floor((y1(kk)-y_submet_sp(1))/abs(dy_met_const)) + 1
            ! Skip over points that leave the domain
          if(ix.le.0.or.ix.ge.nx_comp)cycle
          if(iy.le.0.or.iy.ge.ny_comp)cycle

          ! Get the fractional position within the cell
          xfrac = (x1(kk)-x_submet_sp(ix))/dx_met_const
          yfrac = (y1(kk)-y_submet_sp(iy))/dx_met_const
          xc = 1.0_4-xfrac
          yc = 1.0_4-yfrac
          ! Build interpolation coefficients
          a1 = real(xc*yc,kind=4)
          a2 = real(xfrac*yc,kind=4)
          a3 = real(xfrac*yfrac,kind=4)
          a4 = real(yfrac*xc,kind=4)

          ! Corner velocities for current time
          vx1 = Vx_full(ix  ,iy  ,kk,it) + tfrac*dvxdt(ix  ,iy  ,kk)
          vx2 = Vx_full(ix+1,iy  ,kk,it) + tfrac*dvxdt(ix+1,iy  ,kk)
          vx3 = Vx_full(ix+1,iy+1,kk,it) + tfrac*dvxdt(ix+1,iy+1,kk)
          vx4 = Vx_full(ix  ,iy+1,kk,it) + tfrac*dvxdt(ix  ,iy+1,kk)
          vy1 = Vy_full(ix  ,iy  ,kk,it) + tfrac*dvydt(ix  ,iy  ,kk)
          vy2 = Vy_full(ix+1,iy  ,kk,it) + tfrac*dvydt(ix+1,iy  ,kk)
          vy3 = Vy_full(ix+1,iy+1,kk,it) + tfrac*dvydt(ix+1,iy+1,kk)
          vy4 = Vy_full(ix  ,iy+1,kk,it) + tfrac*dvydt(ix  ,iy+1,kk)

          ! Interpolate velocity onto current position and time (in m/s)
          vel_1(1) = (a1*vx1+a2*vx2+a3*vx3+a4*vx4)
          vel_1(2) = (a1*vy1+a2*vy2+a3*vy3+a4*vy4)
          !  now convert to deg/hr
          vel_1(1) = vel_1(1)*mstodeghr/sin((90.0_4-y1(kk))*DEG2RAD)
          vel_1(2) = vel_1(2)*mstodeghr

          ! Now advect via Forward Euler
          xstep = vel_1(1) * dt
          ystep = vel_1(2) * dt
          x_fin = x1(kk) + xstep
          y_fin = y1(kk) + ystep

          x1(kk) = x_fin
          y1(kk) = y_fin
        enddo

        t1 = t1 + dt
        if(mod(ti,60).eq.0)then
          do kk = 1,ntraj
            if(kk.eq.1)write(21,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.2)write(22,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.3)write(23,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.4)write(24,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.5)write(25,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.6)write(26,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.7)write(27,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.8)write(28,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            if(kk.eq.9)write(29,*)real(x1(kk),kind=4),real(y1(kk),kind=4)
            ! Check min/max of trajectories

            if(x1(kk).lt.lonmin)lonmin=x1(kk)
            if(x1(kk).gt.lonmax)lonmax=x1(kk)
            if(y1(kk).lt.latmin)latmin=y1(kk)
            if(y1(kk).gt.latmax)latmax=y1(kk)
          enddo
        endif
      enddo
      open(unit=40,file='map_range_traj.txt')
      write(40,*)lonmin,lonmax,latmin,latmax,inlon,inlat
      close(40)

      if(ntraj.ge.1)close(21)
      if(ntraj.ge.2)close(22)
      if(ntraj.ge.3)close(23)
      if(ntraj.ge.4)close(24)
      if(ntraj.ge.5)close(25)
      if(ntraj.ge.6)close(26)
      if(ntraj.ge.7)close(27)
      if(ntraj.ge.8)close(28)
      if(ntraj.ge.9)close(29)

       end subroutine Integrate_Trajectory

!##############################################################################
!##############################################################################

