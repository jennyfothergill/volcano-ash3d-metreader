!##############################################################################
!
!     MR_Read_Met_DimVars_netcdf
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_NC_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Read_Met_DimVars_netcdf

      use MetReader
      use netcdf
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: i, k
      real(kind=sp) :: xLL_fullmet
      real(kind=sp) :: yLL_fullmet
      real(kind=sp) :: xUR_fullmet
      real(kind=sp) :: yUR_fullmet

      integer :: ncid
      integer :: ivar,in_var_id,var_ndims
      integer :: i_dim,iivar
      integer :: xdim_id, ydim_id, tdim_id
      integer :: dimlen,maxdimlen
      integer :: nSTAT
      character(len=NF90_MAX_NAME)  :: invar,dimname
      integer :: var_xtype,var_id,idx
      integer :: xtype, length, attnum
      real(kind=dp), parameter :: tol = 1.0e-3_dp
      real(kind=dp) :: x_start,y_start

      !integer :: NC_version
      integer,dimension(:),allocatable :: var_dimIDs
      logical :: FoundOldDim
      logical :: IsPressureDimension
      real(kind=dp),dimension(:), allocatable :: dum1d_dp
      real(kind=sp),dimension(:), allocatable :: dum1d_sp
      character(len=31)  :: ustring
      logical :: IsTruncatedDim
      character(len=130)   :: infile

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Read_Met_DimVars_netcdf                  ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(MR_iwind.eq.5)then
        ! For the case where variables are in different files, we will just hard-code the
        ! grid geometry

        if (MR_iwindformat.eq.25) then
           ! NCEP/NCAR reanalysis 2.5 degree files 
           ! https://rda.ucar.edu/datasets/ds090.0
          maxdimlen            = 17
          nlev_coords_detected = 3
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 17
          nlevs_fullmet(2) = 12
          nlevs_fullmet(3) = 8
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          levs_code(2) = 2
          levs_code(3) = 2
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 17
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 925.0_sp, 850.0_sp, 700.0_sp, 600.0_sp, &
               500.0_sp, 400.0_sp, 300.0_sp, 250.0_sp, 200.0_sp, &
               150.0_sp, 100.0_sp,  70.0_sp,  50.0_sp, 30.0_sp, &
                20.0_sp,  10.0_sp /)
          z_inverted = .false.
          levs_fullmet_sp(1,1:17) = p_fullmet_sp(1:17)*Pressure_Conv_Fac
          levs_fullmet_sp(2,1:12) = p_fullmet_sp(1:12)*Pressure_Conv_Fac
          levs_fullmet_sp(3,1:8 ) = p_fullmet_sp(1:8 )*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 2
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx(30) = 3

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          nx_fullmet = 144
          ny_fullmet = 73
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          dx_met_const = 2.5_sp
          dy_met_const = 2.5_sp
          x_start =   0.0_dp
          y_start =  90.0_dp
          DO i = 0,nx_fullmet+1
            x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
          ENDDO
          x_inverted = .false.
          DO i = 1,ny_fullmet
            y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
          ENDDO
          y_inverted = .true.
          DO i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          ENDDO
          DO i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          ENDDO
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

          iwf25_scale_facs = 0.0_sp
          iwf25_offsets    = 0.0_sp
          iwf25_scale_facs(1)  = 1.0_sp    ; iwf25_offsets(1)  = 32066.0_sp   ! hgt
          iwf25_scale_facs(2)  = 0.01_sp   ; iwf25_offsets(2)  = 202.66_sp    ! uwnd
          iwf25_scale_facs(3)  = 0.01_sp   ; iwf25_offsets(3)  = 202.66_sp    ! vwnd
          iwf25_scale_facs(4)  = 0.001_sp  ; iwf25_offsets(4)  = 29.765_sp    ! omega
          iwf25_scale_facs(5)  = 0.01_sp   ; iwf25_offsets(5)  = 477.66_sp    ! air (temperature)
          iwf25_scale_facs(6)  = 1.0_sp    ; iwf25_offsets(6)  = 0.0_sp       ! level
          iwf25_scale_facs(20) = 10.0_sp   ; iwf25_offsets(20) = 327650.0_sp  ! pres (lcb)
          iwf25_scale_facs(21) = 10.0_sp   ; iwf25_offsets(21) = 327650.0_sp  ! pres (lct)
          iwf25_scale_facs(30) = 0.01_sp   ; iwf25_offsets(30) = 302.66_sp    ! rhum
          iwf25_scale_facs(31) = 1.0e-6_sp ; iwf25_offsets(31) = 0.032666_sp  ! shum
          iwf25_scale_facs(32) = 1.0e-6_sp ; iwf25_offsets(32) = 0.032666_sp  ! shum
          iwf25_scale_facs(44) = 1.0e-7_sp ; iwf25_offsets(44) = 0.0032765_sp ! prate
          iwf25_scale_facs(45) = 1.0e-7_sp ; iwf25_offsets(45) = 0.0031765_sp ! cprat

        elseif(MR_iwindformat.eq.26)then
         ! JRA-55 reanalysis 1.25 degree files 
         !  https://rda.ucar.edu/datasets/ds628.0/

          maxdimlen            = 37
          nlev_coords_detected = 2
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 37
          nlevs_fullmet(2) = 27
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          levs_code(2) = 2
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 37
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 975.0_sp, 950.0_sp, 925.0_sp, 900.0_sp, &
               875.0_sp, 850.0_sp, 825.0_sp, 800.0_sp, 775.0_sp, &
               750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
               500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
               250.0_sp, 225.0_sp, 200.0_sp, 175.0_sp, 150.0_sp, &
               125.0_sp, 100.0_sp,  70.0_sp,  50.0_sp,  30.0_sp, &
                20.0_sp,  10.0_sp,   7.0_sp,   5.0_sp,   3.0_sp, &
                 2.0_sp,   1.0_sp /)
          z_inverted = .true.

          levs_fullmet_sp(1,1:37) = p_fullmet_sp(1:37)*Pressure_Conv_Fac
          levs_fullmet_sp(2,1:27) = p_fullmet_sp(1:27)*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 1
          Met_var_zdim_idx( 5) = 1
          Met_var_zdim_idx(30) = 2

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          nx_fullmet = 288
          ny_fullmet = 145
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          dx_met_const = 1.25_sp
          dy_met_const = 2.0_sp
          x_start =   0.0_dp
          y_start =  90.0_dp
          DO i = 0,nx_fullmet+1
            x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
          ENDDO
          x_inverted = .false.
          DO i = 1,ny_fullmet
            y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
          ENDDO
          y_inverted = .true.
          DO i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          ENDDO
          DO i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          ENDDO
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

        elseif(MR_iwindformat.eq.27)then
          !  NOAA-CIRES reanalysis 2.0 degree files  :: ds131.2

          maxdimlen            = 24
          nlev_coords_detected = 2
          allocate(nlevs_fullmet(nlev_coords_detected))
          nlevs_fullmet(1) = 24
          nlevs_fullmet(2) = 19
          allocate(levs_code(nlev_coords_detected))
          levs_code(1) = 1
          levs_code(2) = 2
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          np_fullmet = 24
          allocate(p_fullmet_sp(np_fullmet))
          levs_fullmet_sp(:,:) = 0.0_sp
          p_fullmet_sp(1:np_fullmet) = &
            (/1000.0_sp, 950.0_sp, 900.0_sp, 850.0_sp, 800.0_sp, &
               750.0_sp, 700.0_sp, 650.0_sp, 600.0_sp, 550.0_sp, &
               500.0_sp, 450.0_sp, 400.0_sp, 350.0_sp, 300.0_sp, &
               250.0_sp, 200.0_sp, 150.0_sp, 100.0_sp,  70.0_sp, &
                50.0_sp,  30.0_sp,  20.0_sp,  10.0_sp /)
          z_inverted = .false.
          levs_fullmet_sp(1,1:24) = p_fullmet_sp(1:24)*Pressure_Conv_Fac
          levs_fullmet_sp(2,1:19) = p_fullmet_sp(1:19)*Pressure_Conv_Fac
          Met_var_zdim_idx( 1) = 1
          Met_var_zdim_idx( 2) = 1
          Met_var_zdim_idx( 3) = 1
          Met_var_zdim_idx( 4) = 2
          Met_var_zdim_idx( 5) = 1

          IsLatLon_MetGrid  = .true.
          IsGlobal_MetGrid  = .true.
          IsRegular_MetGrid = .true.
          nx_fullmet = 180
          ny_fullmet = 91
          allocate(x_fullmet_sp(0:nx_fullmet+1))
          allocate(y_fullmet_sp(ny_fullmet))
          allocate(MR_dx_met(nx_fullmet))
          allocate(MR_dy_met(ny_fullmet))
          dx_met_const = 2.0_sp
          dy_met_const = 2.0_sp
          x_start =   0.0_dp
          y_start =  90.0_dp
          DO i = 0,nx_fullmet+1
            x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=sp)
          ENDDO
          x_inverted = .false.
          DO i = 1,ny_fullmet
            y_fullmet_sp(i) = real(y_start - (i-1)*dy_met_const,kind=sp)
          ENDDO
          y_inverted = .true.
          DO i = 1,nx_fullmet
            MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
          ENDDO
          DO i = 1,ny_fullmet-1
            MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
          ENDDO
          MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

        elseif(MR_iwindformat.eq.29)then
         write(*,*) MR_iwindformat
        elseif(MR_iwindformat.eq.30)then
         write(*,*) MR_iwindformat
        else
          stop 1
        endif

      else  ! MR_iwind not equal to 5
        if(.not.MR_iwindformat.eq.50)then
          !---------------------------------------------------------------------------------
          ! Checking for dimension length and values for x,y,t,p
          !   Assume all files have the same format
          maxdimlen = 0
          if(MR_iwind.eq.5)then
 111        format(a50,a4,i4,a3)
            write(infile,111)trim(adjustl(MR_windfiles(1))), &
                             "hgt.",MR_Comp_StartYear,".nc"
            nSTAT = nf90_open(trim(ADJUSTL(infile)),NF90_NOWRITE,ncid)
          else
            infile = adjustl(trim(MR_windfiles(1)))
          endif
          nSTAT=nf90_open(adjustl(trim(infile)),NF90_NOWRITE, ncid)
          if(nSTAT.ne.NF90_NOERR) then
            write(MR_global_error,*)'MR ERROR: open NC file: ',nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: open NC file: ',nf90_strerror(nSTAT)
            write(MR_global_error,*)'Exiting'
            stop 1
          endif
          do ivar = 1,MR_MAXVARS
            if (.not.Met_var_IsAvailable(ivar)) cycle  ! Only look at variables that are available
            if (Met_var_ndim(ivar).ne.4) cycle         !  and only ones with a 'level' dimension
            invar = Met_var_NC_names(ivar)
            nSTAT = nf90_inq_varid(ncid,invar,in_var_id)  ! get the var_id for this named variable
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR WARNING: inq_varid: ',invar,nf90_strerror(nSTAT)
              write(MR_global_error,*)'  Cannot find variable ',invar
              write(MR_global_error,*)'  Setting Met_var_IsAvailable to .false.'
              write(MR_global_log  ,*)'MR WARNING: inq_varid: ',invar,nf90_strerror(nSTAT)
              Met_var_IsAvailable(ivar) = .false.
              cycle
            endif
            nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                      xtype = var_xtype, &
                      ndims = var_ndims)   ! get the number of dimensions
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              stop 1
            endif
            if (var_ndims.ne.Met_var_ndim(ivar))then
              write(MR_global_error,*)'MR ERROR: The actual number of dimensions differs from'
              write(MR_global_error,*)'          what is expected'
              write(MR_global_error,*)'      Variable : ',ivar,Met_var_NC_names(ivar)
              write(MR_global_error,*)'      Expected : ',Met_var_ndim(ivar)
              write(MR_global_error,*)'      Found    : ',var_ndims
              stop 1
            endif
            allocate(var_dimIDs(var_ndims))
            nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                      dimids = var_dimIDs(:var_ndims))
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              stop 1
            endif
    
            ! 3-d transient variables should be in the COORDS convention (time, level, y, x)
            ! if ivar = 1 (Geopotential Height), then get the info on x,y and t too
            if(ivar.eq.1)then
              i_dim = 1  ! get x info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                stop 1
              endif
              if(index(dimname,Met_dim_names(4)).ne.0)then
                nx_fullmet = dimlen
                xdim_id    = var_dimIDs(i_dim)
              endif
              nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              allocate(x_fullmet_sp(0:nx_fullmet+1))
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                x_fullmet_sp(1:nx_fullmet) = dum1d_sp(1:nx_fullmet)
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                x_fullmet_sp(1:nx_fullmet) = real(dum1d_dp(1:nx_fullmet),kind=sp)
                deallocate(dum1d_dp)
              else
                write(MR_global_error,*)'MR ERROR: Cannot recognize variable type for x'
                stop 1
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                if(nSTAT.ne.NF90_NOERR) then
                  write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  stop 1
                endif
                if(index(ustring,'km').gt.0.or.&
                   index(ustring,'kilo').gt.0)then
                  ! This is a projected grid
                  IsLatLon_MetGrid  = .false.
                elseif(index(ustring,'deg').gt.0)then
                  ! This is a lon/lat grid
                  IsLatLon_MetGrid  = .true.
                else
                  write(MR_global_error,*)"MR ERROR: Cannot determine if the grid is lon/lat or projected"
                  stop 1
                endif
              endif
              ! Finally, check for orientation
              if(x_fullmet_sp(1).lt.x_fullmet_sp(2))then
                x_inverted = .false.
              else
                x_inverted = .true.
              endif
              allocate(MR_dx_met(nx_fullmet))
              x_start = x_fullmet_sp(1)
              dx_met_const = x_fullmet_sp(2)-x_fullmet_sp(1)
              x_fullmet_sp(0)            = x_fullmet_sp(1)          - dx_met_const
              x_fullmet_sp(nx_fullmet+1) = x_fullmet_sp(nx_fullmet) + dx_met_const
              do i = 1,nx_fullmet
                MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
              enddo 
 
              i_dim = 2  ! get y info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                stop 1
              endif
              if(index(dimname,Met_dim_names(3)).ne.0)then
                ny_fullmet = dimlen
                ydim_id    = var_dimIDs(i_dim)
              endif

              nSTAT = nf90_inq_varid(ncid,dimname,var_id) ! get the variable associated with this dim
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              allocate(y_fullmet_sp(1:ny_fullmet))
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                y_fullmet_sp(1:ny_fullmet) = dum1d_sp(1:ny_fullmet)
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                y_fullmet_sp(1:ny_fullmet) = real(dum1d_dp(1:ny_fullmet),kind=sp)
                deallocate(dum1d_dp)
              else
                write(MR_global_error,*)'MR ERROR: Cannot recognize variable type for x'
                stop 1
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR WARNING: cannot file units ',dimname,nf90_strerror(nSTAT)
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                if(nSTAT.ne.NF90_NOERR) then
                  write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  stop 1
                endif
                if(index(ustring,'km').gt.0.or.&
                   index(ustring,'kilo').gt.0)then
                  ! This is a projected grid
                  IsLatLon_MetGrid  = .false.
                elseif(index(ustring,'deg').gt.0)then
                  ! This is a lon/lat grid
                  IsLatLon_MetGrid  = .true.
                else
                  write(MR_global_error,*)"MR ERROR: Cannot determine if the grid is lon/lat or projected"
                  stop 1
                endif
              endif
              ! check for orientation
              if(y_fullmet_sp(1).lt.y_fullmet_sp(2))then
                y_inverted = .false.
              else
                y_inverted = .true.
              endif
              allocate(MR_dy_met(ny_fullmet))
              y_start = y_fullmet_sp(1)
              dy_met_const = y_fullmet_sp(2)-y_fullmet_sp(1)
              do i = 1,ny_fullmet-1
                MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
              enddo

              ! We need to check if this is a regular grid
              IsRegular_MetGrid = .true.
              do i = 1,nx_fullmet-1
                if(abs(MR_dx_met(i+1)-MR_dx_met(i)).gt.tol*MR_dx_met(i))then
                  IsRegular_MetGrid = .false.
                endif
              enddo
              do i = 1,ny_fullmet-1
                if(abs(MR_dy_met(i+1)-MR_dy_met(i)).gt.tol*MR_dy_met(i))then
                  IsRegular_MetGrid = .false.
                endif
              enddo

              i_dim = 4  ! get t info
              nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                           name =  dimname, &
                           len = dimlen)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
                stop 1
              endif
              if(index(dimname,Met_dim_names(1)).ne.0)then
                nt_fullmet = dimlen
                tdim_id    = var_dimIDs(i_dim)
              endif
  
            endif
            ! Now checking level coordinates (pressure, height, depth); third dimension
            i_dim = 3
            nSTAT = nf90_inquire_dimension(ncid,var_dimIDs(i_dim), &
                         name =  dimname, & 
                         len = dimlen)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
              stop 1
            endif
            if(index(dimname,'lev').ne.0.or.&
               index(dimname,'isobaric').ne.0.or.&
               index(dimname,'pressure').ne.0.or.&
               index(dimname,'height').ne.0.or.&
               index(dimname,'depth').ne.0.or.&
               index(dimname,'lv_ISBL1').ne.0.or.&
               index(dimname,'bottom_top').ne.0.or.&
               index(dimname,'bottom_top_stag').ne.0.or.&
               index(dimname,'soil_layers_stag').ne.0)then
              ! Log this level coordinate if it is the first
              if (nlev_coords_detected.eq.0)then
                nlev_coords_detected = nlev_coords_detected + 1
                Met_var_zdim_idx(ivar)  = nlev_coords_detected
                Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                maxdimlen = dimlen
              else
                ! Otherwise, check if this level coordinate has already been logged
                FoundOldDim = .false.
                do iivar = 1,ivar-1
                  if (Met_var_zdim_ncid(iivar).eq.var_dimIDs(i_dim))then
                    FoundOldDim = .true.
                    Met_var_zdim_idx(ivar)  = Met_var_zdim_idx(iivar)
                    Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                    exit
                  endif
                enddo
                if(.not.FoundOldDim)then
                  nlev_coords_detected = nlev_coords_detected + 1
                  Met_var_zdim_idx(ivar)  = nlev_coords_detected
                  Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)
                  if (maxdimlen.lt.dimlen) maxdimlen = dimlen
                endif
              endif
            else
              write(MR_global_error,*)'MR ERROR: level coordinate is not in pos. 3 for ',invar
              write(MR_global_log  ,*)'MR ERROR: level coordinate is not in pos. 3 for ',invar
              stop 1
            endif
            ! tidy up
            deallocate(var_dimIDs)
          enddo ! ivar
  
          ! We have all the level dimension names and dim_ids; now we need to get the sizes
          allocate(nlevs_fullmet(nlev_coords_detected))
          allocate(levs_code(nlev_coords_detected))
          allocate(levs_fullmet_sp(nlev_coords_detected,maxdimlen))
          do ivar = 1,MR_MAXVARS
            ! Check if this variable has a z-dimension (pressure, height, depth, etc.)
            if(Met_var_zdim_ncid(ivar).gt.0)then
              ! log the length of the dimension for this level coordinat
              nSTAT = nf90_inquire_dimension(ncid,Met_var_zdim_ncid(ivar), &
                         name =  dimname, &
                         len = dimlen)
              idx = Met_var_zdim_idx(ivar)
              nlevs_fullmet(idx) = dimlen
              ! Now inquire and populate the dimension variable info
              nSTAT = nf90_inq_varid(ncid,dimname,var_id)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_varid ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              ! Check if we need to read into a float or a double
              nSTAT = nf90_inquire_variable(ncid, var_id, dimname, xtype = var_xtype)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_variable: ',dimname,nf90_strerror(nSTAT)
                stop 1
              endif
              if(var_xtype.eq.NF90_FLOAT)then
                allocate(dum1d_sp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_sp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = dum1d_sp(1:nlevs_fullmet(idx))
                deallocate(dum1d_sp)
              elseif(var_xtype.eq.NF90_DOUBLE)then
                allocate(dum1d_dp(dimlen))
                nSTAT = nf90_get_var(ncid,var_id,dum1d_dp, &
                       start = (/1/),count = (/dimlen/))
                if(nSTAT.ne.NF90_NOERR)then
                  write(MR_global_error,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_var ',dimname,nf90_strerror(nSTAT)
                  stop 1
                endif
                ! copy to local variable
                levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = real(dum1d_dp(1:nlevs_fullmet(idx)),kind=sp)
                deallocate(dum1d_dp)
              endif
              ! Check the units
              nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                             "units",xtype, length, attnum)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR WARNING: cannot find dim units ',dimname,nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR WARNING: cannot find dim units ',dimname,nf90_strerror(nSTAT)
                IsPressureDimension = .false.
                stop 1
              else
                nSTAT = nf90_get_att(ncid, var_id,"units",ustring)
                if(nSTAT.ne.NF90_NOERR) then
                  write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                  stop 1
                endif
                ! Note: the variables below are single-valued, not arrays on ivar
                !       If a pressure is found, the assumption here is that all pressure coordinates will
                !       be given in the same units (hPa or Pa) and the same orientations (bot to top, or
                !       inverted.
                if(index(ustring,'Pa').gt.0.or.&
                   index(ustring,'millibar').gt.0)then
                  ! This is a pressure level
                  IsPressureDimension = .true.
                  if(index(ustring,'hPa').gt.0.or.&
                     index(ustring,'millibar').gt.0)then
                    Pressure_Conv_Fac = 100.0_sp
                  else
                    Pressure_Conv_Fac = 1.0_sp
                  endif
                else
                  IsPressureDimension = .false.
                endif
              endif
              
              ! Finally, check for orientation
              if(IsPressureDimension)then
                if(levs_fullmet_sp(idx,1).lt.levs_fullmet_sp(idx,2))then
                  z_inverted = .true.
                else
                  z_inverted = .false.
                endif
              endif
            endif
          enddo  ! ivar
          ! Close file
          nSTAT = nf90_close(ncid)
          if(nSTAT.ne.NF90_NOERR)then
             write(MR_global_error,*)'MR ERROR: close file: ',nf90_strerror(nSTAT)
             write(MR_global_log  ,*)'MR ERROR: close file: ',nf90_strerror(nSTAT)
             stop 1
          endif
        endif ! MR_iwind.eq.5
        !-----------------------------------------

        write(MR_global_production,*)" Found these levels"
        write(MR_global_production,*)"  VaribleID    LevelIdx       dimID      length"
        do ivar = 1,MR_MAXVARS
          if (Met_var_IsAvailable(ivar)) &
              write(MR_global_production,*)ivar,Met_var_zdim_idx(ivar),Met_var_zdim_ncid(ivar),&
                                           nlevs_fullmet(Met_var_zdim_idx(ivar))
        enddo
        ! Now invert if necessary and convert to Pa
        allocate(p_fullmet_sp(maxdimlen))
        do idx = 1,nlev_coords_detected
          if(z_inverted)then
            do i = 1,nlevs_fullmet(idx)
              p_fullmet_sp(nlevs_fullmet(idx)+1-i) = levs_fullmet_sp(idx,i)*Pressure_Conv_Fac
            enddo
          else
            p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))*Pressure_Conv_Fac
          endif
          levs_fullmet_sp(idx,:) = 0.0_sp
          levs_fullmet_sp(idx,1:nlevs_fullmet(idx)) = p_fullmet_sp(1:nlevs_fullmet(idx))
        enddo
        deallocate(p_fullmet_sp)
  
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
  
        ! Now assign these levels to the working arrays
        np_fullmet    = nlevs_fullmet(Met_var_zdim_idx(1))  ! Assign fullmet the length of H,U,V
  
        ! Geopotential
        if(.not.allocated(p_fullmet_sp))  allocate(p_fullmet_sp(np_fullmet))
        idx = Met_var_zdim_idx(1)
        p_fullmet_sp(1:nlevs_fullmet(idx)) = levs_fullmet_sp(idx,1:nlevs_fullmet(idx))
        MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)/100.0_sp)
  
      endif ! MR_iwindformat.eq.50
      !---------------------------------------------------------------------------------

      if(MR_iwindformat.eq.0)then
        ! Template windfile (example for nam198)
        !  Need to populate
        call MR_Set_Met_Dims_Template_netcdf
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

      end subroutine MR_Read_Met_DimVars_netcdf

!##############################################################################


!##############################################################################
!
!     MR_Get_WRF_grid
!
!     Called once from MR_Read_Met_DimVars_netcdf
!
!     This subroutine reads the variable and dimension IDs, and fills the
!     coordinate dimension variables, just like MR_Read_Met_DimVars_netcdf,
!     but WRF files are a bit more complicated than most since the grid/projection can
!     vary and state variables (U,V,W,Gph) are each on grids staggered in the
!     relevant direction.
!
!     After this subroutine completes, the following variables will be set:
!       All the projection parameters of NWP grid
!       Met_dim_names, Met_var_NC_names, Met_var_conversion_factor, Met_var_IsAvailable
!       The lengths of all the dimensions of the file
!       p_fullmet_sp (converted to Pa)
!       x_fullmet_sp, y_fullmet_sp
!       IsLatLon_MetGrid, IsGlobal_MetGrid, IsRegular_MetGrid 
!
!##############################################################################

      subroutine MR_Get_WRF_grid

      use MetReader
      use netcdf
      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: nSTAT
      integer :: ncid

      character(len = nf90_max_name) :: name_dum

      integer :: t_dim_id          = 0 ! x or lon
      integer :: x_dim_id          = 0 ! x or lon
      integer :: y_dim_id          = 0 ! y or lat
      integer :: z_dim_id          = 0 ! x or lon

      integer :: lon_var_id
      integer :: lat_var_id
      integer :: PB_var_id
      integer :: Ppert_var_id

      integer :: Map_Proj
      real(kind=sp) :: WRF_dx,WRF_dy
      real(kind=sp) :: Cen_lat,Stand_Lon,Truelat1,Truelat2
      real(kind=sp),dimension(:,:,:)  ,allocatable :: dum3d_sp
      real(kind=sp),dimension(:,:,:,:),allocatable :: dum4d_sp
      integer :: i

      real(kind=dp) :: x_start,y_start

      real(kind=dp) :: lat_in,lon_in

        ! MAP_PROJ - Model projection  1=Lambert, 2=polar stereographic, 
        !                              3=mercator, 6=lat-lon

        ! First set spatial (x/y) grid
        ! Open first windfile and assume all grids are the same

      write(MR_global_info,*)"About to open first WRF file : ",MR_windfiles(1)
      nSTAT=nf90_open(adjustl(trim(MR_windfiles(1))),NF90_NOWRITE, ncid)
      if(nSTAT.ne.NF90_NOERR) then
        write(MR_global_error,*)'MR ERROR: open WRF file: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: open WRF file: ',nf90_strerror(nSTAT)
        write(MR_global_error,*)'Exiting'
        stop 1
      endif
      
      Met_dim_names(1) = "Time"       ! time
      Met_dim_names(2) = "bottom_top"      ! pressure (24 levels 10 -> 1000)
      Met_dim_names(3) = "south_north"      ! y        (90.0 -> -90.0)
      Met_dim_names(4) = "west_east"      ! x        (0.0 -> 378.0)
      Met_dim_names(5) = "bottom_top_stag"      ! Stag in Z (pressure coordinate for Vz)
      Met_dim_names(6) = "west_east_stag"       ! Stag in x
      Met_dim_names(7) = "south_north_stag"     ! Stag in y

      ! Get dim ids and sizes
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(1),t_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_dimid Time: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_dimid Time: ',nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,name=name_dum,len=nt_fullmet)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Inquire_Dimension Time: ', &
                             nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension Time: ', &
                             nf90_strerror(nSTAT)
        stop 1
      endif 
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(4),x_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_dimid x: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_dimid x: ',nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_Inquire_Dimension(ncid,x_dim_id,name=name_dum,len=nx_fullmet)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Inquire_Dimension x: ', &
                             nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension x: ', &
                             nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(3),y_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_dimid y: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_dimid y: ',nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_Inquire_Dimension(ncid,y_dim_id,name=name_dum,len=ny_fullmet)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Inquire_Dimension y: ', &   
                             nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension y: ', &
                             nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(2),z_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_dimid z: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_dimid z: ',nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_Inquire_Dimension(ncid,z_dim_id,name=name_dum,len=neta_fullmet)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Inquire_Dimension z: ', &
                             nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension z: ', &
                             nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ", Map_Proj)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: get_att MAP_PROJ: ', &
                             nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: get_att MAP_PROJ: ', &
                             nf90_strerror(nSTAT)
        stop 1
      endif

      if(Map_Proj.eq.1)then
         ! Lambert
         !   truelat1
         !   truelat2 (optional)
         !   stand_lon
         !proj +proj=lcc +lon_0=-175.0 +lat_0=55.0 +lat_1=50.0 +lat_2=60.0 +R=6371.229

        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DX", WRF_dx)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att DX: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att DX: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DY", WRF_dy)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att DY: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att DY: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT", Cen_Lat)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att CEN_LAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att CEN_LAT: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", Stand_Lon)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att STAND_LON: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att STAND_LON: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", Truelat1)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att TRUELAT1: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att TRUELAT1: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", Truelat2)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att TRUELAT2: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att TRUELAT2: ',nf90_strerror(nSTAT)
          stop 1
        endif

          ! convert dx, dy to km
        IsRegular_MetGrid = .true.
        dx_met_const = WRF_dx*1.0e-3_4
        dy_met_const = WRF_dy*1.0e-3_4

        ! Projected grids have Lon and Lat provided as 2d fields
        allocate(Met_Proj_lat(nx_fullmet,ny_fullmet))
        allocate(Met_Proj_lon(nx_fullmet,ny_fullmet))
        allocate(dum3d_sp(nx_fullmet,ny_fullmet,1))

        nSTAT = nf90_inq_varid(ncid,"XLONG",lon_var_id)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_varid XLONG: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: inq_varid XLONG: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_inq_varid(ncid,"XLAT",lat_var_id)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_varid XLAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: inq_varid XLAT: ',nf90_strerror(nSTAT)
          stop 1
        endif

        nSTAT = nf90_get_var(ncid,lon_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_var XLONG: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_var XLONG: ',nf90_strerror(nSTAT)
          stop 1
        endif
           Met_Proj_lon(:,:) = dum3d_sp(:,:,1)
        nSTAT = nf90_get_var(ncid,lat_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_var XLAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_var XLAT: ',nf90_strerror(nSTAT)
          stop 1
        endif
           Met_Proj_lat(:,:) = dum3d_sp(:,:,1)

        ! In the example WRF files, x and y projected values are not actually
        ! provided, so we recreate them here using the coordinates if the LL
        ! point of the Lon/Lat grid
        allocate(x_fullmet_sp(nx_fullmet))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))

        lon_in = real(Met_Proj_lon(1,1),kind=8)
        lat_in = real(Met_Proj_lat(1,1),kind=8)

          ! Setting the projection parameters as libprojection.a expects
        Met_iprojflag = 4  
        Met_lam0   = real(Stand_Lon,kind=8)
        Met_phi0   = real(Cen_Lat,kind=8)
        Met_phi1   = real(Truelat1,kind=8)
        Met_phi2   = real(Truelat2,kind=8)
        Met_k0     = real(1.0,kind=8)
        Met_Re     = PJ_radius_earth
        call PJ_proj_for(lon_in,lat_in, &
                       Met_iprojflag,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
        do i = 1,nx_fullmet
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=4)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=4)
        enddo
        do i = 1,nx_fullmet-1
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        MR_dx_met(nx_fullmet)    = MR_dx_met(nx_fullmet-1)
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(Map_Proj.eq.2)then
        ! Polar Stereographic
        !   truelat1
        !   stand_lon
        write(MR_global_info,*)&
         "WRF: MAP_PROJ=2 : Polar Stereographic : Not implemented"
        stop 1
      elseif(Map_Proj.eq.3)then
        ! Mercator
        !  truelat1
         ! stand_lon
         !proj +proj=merc 

        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DX", WRF_dx)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att DX: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att DX: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "DY", WRF_dy)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att DY: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att DY: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT", Cen_Lat)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att CEN_LAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att CEN_LAT: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", Stand_Lon)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att STAND_LON: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att STAND_LON: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1", Truelat1)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att TRUELAT1: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att TRUELAT1: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2", Truelat2)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_att TRUELAT2: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_att TRUELAT2: ',nf90_strerror(nSTAT)
          stop 1
        endif

          ! convert dx, dy to km
        IsRegular_MetGrid = .true.
        dx_met_const = WRF_dx*1.0e-3_4
        dy_met_const = WRF_dy*1.0e-3_4

        ! Projected grids have Lon and Lat provided as 2d fields
        allocate(Met_Proj_lat(nx_fullmet,ny_fullmet))
        allocate(Met_Proj_lon(nx_fullmet,ny_fullmet))
        allocate(dum3d_sp(nx_fullmet,ny_fullmet,1))

        nSTAT = nf90_inq_varid(ncid,"XLONG",lon_var_id)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_varid XLONG: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: inq_varid XLONG: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_inq_varid(ncid,"XLAT",lat_var_id)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_varid XLAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: inq_varid XLAT: ',nf90_strerror(nSTAT)
          stop 1
        endif

        nSTAT = nf90_get_var(ncid,lon_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_var XLONG: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_var XLONG: ',nf90_strerror(nSTAT)
          stop 1
        endif
           Met_Proj_lon(:,:) = dum3d_sp(:,:,1)
        nSTAT = nf90_get_var(ncid,lat_var_id,dum3d_sp, &
               start = (/1,1,1/),count = (/nx_fullmet,ny_fullmet,1/))
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: get_var XLAT: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: get_var XLAT: ',nf90_strerror(nSTAT)
          stop 1
        endif
           Met_Proj_lat(:,:) = dum3d_sp(:,:,1)

        ! In the example WRF files, x and y projected values are not actually
        ! provided, so we recreate them here using the coordinates if the LL
        ! point of the Lon/Lat grid
        allocate(x_fullmet_sp(nx_fullmet))
        allocate(y_fullmet_sp(ny_fullmet))
        allocate(MR_dx_met(nx_fullmet))
        allocate(MR_dy_met(ny_fullmet))

        lon_in = real(Met_Proj_lon(1,1),kind=8)
        lat_in = real(Met_Proj_lat(1,1),kind=8)

          ! Setting the projection parameters as libprojection.a expects
        Met_iprojflag = 5
        Met_lam0   = real(Stand_Lon,kind=8)
        Met_phi0   = real(Cen_Lat,kind=8)
        Met_phi1   = real(Truelat1,kind=8)
        Met_phi2   = real(Truelat2,kind=8)
        Met_k0     = real(1.0,kind=8)
        Met_Re     = PJ_radius_earth
        call PJ_proj_for(lon_in,lat_in, &
                       Met_iprojflag,Met_lam0,Met_phi0,Met_phi1,Met_phi2,Met_k0,Met_Re, &
                       x_start,y_start)
        do i = 1,nx_fullmet
          x_fullmet_sp(i) = real(x_start + (i-1)*dx_met_const,kind=4)
        enddo
        do i = 1,ny_fullmet
          y_fullmet_sp(i) = real(y_start + (i-1)*dy_met_const,kind=4)
        enddo
        do i = 1,nx_fullmet-1
          MR_dx_met(i) = x_fullmet_sp(i+1)-x_fullmet_sp(i)
        enddo
        MR_dx_met(nx_fullmet)    = MR_dx_met(nx_fullmet-1)
        do i = 1,ny_fullmet-1
          MR_dy_met(i) = y_fullmet_sp(i+1)-y_fullmet_sp(i)
        enddo
        MR_dy_met(ny_fullmet)    = MR_dy_met(ny_fullmet-1)

      elseif(Map_Proj.eq.6)then
        ! Lon-Lat or cylindrical equidistant
        !   pole_lat
        !   pole_lon
        !   stand_lon
        write(MR_global_info,*)"WRF: MAP_PROJ=6 : Lon-Lat : Not implemented"
        stop 1
      else
        write(MR_global_info,*)&
         "The MAP_PROJ global attribute is either not present or is"
        write(MR_global_info,*)"not a recognized projection."
        stop 1
      endif

      ! Now setting up pressure coordinate
      ! WRF data are provided on eta levels instead of pressure level
      ! Fortunately, it will usually be adequate to pretend that they are on
      ! pressure levels since Ash3d interpolates onto a z-grid using GPH.
      ! The "atmosphere" module, will need access to physical pressure and will
      ! need a special case for WRF files
      np_fullmet    = neta_fullmet
      !np_fullmet_Vz = neta_fullmet ! Vz is actually on a staggered grid, but
                                   ! will be interpolated onto the non-staggered grid
      !np_fullmet_RH = neta_fullmet
      allocate(p_fullmet_sp(np_fullmet))
      !allocate(p_fullmet_Vz_sp(np_fullmet_Vz))
      !allocate(p_fullmet_RH_sp(np_fullmet_RH))
      allocate(dum4d_sp(nx_fullmet,ny_fullmet,np_fullmet,1))

      ! To populate a place-holder p_fullmet_sp, read the full pressure grid and
      ! copy a representative column
      nSTAT = nf90_inq_varid(ncid,"PB",PB_var_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_varid PB: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_varid PB: ',nf90_strerror(nSTAT)
        stop 1
      endif

      nSTAT = nf90_get_var(ncid,PB_var_id,dum4d_sp, &
             start = (/1,1,1,1/),count = (/nx_fullmet,ny_fullmet,neta_fullmet,1/))
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: get_var PB: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: get_var PB: ',nf90_strerror(nSTAT)
        stop 1
      endif
         p_fullmet_sp(:) = dum4d_sp(1,1,:,1)
      nSTAT = nf90_inq_varid(ncid,"P",Ppert_var_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_varid P: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_varid P: ',nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_get_var(ncid,Ppert_var_id,dum4d_sp, &
             start = (/1,1,1,1/),count = (/nx_fullmet,ny_fullmet,neta_fullmet,1/))
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: get_var P: ',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: get_var P: ',nf90_strerror(nSTAT)
        stop 1
      endif
         p_fullmet_sp(:) = p_fullmet_sp(:) + dum4d_sp(1,1,:,1)
      !p_fullmet_Vz_sp = p_fullmet_sp
      !p_fullmet_RH_sp = p_fullmet_sp
      MR_Max_geoH_metP_predicted = MR_Z_US_StdAtm(p_fullmet_sp(np_fullmet)*0.01_sp) 

      !p_fullmet_sp    = p_fullmet_sp    * 100.0_sp   ! convert from hPa to Pa
      !p_fullmet_Vz_sp = p_fullmet_Vz_sp * 100.0_sp   ! convert from hPa to Pa
      !p_fullmet_RH_sp = p_fullmet_RH_sp * 100.0_sp   ! convert from hPa to Pa

       x_inverted = .false.
       y_inverted = .false.
       z_inverted = .false.

       ! Close file
       nSTAT = nf90_close(ncid)
       if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: close WRF file: ',nf90_strerror(nSTAT)
          write(MR_global_log  ,*)'MR ERROR: close WRF file: ',nf90_strerror(nSTAT)
          stop 1
       endif

       nlev_coords_detected = 1
       allocate(nlevs_fullmet(nlev_coords_detected))
       nlevs_fullmet(1) = np_fullmet
       allocate(levs_code(nlev_coords_detected))
       levs_code(1) = 1
       allocate(levs_fullmet_sp(nlev_coords_detected,np_fullmet))
       levs_fullmet_sp(1,1:np_fullmet) = p_fullmet_sp(1:np_fullmet)

       Met_var_zdim_idx(:)  = 1
       !Met_var_zdim_ncid(ivar) = var_dimIDs(i_dim)


       end subroutine MR_Get_WRF_grid

!##############################################################################
!
!     MR_Read_Met_Times_netcdf
!
!     Called once from MR_Read_Met_DimVars 
!
!     This subroutine opens each netcdf file and determine the time of each
!     time step of each file in the number of hours since MR_BaseYear.
!     In most cases, the length of the time variable (nt_fullmet) will be 
!     read directly from the file and overwritten (is was set in MR_Read_Met_DimVars_netcdf
!     above).
!
!     After this subroutine completes, the following variables will be set:
!       MR_windfile_starthour(MR_iwindfiles)
!       MR_windfile_stephour(MR_iwindfiles,nt_fullmet)
!
!##############################################################################

      subroutine MR_Read_Met_Times_netcdf

      use MetReader
      use netcdf

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer :: iw,iws
      integer :: itstart_year,itstart_month
      integer :: itstart_day
      real(kind=sp) :: filestart_hour

      integer :: itstart_hour,itstart_min,itstart_sec

      integer :: nSTAT
      integer :: ncid
      integer :: time_var_id = 0
      integer :: reftime_var_id
      integer :: t_dim_id
      real(kind=sp),dimension(:),allocatable :: filetime_in_sp
      character(len=19) :: Timestr_WRF
      integer :: filetime_in_int

      integer            :: var_xtype
      character(len=NF90_MAX_NAME)  :: invar
      integer            :: xtype, length, attnum
      character(len=20)  :: tstring
      character(len=31)  :: tstring2
      real(kind=8)       :: HS_hours_since_baseyear !,HS_HourOfDay
      character(len=130) :: dumstr
      integer            :: iwstep
      logical            :: TimeHasUnitsAttr = .false.
      integer            :: i,ii
      real(kind=dp),dimension(:), allocatable :: dum1d_dp
      real(kind=sp),dimension(:), allocatable :: dum1d_sp
      integer(kind=4),dimension(:), allocatable :: dum1d_int4
      integer,dimension(8) :: values
      integer              :: Current_Year,nt_tst
      character(len=130)   :: infile

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Read_Met_Times_netcdf                    ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(.not.Met_dim_IsAvailable(1))then
        write(MR_global_error,*)"MR ERROR: Time dimension is required and not listed"
        write(MR_global_error,*)"          in template windfile specification file."
        stop 1
      endif

      allocate(MR_windfile_starthour(MR_iwindfiles))
      if(MR_iwind.eq.5)then
      !if(MR_iwindformat.eq.25.or.MR_iwindformat.eq.27)then
        ! Here the branch for when MR_iwindformat = 25 or 27
        ! First copy path read in to slot 2
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
        ! Note: The nt_fullmet given above is the expected number based on a complete year.
        !       If MR_Comp_StartYear is the current year, then we will not have this many
        !       time steps available.  Double-check
        ! Getting current year
        call date_and_time(VALUES=values)
        Current_Year = values(1)
        ! Now getting the actual nt of the file; trying the hgt file
        if(MR_iwindformat.eq.25)then
 111      format(a50,a4,i4,a3)
          write(infile,111)trim(adjustl(MR_windfiles(1))), &
                           "hgt.",MR_Comp_StartYear,".nc"
        elseif(MR_iwindformat.eq.26)then
 112      format(a50,a17,i4,a20)
          write(infile,112)trim(adjustl(MR_windfiles(1))), &
                           "anl_p125.007_hgt.",MR_Comp_StartYear,"060100_2018063018.nc"
        else
          write(*,*)"Hard code file name for time check"
        endif
        nSTAT = nf90_open(trim(ADJUSTL(infile)),NF90_NOWRITE,ncid)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: nf90_open: ',nf90_strerror(nSTAT)
          write(MR_global_error,*)"    Could not open file: ",trim(ADJUSTL(infile))
          write(MR_global_log  ,*)'MR ERROR: nf90_open: ',nf90_strerror(nSTAT)
          write(MR_global_error,*)'Exiting'
          stop 1
        endif
        nSTAT = nf90_inq_dimid(ncid,Met_dim_names(1),t_dim_id)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
          write(MR_global_error,*)"    Could not find dimension: ",Met_dim_names(1)
          write(MR_global_log  ,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_tst)
        if(nSTAT.ne.NF90_NOERR)then
          write(MR_global_error,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
          write(MR_global_error,*)"    Could not dimension length: "
          write(MR_global_log  ,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
          stop 1
        endif
        nSTAT = nf90_close(ncid)
        if(MR_Comp_StartYear.lt.Current_Year.and.nt_tst.lt.nt_fullmet)then
          write(MR_global_info,*)"WARNING:  The NCEP files are for an archived year yet are incomplete."
          write(MR_global_info,*)"          To get the complete year, run the script "
          write(MR_global_info,*)"            autorun_scripts/get_NCEP_50YearReanalysis.sh",MR_Comp_StartYear
          write(MR_global_info,*)"          Steps available = ",nt_tst
          write(MR_global_info,*)"          Hours into year = ",(nt_tst-1)*6
        endif
        nt_fullmet = nt_tst
        MR_windfiles_nt_fullmet(1)=nt_fullmet
        MR_windfiles_nt_fullmet(2)=nt_fullmet  ! Note: we don't care if the next
                                               !       year is a leap year since
                                               !       the simulation will never
                                               !       be long enough to matter.
        allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))

          ! the interval for both iwf25 and iwf27 is 6 hours
        do iwstep = 1,nt_fullmet
          MR_windfile_stephour(:,iwstep) = (iwstep-1)*6.0_dp
        enddo
      elseif(MR_iwindformat.eq.31)then
        ! Here's the branch for the Catania files
        do iw = 1,MR_iwindfiles
          nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: nf90_open to read header:', &
                           nf90_strerror(nSTAT)
            write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
            write(MR_global_error,*)'Exiting'
            stop 1
          endif
          if(iw.eq.1)then
            nSTAT = nf90_inq_dimid(ncid,Met_dim_names(1),t_dim_id)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)"    Could not find dimension: ",Met_dim_names(1)
              write(MR_global_log  ,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
              stop 1
            endif
            nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_fullmet)
            if(nSTAT.ne.NF90_NOERR) then
              write(MR_global_error,*)'MR ERROR: Inquire_Dimension time: ', &
                                 nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension time: ', &
                                 nf90_strerror(nSTAT)
              stop 1
            endif
            allocate(filetime_in_sp(nt_fullmet))
            if(iw.eq.1)allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
          endif
          MR_windfiles_nt_fullmet(iw)=nt_fullmet

          ! get id for time
          nSTAT = nf90_inq_varid(ncid,Met_dim_names(1),time_var_id)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: inq_varid:',"time",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: inq_varid:',"time",nf90_strerror(nSTAT)
            stop 1
          endif
          ! time is an interger*4
          nSTAT = nf90_get_var(ncid,time_var_id,filetime_in_sp)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: get_var:',"time",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: get_var:',"time",nf90_strerror(nSTAT)
            stop 1
          endif
          filetime_in_int = nint(filetime_in_sp(1))

          nSTAT = nf90_inq_varid(ncid,"reftime",reftime_var_id)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: inq_varid:',"reftime",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: inq_varid:',"reftime",nf90_strerror(nSTAT)
            stop 1
          endif

          nSTAT = nf90_get_var(ncid,reftime_var_id,tstring)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: get_var:',"reftime",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: get_var:',"reftime",nf90_strerror(nSTAT)
            stop 1
          endif

          read(tstring,131)itstart_year,itstart_month,itstart_day, &
                            itstart_hour,itstart_min
          itstart_sec = 0
          filestart_hour = real(itstart_hour+filetime_in_int,kind=sp) + &
                           real(itstart_min,kind=sp)/60.0_sp      + &
                           real(itstart_sec,kind=sp)/3600.0_sp
 131      format(i4,1x,i2,1x,i2,1x,i2,1x,i2)

          nSTAT = nf90_close(ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: Could not close file',nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: Could not close file:',nf90_strerror(nSTAT)
            stop 1
          endif

          MR_windfile_starthour(iw) =  real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                         itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)

          do iwstep = 1,nt_fullmet
            MR_windfile_stephour(iw,iwstep) = MR_windfile_stephour(iw,1) + filetime_in_sp(iwstep)
          enddo
        enddo
      elseif(MR_iwindformat.eq.50)then
        ! Branch for WRF files
        ! Loop through all the windfiles
        do iw = 1,MR_iwindfiles
          nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: nf90_open to read header:', &
                           nf90_strerror(nSTAT)
            write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
            write(MR_global_error,*)'Exiting'
            stop 1
          endif
          if(iw.eq.1)then
            nSTAT = nf90_inq_varid(ncid,"Times",time_var_id)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_varid:',"Times",nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: inq_varid:',"Times",nf90_strerror(nSTAT)
              stop 1
            endif
            nSTAT = nf90_inquire_variable(ncid, time_var_id, invar, &
                xtype = var_xtype)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: inq_variable:',"Times",nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: inq_variable:',"Times",nf90_strerror(nSTAT)
              stop 1
            endif
            allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
            allocate(filetime_in_sp(nt_fullmet))
            If(nt_fullmet.gt.1)then
              write(MR_global_error,*)"MR ERROR: Currently WRF files are expected to only have one"
              write(MR_global_error,*)"       timestep/file"
              stop 1
            endif
            filetime_in_sp = 0.0_sp
          endif
          nSTAT = nf90_get_var(ncid,time_var_id,Timestr_WRF,&
                         start = (/1,1/),       &
                         count = (/19,1/))
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: get_var:',"Times",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: get_var:',"Times",nf90_strerror(nSTAT)
            stop 1
          endif

          nSTAT = nf90_close(ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: Could not close file',nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: Could not close file:',nf90_strerror(nSTAT)
            stop 1
          endif

          read(Timestr_WRF,121)itstart_year,itstart_month,itstart_day, &
                            itstart_hour,itstart_min,itstart_sec
          filestart_hour = real(itstart_hour,kind=sp) + &
                           real(itstart_min,kind=sp)/60.0_sp      + &
                           real(itstart_sec,kind=sp)/3600.0_sp
 121      format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2,1x)
          MR_windfiles_nt_fullmet(iw)=nt_fullmet
          MR_windfile_starthour(iw) = real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                       itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
          do iwstep = 1,nt_fullmet
            MR_windfile_stephour(iw,iwstep) = MR_windfile_stephour(iw,1) + filetime_in_sp(iwstep)
          enddo
        enddo

      else
        ! For all other formats, try to read the GRIB_orgReferenceTime string
        ! Loop through all the windfiles
        do iw = 1,MR_iwindfiles

          ! Each wind file needs a ref-time which in almost all cases is given
          ! in the 'units' attribute of the time variable
          write(MR_global_info,*)iw,trim(ADJUSTL(MR_windfiles(iw)))
          nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: nf90_open to read header:', &
                           nf90_strerror(nSTAT)
            write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
            write(MR_global_error,*)'Exiting'
            stop 1
          endif
          ! Find the id of the time dimension
          nSTAT = nf90_inq_dimid(ncid,trim(ADJUSTL(Met_dim_names(1))),t_dim_id)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
            write(MR_global_error,*)"    Could not find dimension: ",Met_dim_names(1)
            write(MR_global_log  ,*)'MR ERROR: inq_dimid time: ',nf90_strerror(nSTAT)
            stop 1
          endif
          if(iw.eq.1)then
            ! Get length of time dimension and allocate MR_windfile_stephour
            nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,len=nt_fullmet)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: Inquire_Dimension time: ', &
                                 nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension time: ', &
                                 nf90_strerror(nSTAT)
              stop 1
            endif
            write(MR_global_info,*)"  Assuming all NWP files have the same number of steps."
            write(MR_global_info,*)"   Allocating time arrays for ",MR_iwindfiles,"files"
            write(MR_global_info,*)"                              ",nt_fullmet,"step(s) each"
            allocate(MR_windfile_stephour(MR_iwindfiles,nt_fullmet))
          endif

          ! get variable id for time
          nSTAT = nf90_inq_varid(ncid,trim(ADJUSTL(Met_dim_names(1))),time_var_id)
          if(nSTAT.ne.NF90_NOERR) then
            write(MR_global_error,*)'MR ERROR: inq_varid:',"time",nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: inq_varid:',"time",nf90_strerror(nSTAT)
            stop 1
          endif
          ! We need the reftime for this file, check time variable for 'units'
          nSTAT = nf90_Inquire_Attribute(ncid, time_var_id,&
                                         "units",xtype, length, attnum)
          if(nSTAT.eq.0)then
            TimeHasUnitsAttr = .true.
            nSTAT = nf90_get_att(ncid, time_var_id,"units",tstring2)
            if(nSTAT.ne.NF90_NOERR) then
              write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
              stop 1
            endif
          else
            ! Try GRIB_orgReferenceTime
            nSTAT = nf90_Inquire_Attribute(ncid, time_var_id,&
                                           "GRIB_orgReferenceTime",xtype, length, attnum)
            if(nSTAT.ne.NF90_NOERR) then
              write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
              stop 1
            endif
            if(nSTAT.eq.0)then
              nSTAT = nf90_get_att(ncid, time_var_id,"GRIB_orgReferenceTime",tstring2)
              if(nSTAT.ne.NF90_NOERR) then
                write(MR_global_error,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_att:',"time",nf90_strerror(nSTAT)
                stop 1
              endif
              TimeHasUnitsAttr = .true.
            else
              TimeHasUnitsAttr = .false.
            endif
          endif

          if(TimeHasUnitsAttr)then
            do i=1,26
              ! try to parse
              !  time:units = "Hour since 2016-01-11T00:00:00Z" ;
              !  time:units = "days since 0001-01-01 00:00:00" ;
              if(tstring2(i:i+5).eq.'since ')then
                ii = i+6
                read(tstring2(ii:31),103)itstart_year,itstart_month,itstart_day, &
                                  itstart_hour,itstart_min,itstart_sec
                write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                           itstart_hour,itstart_min,itstart_sec
                filestart_hour = real(itstart_hour,kind=sp) + &
                                 real(itstart_min,kind=sp)/60.0_sp      + &
                                 real(itstart_sec,kind=sp)/3600.0_sp
                exit
              elseif(i.eq.26)then
                ! If we got to the end of the string without finding 'since',
                ! this may be an old MERRA file
                ii = 1
                read(tstring2(ii:31),103)itstart_year,itstart_month,itstart_day,&
                                  itstart_hour,itstart_min,itstart_sec
                write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                           itstart_hour,itstart_min,itstart_sec
                filestart_hour = real(itstart_hour,kind=sp) + &
                                 real(itstart_min,kind=sp)/60.0_sp      + &
                                 real(itstart_sec,kind=sp)/3600.0_sp
                exit
              endif
            enddo
          else
            ! Time variable does not have units attribute
            ! Try variable 'reftime'
            nSTAT = nf90_inq_varid(ncid,'reftime',reftime_var_id)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)"MR ERROR:  Could not read time:units or reftime"
              write(MR_global_error,*)"        Windfile start time is not defined."
              stop 1
            endif
            nSTAT = nf90_get_var(ncid,reftime_var_id,tstring2)
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)"MR ERROR:  Could not read reftime"
              write(MR_global_error,*)"        Windfile start time is not defined."
              stop 1
            endif
            do i=1,30
              if(tstring2(i:i+1).eq.'20'.or.tstring2(i:i+1).eq.'19')then
                write(MR_global_info,*)"Found reference time: ",tstring2(i:31)
                read(tstring2(i:31),103)itstart_year,itstart_month,itstart_day, &
                                  itstart_hour,itstart_min,itstart_sec
                write(MR_global_info,2100)"Ref time = ",itstart_year,itstart_month,itstart_day, &
                                           itstart_hour,itstart_min,itstart_sec
                filestart_hour = real(itstart_hour,kind=sp) + &
                                 real(itstart_min,kind=sp)/60.0_sp      + &
                                 real(itstart_sec,kind=sp)/3600.0_sp
                exit
              endif
            enddo
          endif
2100      format(20x,a11,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)
 103      format(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)

          ! Assume we now have the parsed reftime

          ! Now get time data
          ! Check if we need to read into an int, float or a double
          !nSTAT = nf90_inquire_variable(ncid, time_var_id, invar, xtype = var_xtype)
          nSTAT = nf90_inquire_variable(ncid, time_var_id, name = invar, xtype = var_xtype)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
            stop 1
          endif
          if(var_xtype.eq.NF90_FLOAT)then
            allocate(dum1d_sp(nt_fullmet))
            nSTAT = nf90_get_var(ncid,time_var_id,dum1d_sp, &
                   start = (/1/),count = (/nt_fullmet/))
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              stop 1
            endif
            ! copy to local variable
            MR_windfile_stephour(iw,1:nt_fullmet) = dum1d_sp(1:nt_fullmet)* &
                                                         Met_dim_fac(1)
            deallocate(dum1d_sp)
          elseif(var_xtype.eq.NF90_DOUBLE)then
            allocate(dum1d_dp(nt_fullmet))
            nSTAT = nf90_get_var(ncid,time_var_id,dum1d_dp, &
                   start = (/1/),count = (/nt_fullmet/))
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              stop 1
            endif
            ! copy to local variable
            MR_windfile_stephour(iw,1:nt_fullmet) = real(dum1d_dp(1:nt_fullmet),kind=4)* &
                                                         Met_dim_fac(1)
            deallocate(dum1d_dp)
          elseif(var_xtype.eq.NF90_INT)then
            allocate(dum1d_int4(nt_fullmet))
            nSTAT = nf90_get_var(ncid,time_var_id,dum1d_int4, &
                   start = (/1/),count = (/nt_fullmet/))
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              write(MR_global_log  ,*)'MR ERROR: get_var ',Met_dim_names(1),nf90_strerror(nSTAT)
              stop 1
            endif
            ! copy to local variable
            MR_windfile_stephour(iw,1:nt_fullmet) = real(dum1d_int4(1:nt_fullmet),kind=4)* &
                                                         Met_dim_fac(1)
            deallocate(dum1d_int4)
          else
            write(MR_global_error,*)"MR ERROR: Unexpected time variable type ",Met_dim_names(i)
            stop 1
          endif

          nSTAT = nf90_close(ncid)
          if(nSTAT.ne.NF90_NOERR)then
            write(MR_global_error,*)'MR ERROR: Could not close file',nf90_strerror(nSTAT)
            write(MR_global_log  ,*)'MR ERROR: Could not close file:',nf90_strerror(nSTAT)
            stop 1
          endif

          MR_windfiles_nt_fullmet(iw) = nt_fullmet
          MR_windfile_starthour(iw) =  real(HS_hours_since_baseyear(itstart_year,itstart_month, &
                                         itstart_day,real(filestart_hour,kind=8),MR_BaseYear,MR_useLeap),kind=4)
        enddo
      endif
      ! Finished setting up the start time of each wind file in HoursSince : MR_windfile_starthour(iw)
      !  and the forecast (offset from start of file) for each step        : MR_windfile_stephour(iw,iwstep)

      if (MR_iwindformat.ne.25.and.MR_iwindformat.ne.27)then
        write(MR_global_info,*)"File, step, Ref, Offset, HoursSince"
        do iw = 1,MR_iwindfiles
          do iws = 1,nt_fullmet
            write(MR_global_info,*)iw,iws,real(MR_windfile_starthour(iw),kind=4),&
                             real(MR_windfile_stephour(iw,iws),kind=4),&
                             real(MR_windfile_starthour(iw)+MR_windfile_stephour(iw,iws),kind=4)
          enddo
        enddo
      endif

      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_Met_Times_netcdf
!##############################################################################

!##############################################################################
!
!     MR_Set_Met_Dims_Template_netcdf
!
!     Called from MR_Set_MetComp_Grids_netcdf
!
!     Sets Met grid for Template windfiles
!
!##############################################################################


      subroutine MR_Set_Met_Dims_Template_netcdf

      use MetReader
      use netcdf

      use projection

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision
      real(kind=sp), parameter :: tol = 1.0e-7_sp

      integer :: iw,i

      integer :: nSTAT
      integer :: ncid

      !integer            :: var_xtype
      integer            :: xtype
      integer            :: length
      integer            :: attnum
      character(len=40)  :: invar
      character(len = nf90_max_name) :: name_dum
      integer :: t_dim_id
      !integer :: x_dim_id,y_dim_id
      integer :: var_id
      !real(kind=dp),dimension(:), allocatable :: dum1d_dp
      !real(kind=sp),dimension(:), allocatable :: dum1d_sp
      real(kind=sp):: dum_sp
      logical :: TimeHasFillVAttr

      write(MR_global_production,*)"--------------------------------------------------------------------------------"
      write(MR_global_production,*)"----------                MR_Set_Met_Dims_Template_netcdf             ----------"
      write(MR_global_production,*)"--------------------------------------------------------------------------------"

      ! To set up the grid, we assume that the grid is the same for all
      ! windfiles.  There is no checking if this is actually the case.
      ! Just read the first windfile.
      iw = 1
      write(MR_global_info,*)"Opening ",iw,trim(ADJUSTL(MR_windfiles(iw)))
      nSTAT = nf90_open(trim(ADJUSTL(MR_windfiles(iw))),NF90_NOWRITE,ncid)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: nf90_open to read header:', nf90_strerror(nSTAT)
        write(MR_global_error,*)'Could not open ',trim(ADJUSTL(MR_windfiles(iw)))
        write(MR_global_error,*)'Exiting'
        stop 1
      else
        write(MR_global_info,*)"Opened ",trim(ADJUSTL(MR_windfiles(iw)))
      endif

      ! Get dim ids, sizes, and associated dimension variable for dims:
      !  1 = time
      !  2 = pressure used for state variables
      !  3 = y or lat
      !  4 = x or lon
      !  5 = pressure used for Vz
      !  6 = pressure uesed for RH or SH
      !  7 = height above ground
      !  8 = depth below surface
      !  9 = extra pressure dimension

      ! Time
      !  This will be repeated in MR_Set_Met_Times_netcdf where the time values
      !  are read, but for now, we just want the dimension size
      i = 1
      if(.not.Met_dim_IsAvailable(i))then
        write(MR_global_error,*)"MR ERROR: TIME dimension is required and not listed"
        write(MR_global_error,*)"          in template windfile specification file."
        stop 1
      endif
      nSTAT = nf90_inq_dimid(ncid,Met_dim_names(i),t_dim_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_dimid ',Met_dim_names(i),nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_dimid ',Met_dim_names(i),nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_Inquire_Dimension(ncid,t_dim_id,name=name_dum,len=nt_fullmet)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Inquire_Dimension ',Met_dim_names(i),nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Inquire_Dimension ',Met_dim_names(i),nf90_strerror(nSTAT)
        stop 1
      endif

      ! We need to check if this is a regular grid
      !IsRegular_MetGrid = .true.
      !do i = 1,nx_fullmet-1
      !  if(abs(MR_dx_met(i+1)-MR_dx_met(i)).gt.tol*MR_dx_met(i))then
      !    IsRegular_MetGrid = .false.
      !  endif
      !enddo
      !do i = 1,ny_fullmet-1
      !  if(abs(MR_dy_met(i+1)-MR_dy_met(i)).gt.tol*MR_dy_met(i))then
      !    IsRegular_MetGrid = .false.
      !  endif
      !enddo

      ! Need to get fill_value
      ! Try to get it from geopotential height variable
      nSTAT = nf90_inq_varid(ncid,Met_var_NC_names(1),var_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_varid: ',invar,nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_varid: ',invar,nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_Inquire_Attribute(ncid, var_id,&
                                     "_FillValue",xtype, length, attnum)
      if(nSTAT.eq.0)then
        TimeHasFillVAttr = .true.
      else
        TimeHasFillVAttr = .false.
      endif
      if(TimeHasFillVAttr)then
        nSTAT = nf90_get_att(ncid, var_id,"_FillValue",dum_sp)
        fill_value_sp(MR_iwindformat) = dum_sp
        write(MR_global_info,*)"    Found fill value",fill_value_sp(MR_iwindformat)
      else
        fill_value_sp(MR_iwindformat) = -9999.0_sp
      endif

      nSTAT = nf90_close(ncid)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: Could not close file',nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: Could not close file:',nf90_strerror(nSTAT)
        stop 1
      endif

      end subroutine MR_Set_Met_Dims_Template_netcdf

!##############################################################################
!
!     MR_Read_MetP_Variable_netcdf
!
!     Called from Read_HGT_arrays and once from Read_3d_MetP_Variable.
!
!     Sets MR_dum3d_metP, MR_dum2d_met, or MR_dum2d_met_int as appropriate
!
!##############################################################################

      subroutine MR_Read_MetP_Variable_netcdf(ivar,istep)

      use MetReader
      use netcdf

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer,intent(in) :: ivar
      integer,intent(in) :: istep

      integer :: iw,iwstep
      integer :: np_met_loc
      character(len=130) :: infile
      character(len=71)  :: invar

      integer :: ncid       = 0
      integer :: nSTAT      = 0
      integer :: in_var_id  = 0
      integer :: in_var_id1 = 0
      integer :: in_var_id2 = 0

      real(kind=sp) :: del_H,del_P,dpdz

      integer :: i,j,k,ii,jj,kk,kkk,itmp
      integer :: ict, ileft(2),iright(2)   !if wrapgrid=.true. ict=2 and left & iright have 2 values, otherwise 1
      integer :: iistart(2),iicount(2)     !if (wrapgrid), iistart(1)=istart, iistart(2)=1

      integer :: Dimension_of_Variable
      logical :: IsCatagorical

      integer :: var_xtype
      integer :: NC_version

      real(kind=sp),dimension(:,:,:,:),allocatable :: dum3d_metP_aux
      real(kind=sp) :: theta,cofac

      real(kind=sp) :: Z_top, T_top
      real(kind=sp) :: pp
      integer       :: idx

      !write(MR_global_production,*)"--------------------------------------------------------------------------------"
      !write(MR_global_production,*)"----------                MR_Read_MetP_Variable_netcdf                ----------"
      !write(MR_global_production,*)"--------------------------------------------------------------------------------"

      if(.not.Met_var_IsAvailable(ivar))then
        write(MR_global_error,*)"MR ERROR:  Variable not available for this windfile"
        write(MR_global_error,*)"             ivar = ",ivar
        write(MR_global_error,*)"            vname = ",Met_var_NC_names(ivar)
        write(MR_global_error,*)"             iwf  = ",MR_iwindformat
        stop 1
      endif

      iw     = MR_MetStep_findex(istep)
      iwstep = MR_MetStep_tindex(istep)

      if(Met_var_NC_names(ivar).eq."")then
        write(MR_global_info,*)"Variable ",ivar," not available for MR_iwindformat = ",&
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
      if(MR_iwindformat.eq.24)then
          ! NASA MERRA has 3d precipitation
        if(ivar.eq.44) Dimension_of_Variable = 3 ! Precipitation rate large-scale (liquid)
        if(ivar.eq.45) Dimension_of_Variable = 3 ! Precipitation rate convective (liquid)
      else
          ! All other met files use surface precip
        if(ivar.eq.44) Dimension_of_Variable = 2 ! Precipitation rate large-scale (liquid)
        if(ivar.eq.45) Dimension_of_Variable = 2 ! Precipitation rate convective (liquid)
      endif
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

      if(MR_iwind.eq.5)then
        ! For iwind=5, variables are given in different files.
        ! We need to explicitly create the file name, given the input directory.
        if(MR_iwindformat.eq.25)then
          ! Get correct file
          if(ivar.eq.1)then
            write(infile,115)trim(adjustl(MR_MetStep_File(istep))), &
                             "hgt.",MR_iwind5_year(istep),".nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.2)then
            write(infile,116)trim(adjustl(MR_MetStep_File(istep))), &
                             "uwnd.",MR_iwind5_year(istep),".nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.3)then
            write(infile,116)trim(adjustl(MR_MetStep_File(istep))), &
                             "vwnd.",MR_iwind5_year(istep),".nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.4)then
            write(infile,117)trim(adjustl(MR_MetStep_File(istep))), &
                             "omega.",MR_iwind5_year(istep),".nc"
            np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          elseif(ivar.eq.5)then
            write(infile,115)trim(adjustl(MR_MetStep_File(istep))), &
                             "air.",MR_iwind5_year(istep),".nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.30)then
            write(infile,116)trim(adjustl(MR_MetStep_File(istep))), &
                             "rhum.",MR_iwind5_year(istep),".nc"
            np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          elseif(ivar.eq.31)then
            write(infile,116)trim(adjustl(MR_MetStep_File(istep))), &
                             "shum.",MR_iwind5_year(istep),".nc"
            np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          else
            write(MR_global_info,*)"Requested variable not available."
            stop 1
          endif
          infile = trim(adjustl(infile))

 115      format(a50,a4,i4,a3)
 116      format(a50,a5,i4,a3)
 117      format(a50,a6,i4,a3)
! 118      format(a50,a16,i4,a3)
! 119      format(a50,a15,i4,a3)
        elseif(MR_iwindformat.eq.26)then
         ! JRA-55 reanalysis 1.25 degree files 
          ! Get correct file
          if(ivar.eq.1)then
            write(infile,125)trim(adjustl(MR_MetStep_File(istep))), &
                             "anl_p125.007_hgt.",MR_iwind5_year(istep),"060100_2018063018.nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.2)then
            write(infile,126)trim(adjustl(MR_MetStep_File(istep))), &
                             "anl_p125.033_ugrd.",MR_iwind5_year(istep),"060100_2018063018.nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.3)then
            write(infile,126)trim(adjustl(MR_MetStep_File(istep))), &
                             "anl_p125.034_vgrd.",MR_iwind5_year(istep),"060100_2018063018.nc"
            np_met_loc = np_fullmet
          elseif(ivar.eq.4)then
            write(infile,126)trim(adjustl(MR_MetStep_File(istep))), &
                             "anl_p125.039_vvel.",MR_iwind5_year(istep),"060100_2018063018.nc"
            np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          elseif(ivar.eq.5)then
            write(infile,125)trim(adjustl(MR_MetStep_File(istep))), &
                             "anl_p125.011_tmp.",MR_iwind5_year(istep),"060100_2018063018.nc"
            np_met_loc = np_fullmet
          !elseif(ivar.eq.31)then
          !  write(infile,126)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "shum.",MR_iwind5_year(istep),".nc"
          !  np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          else
            write(MR_global_info,*)"Requested variable not available."
            stop 1
          endif
          infile = trim(adjustl(infile))

 125      format(a50,a17,i4,a20)
 126      format(a50,a18,i4,a20)
! 128      format(a50,a16,i4,a3)
! 129      format(a50,a15,i4,a3)

        elseif(MR_iwindformat.eq.27)then
          ! Get correct file
          !if(ivar.eq.1)then
          !  write(infile,135)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_HGT_pres.nc"
          !  np_met_loc = np_fullmet
          !elseif(ivar.eq.2)then
          !  write(infile,136)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_UGRD_pres.nc"
          !  np_met_loc = np_fullmet
          !elseif(ivar.eq.3)then
          !  write(infile,136)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_VGRD_pres.nc"
          !  np_met_loc = np_fullmet
          !elseif(ivar.eq.4)then
          !  write(infile,136)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_VVEL_pres.nc"
          !  !np_met_loc = np_fullmet_Vz
          !  np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          !elseif(ivar.eq.5)then
          !  write(infile,135)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_TMP_pres.nc"
          !  np_met_loc = np_fullmet
          !elseif(ivar.eq.10)then
          !  write(infile,138)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "sflxgrbfg_mean_",MR_iwind5_year(istep), &
          !                   "_HPBL_sfc.nc"
          !elseif(ivar.eq.22)then
          !  write(infile,140)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "sflxgrbfg_mean_",MR_iwind5_year(istep), &
          !                   "_TMP_low-cldtop.nc"
          !elseif(ivar.eq.23)then
          !  write(infile,141)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "sflxgrbfg_mean_",MR_iwind5_year(istep), &
          !                   "_TCDC_low-cldlay.nc"
          !elseif(ivar.eq.30)then
          !  write(infile,137)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "pgrbanl_mean_",MR_iwind5_year(istep), &
          !                   "_RH_pres.nc"
          !  !np_met_loc = np_fullmet_RH
          !  np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          !elseif(ivar.eq.44)then
          !  write(infile,139)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "sflxgrbfg_mean_",MR_iwind5_year(istep), &
          !                   "_PRATE_sfc.nc"
          !elseif(ivar.eq.45)then
          !  write(infile,139)trim(adjustl(MR_MetStep_File(istep))), &
          !                   "sflxgrbfg_mean_",MR_iwind5_year(istep), &
          !                   "_CPRAT_sfc.nc"
          !else
          !  write(MR_global_info,*)"Requested variable not available."
          !  stop 1
          !endif
          infile = trim(adjustl(infile))

! 135      format(a50,a13,i4,a12)
! 136      format(a50,a13,i4,a13)
! 137      format(a50,a13,i4,a11)
! 138      format(a50,a15,i4,a12)
! 139      format(a50,a15,i4,a13)
! 140      format(a50,a15,i4,a18)
! 141      format(a50,a15,i4,a19)
        elseif(MR_iwindformat.eq.29)then
         ! ECMWF ERA5
        elseif(MR_iwindformat.eq.30)then
         ! ECMWF ERA-20C
        else
          write(MR_global_info,*)"MR ERROR: MR_iwind=5, but MR_iwindformat is not a recognized type."
          write(MR_global_info,*)"            MR_iwindformat=25  NCEP 20th Cent. Reanalysis"
          write(MR_global_info,*)"                           27  NOAA-CIRES"
          write(MR_global_info,*)"                           29  ERA5"
          write(MR_global_info,*)"                           29  ERA-20C"
          write(MR_global_info,*)"           MR_iwindformat provided = ",MR_iwindformat
          stop 1
        endif
      else  ! cases for MR_iwind not equal to 5
        np_met_loc = nlevs_fullmet(Met_var_zdim_idx(ivar))
          ! Files are listed directly, not through directories (as in MR_iwindformat=25,27)
        infile = trim(adjustl(MR_MetStep_File(istep)))
      endif

      invar = Met_var_NC_names(ivar)

      write(MR_global_info,*)istep,ivar,"Reading ",trim(adjustl(invar))," from file : ",&
                trim(adjustl(infile))!,nx_submet,ny_submet,np_met_loc
      nSTAT = nf90_open(trim(adjustl(infile)),NF90_NOWRITE,ncid)

      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR open file:',infile,nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR open file:',infile,nf90_strerror(nSTAT)
        write(MR_global_error,*)trim(adjustl(infile)),NF90_NOWRITE,ncid,nSTAT
        write(MR_global_error,*)'Exiting'
        stop 1
      endif

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

      nSTAT = nf90_inq_varid(ncid,invar,in_var_id)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_varid: ',invar,nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_varid: ',invar,nf90_strerror(nSTAT)
        stop 1
      endif
      nSTAT = nf90_inquire_variable(ncid, in_var_id, invar, &
                xtype = var_xtype)
      if(nSTAT.ne.NF90_NOERR)then
        write(MR_global_error,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
        write(MR_global_log  ,*)'MR ERROR: inq_variable: ',invar,nf90_strerror(nSTAT)
        stop 1
      endif

      ! Test for version 3 vs 4 NCEP files
      if(var_xtype.eq.NF90_FLOAT)then
        NC_version = 4
      elseif(var_xtype.eq.NF90_SHORT)then
        NC_version = 3
      endif

      if(Dimension_of_Variable.eq.3)then
        MR_dum3d_metP = 0.0_sp
        If(MR_iwindformat.ne.50)then
          allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
        else
            ! For MR_iwindformat = 50 (WRF), we need an extra point in p
            ! Allocate auxillary array
          if(ivar.eq.1)then
                ! Geopotential lives on z-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc+1,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc+1,1))
          elseif(ivar.eq.2)then
                ! U wind lives on x-staggered grid
            allocate(temp3d_sp(nx_submet+1,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet+1,ny_submet,np_met_loc,1))
          elseif(ivar.eq.3)then
                ! V wind lives on y-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet+1,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet+1,np_met_loc,1))
          elseif(ivar.eq.4)then
                ! W wind lives on z-staggered grid
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc+1,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc+1,1))
          elseif(ivar.eq.5)then
                ! Temp lives on non-staggered grid, but we need pres. base and pert.
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc,1))
          elseif(ivar.eq.6)then
                ! pressure lives on non-staggered grid, but we need base and pert.
            allocate(temp3d_sp(nx_submet,ny_submet,np_met_loc,1))
            allocate(dum3d_metP_aux(nx_submet,ny_submet,np_met_loc,1))
          endif
        endif ! MR_iwindformat.ne.50

        do i=1,ict        !read subgrid at current time step
          if(MR_iwind.eq.5.and.MR_iwindformat.eq.25)then
            ! NCEP reanalysis files are now NCv4 (stored as float), but the
            ! older version, NCv3 (stored as short) might still be around.
            if(i.eq.1)allocate(temp3d_short(nx_submet,ny_submet,np_met_loc,1))
            if(var_xtype.eq.NF90_FLOAT)then
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT),iicount(i),ny_submet,np_met_loc
                write(MR_global_log  ,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT),iicount(i),ny_submet,np_met_loc
                stop 1
              endif
            elseif(var_xtype.eq.NF90_SHORT)then
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_short(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
                stop 1
              endif
            endif
          elseif(MR_iwind.eq.5.and.MR_iwindformat.eq.26)then
            !  JRA-55
          elseif(MR_iwind.eq.5.and.MR_iwindformat.eq.27)then
            !  NOAA-CIRES
          elseif(MR_iwind.eq.5.and.MR_iwindformat.eq.29)then
            !  ECMWF ERA5
          elseif(MR_iwind.eq.5.and.MR_iwindformat.eq.30)then
            !  ECMWF ERA-20C
          elseif(MR_iwindformat.eq.50)then

            ! Now read the data and convert if necessary
            if(ivar.eq.1)then
                ! Geopotential
                ! Get PHB
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: PHB',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: PHB',nf90_strerror(nSTAT)
                stop 1
              endif
              write(MR_global_info,*)istep,"Reading ","PH"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"PH",in_var_id2)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_var: PH',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_var: PH',nf90_strerror(nSTAT)
                stop 1
              endif
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: PH',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: PH',nf90_strerror(nSTAT)
                stop 1
              endif
              temp3d_sp(:,:,:,:) = temp3d_sp(:,:,:,:) + dum3d_metP_aux(:,:,:,:)
              do kk=1,np_met_loc
                MR_dum3d_metP(:,:,kk) = 0.5_sp*(temp3d_sp(:,:,kk  ,1) + &
                                          temp3d_sp(:,:,kk+1,1))
              enddo

            elseif(ivar.eq.2)then
                ! U wind lives on x-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i)+1,:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i)+1,ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then 
                write(MR_global_error,*)'MR ERROR: get_var: U',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: U',nf90_strerror(nSTAT)
                stop 1
              endif
              do ii=1,iicount(i)
                MR_dum3d_metP(ii,:,:) = 0.5_sp*(temp3d_sp(ii  ,:,:,1) + &
                                          temp3d_sp(ii+1,:,:,1))
              enddo
            elseif(ivar.eq.3)then
                ! V wind lives on y-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet+1,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: V',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: V',nf90_strerror(nSTAT)
                stop 1
              endif
              do jj=1,ny_submet
                MR_dum3d_metP(:,jj,:) = 0.5_sp*(temp3d_sp(:,jj,:,1) + &
                                          temp3d_sp(:,jj,:,1))
              enddo
            elseif(ivar.eq.4)then
                ! W wind lives on z-staggered grid
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc+1,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: W',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: W',nf90_strerror(nSTAT)
                stop 1
              endif
              do kk=1,np_met_loc
                MR_dum3d_metP(:,:,kk) = 0.5_sp*(temp3d_sp(:,:,kk  ,1) + &
                                          temp3d_sp(:,:,kk+1,1))
              enddo
            elseif(ivar.eq.5)then
                ! Temperature is actually stored as potential temperature
                ! perturbation: we convert to potential temperature via
                !  theta = (pot.temp.pert + 300.0)
                ! and convert to thermodynamic temperature via
                !  Temp = theta*(p/p_0)^kappa
                !    where p_0   = 1000.0mb (or 1.0e5 Pa)
                !    and   kappa = R/c_p = 0.2854 (for dry air)
                ! Temp lives on non-staggered grid, but we need pres. base and pert.
                ! First get PB (base pressure)
              write(MR_global_info,*)istep,"Reading ","PB"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"PB",in_var_id1)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_var: PB',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_var: PB',nf90_strerror(nSTAT)
                stop 1
              endif
              nSTAT = nf90_get_var(ncid,in_var_id1,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: PB',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: PB',nf90_strerror(nSTAT)
                stop 1
              endif
                ! Now get P (pertubation pressure)
              write(MR_global_info,*)istep,"Reading ","P"," from file : ",trim(adjustl(infile))
              nSTAT = nf90_inq_varid(ncid,"P",in_var_id2)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_var: P',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_var: P',nf90_strerror(nSTAT)
                stop 1
              endif
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: P',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: P',nf90_strerror(nSTAT)
                stop 1
              endif
                ! Now get total pressure in Pa
              MR_dum3d_metP(:,:,:) = temp3d_sp(:,:,:,1) + dum3d_metP_aux(:,:,:,1)

                ! Read potential temperature perturbation
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: T',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: T',nf90_strerror(nSTAT)
                stop 1
              endif
              do ii=1,iicount(i)
                do jj=1,ny_submet
                  do kk=1,np_met_loc
                    theta = temp3d_sp(ii,jj,kk,1) + 300.0_sp
                    cofac = (MR_dum3d_metP(ii,jj,kk)*1.0e-5_sp)**(0.2854_sp)
                    MR_dum3d_metP(ii,jj,kk) = theta * cofac
                  enddo
                enddo
              enddo

            elseif(ivar.eq.6)then
                ! pressure
                ! Get PB (base pressure)
              nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: PB',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: PB',nf90_strerror(nSTAT)
                stop 1
              endif
                ! Get P (pertubation pressure)
              nSTAT = nf90_inq_varid(ncid,"P",in_var_id2)
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: inq_var: P',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: inq_var: P',nf90_strerror(nSTAT)
                stop 1
              endif
              nSTAT = nf90_get_var(ncid,in_var_id,dum3d_metP_aux(ileft(i):iright(i),:,:,:), &
                       start = (/iistart(i),jstart,1,iwstep/),       &
                       count = (/iicount(i),ny_submet,np_met_loc,1/))
              if(nSTAT.ne.NF90_NOERR)then
                write(MR_global_error,*)'MR ERROR: get_var: P',nf90_strerror(nSTAT)
                write(MR_global_log  ,*)'MR ERROR: get_var: P',nf90_strerror(nSTAT)
                stop 1
              endif
                ! Now get total pressure in Pa
              temp3d_sp(:,:,:,:) = temp3d_sp(:,:,:,:) + dum3d_metP_aux(:,:,:,:)
              MR_dum3d_metP(:,:,:) = temp3d_sp(:,:,:,1)
            else
            ! for any other 3d WRF variable, assume non-staggered grid

            endif

          else ! end of MR_iwind=5 and iwf=50 (WRF) sections

            ! for any other 3d variable (non-WRF, non-NCEP/2.5 reannalysis)
            nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                     start = (/iistart(i),jstart,1,iwstep/),       &
                     count = (/iicount(i),ny_submet,np_met_loc,1/))
            if(nSTAT.ne.NF90_NOERR)then
              write(MR_global_error,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
              write(MR_global_error,*)i
              write(MR_global_error,*)ileft(i),iright(i)
              write(MR_global_error,*)iistart(i),jstart,1,iwstep
              write(MR_global_error,*)iicount(i),ny_submet,np_met_loc,1
              write(MR_global_log  ,*)'MR ERROR: get_var: ',nf90_strerror(nSTAT)
              stop 1
            endif
          endif
        enddo
        !  At this point, we have temp3d_sp with the direct read from the file
        !  Next, we need to copy it to MR_dum3d_metP

        If(MR_iwindformat.ne.50)then
          do j=1,ny_submet
            itmp = ny_submet-j+1
            !reverse the j indices (since they increment from N to S)
            if(MR_iwindformat.eq.25)then
              if(var_xtype.eq.NF90_FLOAT)then
                  ! No scaling/offset needed
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc) = &
                                      real(temp3d_sp(:,itmp,:,1),kind=sp)
              elseif(var_xtype.eq.NF90_SHORT)then
                  ! NC_version = 3 :: need to scale and offset shorts to get float
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc) = &
                                      real(temp3d_short(1:nx_submet,itmp,1:np_met_loc,1),kind=sp) * &
                                    iwf25_scale_facs(ivar) + iwf25_offsets(ivar)
              endif
            else
              if(y_inverted)then
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,itmp,1:np_met_loc,1)
              else
                MR_dum3d_metP(1:nx_submet,j,1:np_met_loc)  = temp3d_sp(1:nx_submet,j,1:np_met_loc,1)
              endif
            endif
          enddo
          if(z_inverted)then ! reverse the vertical coordinate
            temp3d_sp(1:nx_submet,1:ny_submet,1:np_met_loc,1) = MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc)
            do i=1,np_met_loc
              itmp = np_met_loc-i+1
              MR_dum3d_metP(1:nx_submet,1:ny_submet,i) = temp3d_sp(1:nx_submet,1:ny_submet,itmp,1)
            enddo
          endif
        endif !MR_iwindformat.eq.50, MR_iwindformat.eq.25, else

        if(MR_iwind.eq.5.and.MR_iwindformat.eq.25)then
          deallocate(temp3d_short)
        endif
        if(MR_iwindformat.eq.50)then
          deallocate(dum3d_metP_aux)
        endif
        deallocate(temp3d_sp)

      elseif(Dimension_of_Variable.eq.2)then
        if(IsCatagorical)then
          allocate(temp2d_int(nx_submet,ny_submet,1))
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then
              ! No catagorical variables for MR_iwindformat = 25
            else
              nSTAT = nf90_get_var(ncid,in_var_id,temp2d_int(ileft(i):iright(i),:,:), &
                         start = (/iistart(i),jstart,iwstep/),       &
                         count = (/iicount(i),ny_submet,1/))
              do j=1,ny_submet
                itmp = ny_submet-j+1
                if(y_inverted)then
                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,itmp,1)
                else
                  MR_dum2d_met_int(1:nx_submet,j)  = temp2d_int(1:nx_submet,j,1)
                endif
              enddo
            endif
            if(nSTAT.ne.NF90_NOERR)then
               write(MR_global_error,*)'MR ERROR: get_var:Vx ',invar,nf90_strerror(nSTAT)
               write(MR_global_log  ,*)'MR ERROR: get_var:Vx ',invar,nf90_strerror(nSTAT)
               stop 1
             endif
          enddo
          deallocate(temp2d_int)
        else
          allocate(temp2d_sp(nx_submet,ny_submet,1))
          if(ivar.eq.11.or.ivar.eq.12)then
              ! Surface winds usually have a z coordinate as well
            allocate(temp3d_sp(nx_submet,ny_submet,1,1))
          endif
  
          do i=1,ict        !read subgrid at current time step
            if(MR_iwindformat.eq.25)then
              allocate(tmpsurf2d_short(192,94,1))
              if(NC_version.eq.4)then
                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(:,:,1), &
                         start = (/1,1,iwstep/),       &
                         count = (/192,94,1/))
                if(nSTAT.ne.NF90_NOERR)then
                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   stop 1
                endif
                write(MR_global_info,*)"Need to write interp_iwf25_grid for float"
                stop 1
                call MR_interp_iwf25_grid(nx_submet,ny_submet,tmpsurf2d_short,temp2d_sp,&
                                    iwf25_scale_facs(ivar),iwf25_offsets(ivar))
                MR_dum2d_met(1:nx_submet,:) = temp2d_sp(1:nx_submet,:,1)
              else
                nSTAT = nf90_get_var(ncid,in_var_id,tmpsurf2d_short(:,:,1), &
                         start = (/1,1,iwstep/),       &
                         count = (/192,94,1/))
                if(nSTAT.ne.NF90_NOERR)then
                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   stop 1
                endif
                call MR_interp_iwf25_grid(nx_submet,ny_submet,tmpsurf2d_short,temp2d_sp,&
                                    iwf25_scale_facs(ivar),iwf25_offsets(ivar))
                MR_dum2d_met(1:nx_submet,:) = temp2d_sp(1:nx_submet,:,1)
              endif
            else
              ! 2d variables for iwf .ne. 25
              if(ivar.eq.11.or.ivar.eq.12)then
                ! Surface velocities do have a z dimension
                nSTAT = nf90_get_var(ncid,in_var_id,temp3d_sp(ileft(i):iright(i),:,:,:), &
                         start = (/iistart(i),jstart,1,iwstep/),       &
                         count = (/iicount(i),ny_submet,1,1/))
                if(nSTAT.ne.NF90_NOERR)then
                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   stop 1
                endif
                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,itmp,1,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp3d_sp(1:nx_submet,j,1,1)
                  endif
                enddo
              else
                nSTAT = nf90_get_var(ncid,in_var_id,temp2d_sp(ileft(i):iright(i),:,:), &
                         start = (/iistart(i),jstart,iwstep/),       &
                         count = (/iicount(i),ny_submet,1/))
                if(nSTAT.ne.NF90_NOERR)then
                   write(MR_global_error,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   write(MR_global_log  ,*)'MR ERROR: get_var: ',invar,nf90_strerror(nSTAT)
                   stop 1
                endif
                do j=1,ny_submet
                  itmp = ny_submet-j+1
                  if(y_inverted)then
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,itmp,1)
                  else
                    MR_dum2d_met(1:nx_submet,j)  = temp2d_sp(1:nx_submet,j,1)
                  endif
                enddo
              endif
            endif
          enddo
          deallocate(temp2d_sp)
          if(ivar.eq.11.or.ivar.eq.12) deallocate(temp3d_sp)
          if(MR_iwindformat.eq.25) deallocate(tmpsurf2d_short)
        endif ! IsCatagorical
      endif ! Dimension_of_Variable.eq.2


      ! Quality control checks
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
        If(ivar.eq.2.or.ivar.eq.3.or.ivar.eq.4)then
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
          idx = Met_var_zdim_idx(ivar)
          pp = levs_fullmet_sp(idx,nlevs_fullmet(idx))/real(100.0,kind=sp)
          Z_top = MR_Z_US_StdAtm(pp)
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
        if(MR_iwindformat.ne.50)then
            ! For pressure vertical velocity, convert from Pa s to m/s by dividing
            ! by pressure gradient
          idx = Met_var_zdim_idx(ivar)
          do k=1,np_met_loc
            do i=1,nx_submet
              do j=1,ny_submet
                if(k.eq.1)then
                  ! Use one-sided gradients for bottom
                  del_p = levs_fullmet_sp(idx,2) - levs_fullmet_sp(idx,1)
                  !del_P = p_fullmet_Vz_sp(2)-p_fullmet_Vz_sp(1)
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
                del_H = del_H * 1000.0_sp ! convert to m
                if(abs(del_H).gt.MR_EPS_SMALL)then
                  dpdz  = del_P/del_H
                else
                  write(MR_global_error,*)'MR ERROR: failed to calculate dpdz'
                  write(MR_global_error,*)i,j,k,del_P,del_H
                  write(MR_global_error,*)MR_geoH_metP_last(i,j,:)
                  stop 1
                endif
                MR_dum3d_metP(i,j,k) = MR_dum3d_metP(i,j,k) / dpdz
              enddo ! j
            enddo ! i
          enddo ! k
        endif
      endif

      MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc) =  &
      MR_dum3d_metP(1:nx_submet,1:ny_submet,1:np_met_loc) * Met_var_conversion_factor(ivar)

      !write(MR_global_production,*)"--------------------------------------------------------------------------------"

      end subroutine MR_Read_MetP_Variable_netcdf

!##############################################################################
!
!    MR_interp_iwf25_grid
!
!##############################################################################

      subroutine MR_interp_iwf25_grid(imax,jmax,invar,outvar,scale_fac,offset)

      use MetReader

      implicit none

      integer, parameter :: sp        = 4 ! single precision
      integer, parameter :: dp        = 8 ! double precision

      integer         ,intent(in)  :: imax,jmax
      integer(kind=sp),intent(in)  :: invar(192,94,1)
      real(kind=sp)   ,intent(out) :: outvar(imax,jmax)
      real(kind=sp)   ,intent(in)  :: scale_fac
      real(kind=sp)   ,intent(in)  :: offset

      real(kind=sp)    :: a1,a2,a3,a4
      real(kind=sp)    :: v1,v2,v3,v4

      integer :: ilon,ilat

      do ilon = 1,imax
        do ilat = 1,jmax
          a1 = amap_iwf25(ilon,ilat,1)
          a2 = amap_iwf25(ilon,ilat,2)
          a3 = amap_iwf25(ilon,ilat,3)
          a4 = amap_iwf25(ilon,ilat,4)
          v1 = invar(imap_iwf25(ilon,ilat,1),imap_iwf25(ilon,ilat,3),1) &
                 * scale_fac + offset
          v2 = invar(imap_iwf25(ilon,ilat,2),imap_iwf25(ilon,ilat,3),1) &
                 * scale_fac + offset
          v3 = invar(imap_iwf25(ilon,ilat,2),imap_iwf25(ilon,ilat,4),1) &
                 * scale_fac + offset
          v4 = invar(imap_iwf25(ilon,ilat,1),imap_iwf25(ilon,ilat,4),1) &
                 * scale_fac + offset

          outvar(ilon,ilat) = a1*v1 + a2*v2 + a3*v3 + a4*v4
        enddo
      enddo

      end subroutine MR_interp_iwf25_grid


