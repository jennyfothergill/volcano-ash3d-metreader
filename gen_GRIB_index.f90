      program gen_grib_index

      use MetReader

      use grib_api

      implicit none

      integer :: nargs
      character(len=130) :: lllinebuffer
      character(len=130)  :: grib_file

      nargs = iargc()
      if (nargs.ne.1) then
        write(0,*)"MR ERROR: no grib file given"
        stop 1
      else
        call getarg(1,lllinebuffer)
        read(lllinebuffer,*)grib_file
        call MR_Set_Gen_Index_GRIB(grib_file)
      endif

      end program
