      program gen_grib2_index

      use MetReader
      use grib_api

      implicit none

      integer :: nargs
      character(len=130) :: lllinebuffer

      character(len=130)  :: grib2_file !='nam.t00z.hawaiinest.hiresf01.tm00.grib2'

      nargs = iargc()
      if (nargs.ne.1) then
        !write(MR_global_error,*)"MR ERROR: no grib2 file given"
        write(0,*)"MR ERROR: no grib2 file given"
        stop 1
      else
        call getarg(1,lllinebuffer)
        read(lllinebuffer,*)grib2_file
        call MR_Set_Gen_Index_GRIB(grib2_file)
      endif

      end program
