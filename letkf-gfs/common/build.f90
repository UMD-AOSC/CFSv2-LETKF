program build
  !! this is a dummy module,
  !! used only so that fmkmf can find all of the dependencies required
  use sfmt
  use common
  use common_gfs
  use common_letkf
  use common_mpi
  use common_mpi_gfs
  use common_mtx
  use common_obs
  use common_obs_gfs
  use kdtree
  use letkf_params

  write (*,*) "THis is just a dummy program needed for compiler configuration, why are you running this!?!?"
end program build
