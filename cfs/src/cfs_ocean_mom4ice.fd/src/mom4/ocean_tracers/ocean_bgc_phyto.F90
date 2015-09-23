!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                   !!
!!                   GNU General Public License                      !!
!!                                                                   !!
!! This file is part of the Flexible Modeling System (FMS).          !!
!!                                                                   !!
!! FMS is free software; you can redistribute it and/or modify       !!
!! it and are expected to follow the terms of the GNU General Public !!
!! License as published by the Free Software Foundation.             !!
!!                                                                   !!
!! FMS is distributed in the hope that it will be useful,            !!
!! but WITHOUT ANY WARRANTY; without even the implied warranty of    !!
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     !!
!! GNU General Public License for more details.                      !!
!!                                                                   !!
!! You should have received a copy of the GNU General Public License !!
!! along with FMS; if not, write to:                                 !!
!!          Free Software Foundation, Inc.                           !!
!!          59 Temple Place, Suite 330                               !!
!!          Boston, MA  02111-1307  USA                              !!
!! or see:                                                           !!
!!          http://www.gnu.org/licenses/gpl.txt                      !!
!!                                                                   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!<OVERVIEW>
! Ocean Carbon Model Intercomparison Study II: Biotic module
! NOTE: This module is NOT available in the current release
!</OVERVIEW>


module  ocean_bgc_phyto_mod 
use ocean_types_mod,    only: ocean_thickness_type
logical, public :: do_ocean_bgc_phyto
public  :: ocean_bgc_phyto_bbc
public  :: ocean_bgc_phyto_end
public  :: ocean_bgc_phyto_init
public  :: ocean_bgc_phyto_sbc
public  :: ocean_bgc_phyto_source
public  :: ocean_bgc_phyto_start
public  :: ocean_bgc_phyto_tracer
contains

subroutine ocean_bgc_phyto_bbc 
return
end subroutine  ocean_bgc_phyto_bbc 

subroutine ocean_bgc_phyto_end(Thickness)  
type(ocean_thickness_type), intent(in) :: Thickness
return
end subroutine  ocean_bgc_phyto_end 

subroutine ocean_bgc_phyto_sbc(robert) 
real, intent(in)        :: robert 
return
end subroutine  ocean_bgc_phyto_sbc 

subroutine ocean_bgc_phyto_init
return
end subroutine ocean_bgc_phyto_init 

subroutine ocean_bgc_phyto_source(Thickness) 
type(ocean_thickness_type), intent(in) :: Thickness
return
end subroutine  ocean_bgc_phyto_source 

subroutine ocean_bgc_phyto_start  
return
end subroutine  ocean_bgc_phyto_start  

subroutine ocean_bgc_phyto_tracer 
return
end subroutine  ocean_bgc_phyto_tracer 

end module  ocean_bgc_phyto_mod 
