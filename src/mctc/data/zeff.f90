! This file is part of mctc-lib.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file mctc/data/zeff.f90
!> Provides effective nuclear charges for all elements

module mctc_data_zeff
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_effective_charge


   interface get_effective_charge
      module procedure :: get_effective_charge_num
      module procedure :: get_effective_charge_sym
   end interface get_effective_charge


   integer, parameter :: max_elem = 118


  !> Effective nuclear charges from the def2-ECPs used for calculating the 
  !> reference polarizabilities for DFT-D4.
  real(wp), parameter :: effective_nuclear_charge(max_elem) = [ &
    &   1,                                                 2,  & ! H-He
    &   3, 4,                               5, 6, 7, 8, 9,10,  & ! Li-Ne
    &  11,12,                              13,14,15,16,17,18,  & ! Na-Ar
    &  19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,  & ! K-Kr
    &   9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,  & ! Rb-Xe
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Cs-Lu
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26, & ! Hf-Rn
    !  just copy & paste from above
    &   9,10,11,30,31,32,33,34,35,36,37,38,39,40,41,42,43,  & ! Fr-Lr
    &  12,13,14,15,16,17,18,19,20,21,22,23,24,25,26] ! Rf-Og

contains


!> Get effective nuclear charge for a given element symbol
elemental function get_effective_charge_sym(sym) result(zeff)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Effective nuclear charge
   real(wp) :: zeff

   zeff = get_effective_charge(to_number(sym))

end function get_effective_charge_sym


!> Get effective nuclear charge for a given atomic number
elemental function get_effective_charge_num(num) result(zeff)

   !> Atomic number
   integer, intent(in) :: num

   !> Effective nuclear charge
   real(wp) :: zeff

   if (num > 0 .and. num <= size(effective_nuclear_charge)) then
      zeff = effective_nuclear_charge(num)
   else
      zeff = 0.0_wp
   end if

end function get_effective_charge_num


end module mctc_data_zeff