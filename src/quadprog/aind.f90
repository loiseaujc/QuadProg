! Copyright (C) 1997-2010 Berwin A. Turlach <Berwin.Turlach@gmail.com>
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
! USA.
!
!
! this routine checks whether Aind has valid entries, i.e.,
!   1) 1<= Aind(1,i) <= n for i=1,...,q (number of constraints)
!   2) 1<= Aind(j,i) <= n for j=2,...,Aind(1,i)+1, i=1,...,q
!
! Aind is a m times q matrix constructed in Splus
subroutine aind(ind, m, q, n, ok)
   implicit none
   integer, intent(in) :: m, ind(m, *), q, n
   integer, intent(out) :: ok
   integer :: i, j

   ok = 0
   do i = 1, q
      if (ind(1, i) < 1 .or. ind(1, i) > n) return
      do j = 2, ind(1, i) + 1
         if (ind(j, i) < 1 .OR. ind(j, i) > n) return
      end do
   end do
   ok = 1
   return
end
