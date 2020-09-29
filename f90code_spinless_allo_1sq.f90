subroutine eigenvalz(n,mat,evs)
   implicit none
   integer :: n,lwork,info
   complex(8) :: mat(n,n),evs(n),ma(n,n)
   complex(8) :: vl(n,n),work(3*n),vr(n,n)
   real(8) :: rwork(2*n)
   ma=mat
   lwork=3*n
   call zgeev('Na','Na',n,ma,n,evs,vl,n,vr,n,work,3*n,rwork,info)
!~ 	write(*,*) info
end subroutine

program main
use exact_diag_simple
implicit none
call random_seed()
   allocate(vij(Nv))
   call random_number(vij)
   vij=8.0*(vij-0.5)
   
   call exact_4sites!(vij)
end program
