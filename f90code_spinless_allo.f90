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
subroutine eigenvalh(n,mat,evs)
   implicit none
   integer :: n,lwork,info
   complex(8) :: mat(n,n)!,ma(n,n)
   complex(8) :: vl(n,n),work(3*n),vr(n,n)
   real(8) :: rwork(3*n-2),evs(n)
!~    ma=mat
   lwork=3*n
   call zheev('V','L',n,mat,n,evs,work,3*n,rwork,info)
!~ 	write(*,*) info
end subroutine
program main
   use wmpi
   use exact_diag_simple
!~    use para
!~    use
   implicit none
!~    integer,parameter :: Nv=4
!~    integer :: atom1(Na,2),atom2(Na,2)
   integer :: i,j,k,z, m,n,as(-3:3),Na
   real(8) :: t,u2,u1

   real(8),allocatable :: evs(:)
   complex(8),allocatable :: ham_hub_2(:,:),gspsi(:)
   integer,allocatable :: elecs(:,:),mm(:),btom1(:),btom2(:)

!~    integer :: indom,oudom
   integer :: ierr,mpii,nsets


   read(*,*) Na,nsets
   Na=Nv*Na

   allocate(vij(Nv))
   allocate(evs(2**Na))
   allocate(ham_hub_2(2**Na,2**Na),gspsi(2**Na))
   allocate(elecs(2**Na,Na))
   allocate(mm(Na))
   allocate(btom1(0:Na+1),btom2(0:Na+1))

   indom=10000;oudom=20000

   call random_seed()
   call mpi_init(ierr)
   call mpi_comm_rank(mpi_cmw,cpuid,ierr)
   call mpi_comm_size(mpi_cmw,num_cpu,ierr)
   do mpii= 1+cpuid, nsets, num_cpu
      call random_number(vij)
      vij=4.0*(vij-0.5)
      call exact_4sites(mpii)
   enddo

end program
