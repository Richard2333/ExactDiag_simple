
module exact_diag_simple

   real(8),allocatable :: vij(:)
   integer,parameter :: Nv=4
   integer :: indom,oudom,Na

contains

   subroutine eigenvalh(n,mat,evs)
      implicit none
      integer :: n,lwork,info
      complex(8) :: mat(n,n)!,ma(n,n)
      complex(8) :: work(3*n)
      real(8) :: rwork(3*n-2),evs(n)
!~    ma=mat
      lwork=3*n
      call zheev('V','L',n,mat,n,evs,work,3*n,rwork,info)
!~ 	write(*,*) info
   end subroutine
   subroutine exact_4sites!(indx_s)
!~    use wmpi
!~    use para
!~    use
      implicit none

      integer :: indx_s
!~    integer :: atom1(Na,2),atom2(Na,2)
      integer :: i,j,k,z, m,n
      real(8) :: t,u2,u1
      real(8) :: correlator(Nv,Nv),rho(Nv),psi2,filling

      real(8),allocatable :: evs(:)
      complex(8),allocatable :: ham_hub_2(:,:),gspsi(:)
      integer,allocatable :: elecs(:,:),mm(:),btom1(:),btom2(:),as(:)

      integer :: indom,oudom
!~    integer :: ierr,mpii,nsets


!~    read(*,*) Na
!~       Na=2
!~       Na=Nv*Na


      allocate(evs(2**Na))
      allocate(ham_hub_2(2**Na,2**Na),gspsi(2**Na))
      allocate(elecs(2**Na,Na))
      allocate(mm(Na),as(-3:Na+3))
      allocate(btom1(1:Na),btom2(1:Na))

      indom=10000;oudom=20000

!~    call random_seed()
!~    call mpi_init(ierr)
!~    call mpi_comm_rank(mpi_cmw,cpuid,ierr)
!~    call mpi_comm_size(mpi_cmw,num_cpu,ierr)
!~    do mpii= 1+cpuid, 8, num_cpu

      elecs=0
      k=0
      do k=1,2**Na
         z=k-1
         mm=0
         do i=Na-1,0,-1
            if(z>=2**(i)) then
               mm(i+1)=1
               z=z-2**(i)
            else
               cycle
            endif
            if(z<=0) cycle
         enddo
         elecs(k,:)= mm
!~          write(933,*) k,elecs(k,:)
      enddo
      ham_hub_2=0.0d0
      t=-1.0
      u1=1.0
      u2=0.5

      do i=1,2**Na
         do j=1,2**Na
            ham_hub_2(i,j)=(0.0d0,0.0d0)
            btom1=elecs(i,:)
            btom2=elecs(j,:)
            ham_hub_2(i,j)=matrix_elem_hubbard(btom1,btom2,Na,vij,t,u1,u2)
            
         enddo
      enddo


      i=2**Na
!~       do j=1,i
!~          write(299,*) (real(ham_hub_2(j,m)),m=1,2**Na)
!~       enddo
      call eigenvalh(i,ham_hub_2,evs)
!~       call eigensystem_c('V','L',i,ham_hub_2,evs)
!~       gspsi=ham_hub_2(1,:)
      gspsi=ham_hub_2(:,1)
!~       write(*,*) 'diff=',sum(abs(ham_hub_2(:,1)-conjg(ham_hub_2(1,:))))
!~    stop
      correlator=0
      rho=0
      filling=0
      do i=1,2**Na
         psi2=abs(gspsi(i))**2
         do j=1,Na
            filling=filling+psi2*elecs(i,j)
            rho(mod(j-1,Nv)+1)=rho(mod(j-1,Nv)+1)+psi2*elecs(i,j)
            do k=1,Na
               correlator(mod(j-1,Nv)+1,mod(k-1,Nv)+1)=correlator(mod(j-1,Nv)+1,mod(k-1,Nv)+1)+psi2*elecs(i,j)*elecs(i,k)
            enddo
!~             filling=filling+psi2*elecs(i,j)
         enddo
      enddo
      write(*,*) evs(1),filling/sum(abs(gspsi)**2),Nv*filling/Na
      
      deallocate(evs)
      deallocate(ham_hub_2,gspsi)
      deallocate(elecs)
      deallocate(mm,as)
      deallocate(btom1,btom2)

!~    enddo

   end subroutine

   function matrix_elem_hubbard(btom1,btom2,Na,vij,t,u1,u2) result(hij)
      implicit none
      complex(8) :: hij
      integer :: btom1(Na),btom2(Na),Na,as(-3:Na+3),n
      real(8),intent(in) :: vij(Nv),t,u1,u2
      integer :: i,j,k,bdiff(Na),bdif1,bdif2

      do n=-3,Na+3
         as(n)=mod(n+Na*5-1,Na)+1
!~       write(*,*) n,as(n)
      enddo
      hij=0.00d0
      bdiff=btom1-btom2
      if (sum(abs(bdiff))==0) then
         do i=1,Na
            if(btom1(i)==1) then
               hij=hij+vij(mod(i-1,Nv)+1)
               if(btom1(as(i+1))==1) hij=hij+u1
               if(btom1(as(i+2))==1) hij=hij+u2
            endif
         enddo
      else
         do i=1,Na
            bdif1=sum(abs(bdiff))-abs(bdiff(i))-abs(bdiff(as(i+1)))
            bdif2=sum(abs(bdiff))-abs(bdiff(i))-abs(bdiff(as(i-1)))
            if(btom1(i)*btom2(as(i+1))==1 .and. btom2(i)*btom1(as(i+1))==0 .and. bdif1==0) hij=hij+t
            if(btom1(i)*btom2(as(i-1))==1 .and. btom2(i)*btom1(as(i-1))==0 .and. bdif2==0) hij=hij+t
         enddo
      endif

      return
   end function

end module
