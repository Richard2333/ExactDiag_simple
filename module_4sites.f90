
module exact_diag_simple

   real(8),allocatable :: vij(:)
   integer,parameter :: Nv=4
   integer :: indom,oudom

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
   subroutine exact_4sites(indx_s)
!~    use wmpi
!~    use para
!~    use
      implicit none

      integer :: indx_s
!~    integer :: atom1(Na,2),atom2(Na,2)
      integer :: i,j,k,z, m,n,Na
      real(8) :: t,u2,u1
      real(8) :: correlator(Nv,Nv),rho(Nv),psi2,filling

      real(8),allocatable :: evs(:)
      complex(8),allocatable :: ham_hub_2(:,:),gspsi(:)
      integer,allocatable :: elecs(:,:),mm(:),btom1(:),btom2(:),as(:)

      integer :: indom,oudom
!~    integer :: ierr,mpii,nsets


!~    read(*,*) Na
      Na=2
      Na=Nv*Na


      allocate(evs(2**Na))
      allocate(ham_hub_2(2**Na,2**Na),gspsi(2**Na))
      allocate(elecs(2**Na,Na))
      allocate(mm(Na),as(-3:Na+3))
      allocate(btom1(0:Na+1),btom2(0:Na+1))

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
!~       write(*,*) k,elecs(k,:)
      enddo
      ham_hub_2=0.0d0
      t= -1.0
      u1=1.0
      u2=0.5
      do n=-3,Na+3
         as(n)=mod(n+Na*5-1,Na)+1
!~       write(*,*) n,as(n)
      enddo
      do i=1,2**Na
         do j=1,2**Na
            ham_hub_2(i,j)=(0.0d0,0.0d0)
            btom1=[elecs(i,Na),elecs(i,1:Na),elecs(i,1)]
            btom2=[elecs(j,Na),elecs(j,1:Na),elecs(j,1)]

            do m=1,Na
               if(btom1(m)==1 .and. i==j) ham_hub_2(i,j)=ham_hub_2(i,j)+vij(mod(m-1,Nv)+1)

               if(btom1(m)*btom2(as(m-1))==1) then
!~                ham_hub_2(i,j)=ham_hub_2(i,j)+u1
                  if(btom1(as(m-1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
               endif
               if(btom1(m)*btom2(as(m+1))==1) then
                  ham_hub_2(i,j)=ham_hub_2(i,j)+u1
                  if(btom1(as(m+1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
               endif
!~             if(btom1(m)*btom2(as(-2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2
               if(btom1(m)*btom2(as(m+2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2

!~             enddo
            enddo
         enddo
      enddo
      do i=1,2**Na
         do j=1,2**Na
            if(abs( ham_hub_2(i,j))>0.1) k=k+1
         enddo
      enddo
!~       write(*,*) k,1.0*k/4**Na
!~       write(*,*) sum(ham_hub_2-conjg(transpose(ham_hub_2)))
!~    stop
      i=2**Na
      call eigenvalh(i,ham_hub_2,evs)
      gspsi=ham_hub_2(:,1)
!~    write(*,*) 'EGS=',evs(1),evs(2)
!~    do i=1,2**Na
!~       write(*,*) ham_hub_2(i,1),ham_hub_2(i,2)
!~    enddo
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
         enddo
      enddo
      write(oudom+indx_s,*) 'EGS',evs(1)
      write(oudom+indx_s,*) 'correlator'
      write(indom+indx_s,*) 'ne',Nv*filling/Na
      write(indom+indx_s,*) 'rho           ','vij'
      do i=1,Nv
         write(indom+indx_s,*) rho(i),vij(i)
         write(oudom+indx_s,*) (correlator(i,j),j=1,Nv)
      enddo
      write(*,*) filling/sum(abs(gspsi)**2),Nv*filling/Na
!~    enddo

   end subroutine
end module
