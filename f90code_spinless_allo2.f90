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
program exact_4sites
   use wmpi
!~    use para
!~    use
   implicit none
   integer,parameter :: Nv=4
!~    integer :: atom1(Na,2),atom2(Na,2)
   integer :: i,j,k,z, m,n,as(-3:3),Na
   real(8) :: t,u2,u1

   real(8),allocatable :: vij(:),evs(:)
   complex(8),allocatable :: ham_hub_2(:,:),gspsi(:)
   integer,allocatable :: elecs(:,:),mm(:),btom1(:),btom2(:)

   integer :: indom,oudom
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
      vij=-4.0*(vij)
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
      t=-1.0
      u1=0.0
      u2=0.0
      do i=1,2**Na
         do j=1,2**Na
            ham_hub_2(i,j)=(0.0d0,0.0d0)
            btom1=[elecs(i,Na),elecs(i,1:Na),elecs(i,1)]
            btom2=[elecs(j,Na),elecs(j,1:Na),elecs(j,1)]

            do m=1,Na
               if(btom1(m)==1 .and. i==j) ham_hub_2(i,j)=ham_hub_2(i,j)+vij(mod(m,Nv)+1)
               do n=-3,3
                  as(n)=mod(m+n+Na*5-1,Na)+1
               enddo
               if(btom1(m)*btom2(as(-1))==1) then
!~                ham_hub_2(i,j)=ham_hub_2(i,j)+u1
                  if(btom1(as(-1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
               endif
               if(btom1(m)*btom2(as(1))==1) then
                  ham_hub_2(i,j)=ham_hub_2(i,j)+u1
                  if(btom1(as(1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
               endif
!~             if(btom1(m)*btom2(as(-2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2
               if(btom1(m)*btom2(as(2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2

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
      write(oudom+mpii,*) evs(1)
      do i=1,2**Na   
         write(233+mpii,*) gspsi(i)
      enddo
   enddo

end program
