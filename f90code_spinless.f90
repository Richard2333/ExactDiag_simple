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
   complex(8) :: mat(n,n),ma(n,n)
   complex(8) :: vl(n,n),work(3*n),vr(n,n)
   real(8) :: rwork(3*n-2),evs(n)
   ma=mat
   lwork=3*n
   call zheev('V','L',n,ma,n,evs,work,3*n,rwork,info)
!~ 	write(*,*) info
end subroutine
program exact_4sites
!~    use
   implicit none
   integer,parameter :: Na=9
!~    integer :: atom1(Na,2),atom2(Na,2)
   integer :: i,j,k,z, m,n,as(-3:3)
   real(8) :: t,u2,u1
   
   real(8) :: vij(Na),evs(2**Na)
   complex(8) :: ham_hub_2(2**Na,2**Na)
   integer :: elecs(2**Na,Na),mm(Na),btom1(0:Na+1),btom2(0:Na+1)

   
   
   
   
   call random_seed()

   call random_number(vij)
   vij=vij-0.5
   elecs=0
   k=0
   do k=1,2**Na
      z=k-1
      mm=0
      do i=7,0,-1
         if(z>=2**(i)) then
            mm(i+1)=1
            z=z-2**(i)
         else
            cycle
         endif
         if(z<=0) cycle
      enddo
      elecs(k,:)= mm
      write(*,*) k,elecs(k,:)
   enddo
   ham_hub_2=0.0d0
   t=-1.0
   u1=0.1
   u2=0.05
   do i=1,2**Na
      do j=1,2**Na
         ham_hub_2(i,j)=(0.0d0,0.0d0)
         btom1=[elecs(i,Na),elecs(i,1:Na),elecs(i,1)]
         btom2=[elecs(j,Na),elecs(j,1:Na),elecs(j,1)]

         do m=1,Na
            if(btom1(m)==1 .and. i==j) ham_hub_2(i,j)=ham_hub_2(i,j)+vij(m)
            do n=-3,3
               as(n)=mod(m+n+Na*5-1,Na)+1
            enddo
            if(btom1(m)*btom2(as(-1))==1) then
               ham_hub_2(i,j)=ham_hub_2(i,j)+u1
               if(btom1(as(-1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
            endif
            if(btom1(m)*btom2(as(1))==1) then
               ham_hub_2(i,j)=ham_hub_2(i,j)+u1
               if(btom1(as(1))+btom2(m)==0) ham_hub_2(i,j)=ham_hub_2(i,j)+t
            endif
            if(btom1(m)*btom2(as(-2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2
            if(btom1(m)*btom2(as(2))==1) ham_hub_2(i,j)=ham_hub_2(i,j)+u2

!~             enddo
         enddo
      enddo
   enddo
!~    do i=1,16
!~    do j=1,16
!~       write(*,*) i,j, ham_hub_2(i,j)
!~    enddo
!~    enddo
   write(*,*) sum(ham_hub_2-conjg(transpose(ham_hub_2)))
   i=2**Na
   call eigenvalh(i,ham_hub_2,evs)
   do i=1,2**Na
      write(*,*) evs(i)
      write(233,*) i,evs(i)
   enddo


end program
