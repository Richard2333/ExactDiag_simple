program exact_4sites
   implicit none
   
   integer :: elecs(256,4,2),mm(8),btom1(0:5,2),btom2(0:5,2)
   integer :: atom1(4,2),atom2(4,2)
   integer :: i,j,k,z, m,n,a1,a2
   real :: t,u0,u1
   complex(8)  :: hrs(-2:2,-2:2,-2:2)
   real(8) :: ham_hub_2(256,256)
   
   elecs=0
   k=0
   do k=1,256
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
      elecs(k,:,1)=mm(1:4)
      elecs(k,:,2)=mm(5:8)
   enddo
   ham_hub_2(i,j)=0
   do i=1,256
      do j=1,256
         btom1(1:4,:)=elecs(i,1:4,:)
         btom1(0,:)=elecs(i,4,:)
         btom1(5,:)=elecs(i,1,:)
         btom2(1:4,:)=elecs(j,1:4,:)
         btom2(0,:)=elecs(j,4,:)
         btom2(5,:)=elecs(j,1,:)
         do m=1,4
            do n=1,2
               if(btom1(m,n)*btom2(m-1,n)==1) ham_hub_2(i,j)=ham_hub_2(i,j)+t
               if(btom1(m,n)*btom2(m+1,n)==1) ham_hub_2(i,j)=ham_hub_2(i,j)+t
            enddo
         enddo
      enddo
   enddo
         
      
end program
