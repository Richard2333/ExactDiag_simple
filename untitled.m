atoms=zeros(256,4,2);
for i=1:256
    dec2bin(i-1)
    A=bitget(uint8(i-1),1:8);
    atoms(i,:,1)=A(1:4);
    atoms(i,:,2)=A(5:8);
end

Ham_hub1=zeros(256,256);
u0=0.2;u1=0.4;t=1.0;
btom1=zeros(6,2)
btom2=zeros(6,2)

for i=1:256
    for j=1:256
        btom1(1,1:2)=atoms(i,4,:)
        btom1(2:5,:)=atoms(i,:,:)
        btom1(6,:)=atoms(i,1,:)
        btom2(1,1:2)=atoms(j,4,:)
        btom2(2:5,:)=atoms(j,:,:)
        btom2(6,:)=atoms(j,1,:)
        for a1=1:4
            for s1=1:2
                if(atoms(i,a1,s1)==1 && atoms(j,a1-1,s1)==0) 
                    Ham_hub1(i,j)=Ham_hub1(i,j)+t;
                end
                if(atoms(i,a1,s1)==1 && atoms(j,a1+1,s1)==0) 
                    Ham_hub1(i,j)=Ham_hub1(i,j)+t;
                end
                if(atoms(i,a1,s1)==0 && atoms(j,a1-1,s1)==1) 
                    Ham_hub1(i,j)=Ham_hub1(i,j)+t;
                end
                if(atoms(i,a1,s1)==0 && atoms(j,a1+1,s1)==1) 
                    Ham_hub1(i,j)=Ham_hub1(i,j)+t;
                end
            end
        end
    end
end

