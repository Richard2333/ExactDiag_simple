xs=0;
ys=0;
xx=0;
yx=0;
for i=-40:40
    for j=-40:40
        r=sqrt(i*i+j*j+1);
        x=i/(r+1);
        y=j/(r+1);
        xs=[xs,x];
        ys=[ys,y];
%         q=sqrt(4*r*r)-1;
%         q=q/(2*r*r);
%         x=i*q;
%         y=j*q;
%         xx=[xx,x];
%         yx=[yx,y];
    end
end
plot(xs,ys,'.')
% plot(xs,ys,'.r',xx,yx,'.b')