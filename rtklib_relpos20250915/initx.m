function [rtkX,rtkP]=initx(rtknx,rtkX_,rtkP_,xi,var,i)
rtkX=rtkX_;
rtkP=rtkP_;
rtkX(i)=xi;
for j=1:rtknx
    if i==j
        rtkP(i,j)=var;
    else
        rtkP(i,j)=0;
        rtkP(j,i)=0;
    end
end
