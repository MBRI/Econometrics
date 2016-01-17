function [Y1,Y2,Y3,Y4]=remnan(Y1,Y2,Y3,Y4)
% This function would remove all rows of all matrices wich does not contian
% a finite number.
% Check the compatibility of inputs
[i1 ,~]=find(isfinite(Y1)==0);
[i2 ,~]=find(isfinite(Y2)==0);
[i3 ,~]=find(isfinite(Y3)==0);
[i4 ,~]=find(isfinite(Y4)==0);

ii=unique([i1;i2;i3;i4]);

Y1(ii,:)=[];
Y2(ii,:)=[];
Y3(ii,:)=[];
Y4(ii,:)=[];

end
