function [x,iseed] = usar(iseed, nd);

nd=nd+500;

xx=zeros(nd,1);
%randn('state',iseed);
randn('state',sum(100*clock));
e=randn(nd,1);

u=e*0.0096;

i = 3;
xx(1) = 7.3;
xx(2) = 7.3;

while i <= nd;

xx(i) = 0.0055 + 1.3438*xx(i-1) - .3438*xx(i-2) +  u(i);

i = 1+i;
end;


x=xx(501:nd);
nd=nd-500;


% end of function