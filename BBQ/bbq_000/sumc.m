function m=sumc(x);
% get sum of each column
m=sum(x);
if size(m,2)>1;
   m=m';
end;
return;
