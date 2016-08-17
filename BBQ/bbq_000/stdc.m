function m=stdc(x);
% get standard deviation of each column
m=std(x);
if size(m,2)>1;
   m=m';
end;
return;
