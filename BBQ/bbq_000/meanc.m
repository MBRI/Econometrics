function m=meanc(x);
% get mean of each column
m=mean(x);
if size(m,2)>1;
   m=m';
end;
return;
