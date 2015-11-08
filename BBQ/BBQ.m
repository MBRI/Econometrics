function Dt=BBQ()%Dt,Equator)
% Put the Data in Data.xls
% the First Column must be data in eviews format and entitled "Date"
% the data must be sorted according to Date
lg=1; % if the data come in logarithmic scale

Dt=dataset('xlsfile', 'Data.xlsx');
addpath([cd '/bbq_000']);


%% find date
QQ1=char(Dt.Date(1));
QQ2=char(Dt.Date(end));

if strcmp(QQ1(5),'Q')
    freq=1;
else
    freq=2;
end

Fy=cell2num(cellstr(QQ1(1:4)));
FD=cell2num(cellstr(QQ1(6:end)));
Ly=cell2num(cellstr(QQ2(1:4)));
LD=cell2num(cellstr(QQ2(6:end)));

%%
for j=2:size(Dt,2)
X=double(Dt(:,j));
st=mbbq(X,freq,Fy,FD,Ly,LD,lg); % I changed the mfile from the Original Source

Dt.(['Booms_' Dt.Properties.VarNames{j}])=st;

%% Ploting
ss=st-lagmatrix(st,1);
Pk=find(ss==1);
tr=find(ss==-1);
a=size(Pk,1)-size(tr,1);
if a==1
    tr(end+1)=size(X,1)+1;
elseif a==-1
    Pk(end+1)=2;
    Pk=sort(Pk);
end
z=[Pk-1,tr-1];
mX=min(X);
MX=max(X);
mm=repmat([mX, MX , MX,mX],size(z,1),1);
zz=[z(:,1),z(:,1),z(:,2),z(:,2)];
grbkgrnd = [.8 .8 .8];
%
figure;
hold on
h = patch(zz.',mm.',grbkgrnd);
set(h,'linestyle','none')
plot(X,'black');
% %xlabel(Dt.date);
axis([0 length(X) mX MX])

hold off
end
export(Dt,'xlsfile','dd')
end
function [outputmat]=cell2num(inputcell)
% this function is modified from the first usage
% Function Force to convert all numeric cell array to a double precision array
% ********************************************
% Usage: outputmatrix=cell2num(inputcellarray)
% ********************************************
% Output matrix will have the same dimensions as the input cell array
% Non-numeric cell contest will become 0 outputs in outputmat
% This function only works for 1-2 dimensional cell arrays

if ~iscell(inputcell), error('Input cell array is not.'); end

outputmat=zeros(size(inputcell));

for c=1:size(inputcell,2)
    for r=1:size(inputcell,1)
        %if isnumeric(inputcell{r,c})
        try
            outputmat(r,c)=str2num(inputcell{r,c});
        catch  %#ok<CTCH>
            a=inputcell{r,c};% get the unconvertable string
            if isa(a,'char')
                if isempty(a)
                    c1=0;
                else
                    for i=1:length(a)
                        b(i)=unicode2native(a(i));%#ok<*AGROW> %convert to ascii code
                    end
                    b(b<48)=48; %ignore the non numberic letter
                    b(b>57)=48; %ignore the non numberic letter
                    b=b-48; % set to nubmer
                    c1=double(0);
                    for i=1:length(b)
                        c1=c1*10+double(b(i));% sum to one nubmer
                    end
                end
            else
                c1=a;
            end
            %   disp(['Data in row ' num2str(r) ': col ' num2str(c) ' removed by ' num2str(c1)]);
            outputmat(r,c)=c1;
        end
    end
end

end
% function CP=UT(X,Equator)
%
% CP=[]; % critical point
% X(X<Equator)=0;
% while(max(X)>0)
%     A=max(X);
%     [A, B]=find(X==A);
%     CP=[CP; A];
%     for l=1:length(A)
%         i=A(l);
%         if i<=size(X,1)
%             while(X(i)>0 )%|| i-A(l)<=MinLeng )
%                 X(i)=0;
%                 i=i+1;
%                 if i>size(X,1)
%                     break;
%                 end
%             end
%         end
%
%         i=A(l)-1;
%         if i>=1
%             while(X(i)>0)%|| A(l)-i<=MinLeng )
%                 X(i)=0;
%
%                 i=i-1;
%                 if i<1
%                     break;
%                 end
%             end
%         end
%     end
%
% end
% CP=sort(CP);
% end
%
% function Z=Sorter(Peaks,Troughs)
% MinLeng=2;
% Troughs(Troughs<Peaks(1))=[];
% Peaks(Peaks>Troughs(end))=[];
% A= length(Peaks)-length(Troughs);
%
% if  length(Peaks)~=length(Troughs)
%
%     error('Problem A~=0');
% end
%
% Z=Troughs-Peaks;
% if sum(Z<0)>0
%     warning('Problem in Peaks order');
% end
% if sum(Z<MinLeng)>0
%     warning('Problem in Cycle Length');
%     %Troughs(Z<3)=[];
%     %Peaks(Z<3)=[];
% end
% Z=[Peaks,Troughs];
%
% end
%
%
% function [Y,Dat]=modifyData(Dt, DeS, Lg, Hp)
%
% addpath('E:/MINE/P11_Monetary Mechanism/Matlab/SVECX2/IRIS_Tbx');
% irisstartup
% %  Date
% Gd=1; Dd=2;
% %irisstartup
% Mn=min(double(Dt(:,Gd)));
% Mx=max(double(Dt(:,Gd)));
% q=qq(fix(Mn),(Mn-fix(Mn))*4+1):qq(fix(Mx),(Mx-fix(Mx))*4+1);
% q=q';
% %---------------------------Identify NaN-----------------------------------
% %[~,N]=size(data);
% Dt=tseries(q,double(Dt(:,Dd)));
% Dat=Dt;
% if DeS
% Dt=x12(Dt);
% end
%
% if Lg
% Dt=log(Dt);
% end
% %----------------
% %-------------------
% if Hp
% [~,Dt,~,~]=hpf(Dt,inf,'Lambda',677);
% end
%
% for k=1:1
% [Dt,~,~,~]=hpf(Dt,inf,'Lambda',1);
% end
% %Y=double(Dt(:,2))-double(Dt_C(:,2));
% % the
% plot(Dt)
% Y=Dt.data;
% irisstartup -shutup
%
% %Y=Dt(:,2);
% end
%
