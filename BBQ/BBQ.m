function Dt=BBQ()%Dt,Equator)
% Put the Data in Data.xls
% the First Column must be data in eviews format and entitled "Date"
% the data must be sorted according to Date
lg=1; % if the data come in logarithmic scale

Dt=dataset('xlsfile', 'Data.xlsx');
addpath([cd '/bbq_000']);


%% find date
QQ1=Dt.Date{1};
QQ2=Dt.Date{end};

if strcmp(QQ1(5),'Q')
    freq=1;
else
    freq=2;
end

Fy=str2num(QQ1(1:4));
FD=str2num(QQ1(6:end));
Ly=str2num(QQ2(1:4));
LD=str2num(QQ2(6:end));

%% init figure
%parameters for figure and panel size
NumberofPlot=size(Dt,2)-1;
plotheight=20;
plotwidth=16;
subplotsx=floor(NumberofPlot^0.5);
subplotsy=ceil(NumberofPlot/subplotsx);
leftedge=1.2;
rightedge=0.4;
topedge=1;
bottomedge=1.5;
spacex=1;
spacey=1;
fontsize=10;
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

%setting the Matlab figure
f=figure('visible','on');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);



%%

for i=1:subplotsx
    for ii=1:subplotsy
        j=(i-1)*subplotsx+ii+1;% the one is for Date Column in the begining of Data
        if j>NumberofPlot+1% the one is for Date Column in the begining of Data
            break;
        end
        %for j=2:size(Dt,2)
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
        %figure;
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top');
        hold on
        h = patch(zz.',mm.',grbkgrnd);
        set(h,'linestyle','none')
        plot(X,'black');
        title(Dt.Properties.VarNames{j});
        % %xlabel(Dt.date);
        axis([0 length(X) mX MX])
        if ii==subplotsy || j==NumberofPlot+1
           set(gca,'XTicklabel',Dt.Date.','fontsize',8,'XTickLabelRotation',90);  
        else
           set(gca,'XTick',[]); 
        end
        hold off
    end
end
Printer('Res')
export(Dt,'xlsfile','out\Res')
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

function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,nbx,nby,spacex,spacey)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subxsize=(plotwidth-leftmargin-rightmargin-spacex*(nbx-1.0))/nbx;
subysize=(plotheight-topmargin-bottommargin-spacey*(nby-1.0))/nby;

for i=1:nbx
    for j=1:nby
        
        xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
        yfirst=bottommargin+(j-1.0)*(subysize+spacey);
        
        positions{i,1+nby-j}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
        
    end
end


%setting the Matlab figure
f=figure('visible','off');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);
end

function Printer(filename)
Adr='OUT\';
%Adr='E:\MINE\PrjData_94_1\Code_matlab\OUT\img\';
if ~exist(Adr,'dir')
    mkdir(Adr);
end
filename=[Adr filename];
print(gcf, '-depsc2','-loose',[filename,'.eps']);
system(['epstopdf ',filename,'.eps']);
%system(['convert -density 300 ',filename,'.eps ',filename,'.png'])
close gcf;
end