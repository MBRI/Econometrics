function st=mbbq(x,freq,ndfy,ndfm,ndly,ndlm,lg)
% x data
% lg is data in logarithem

%frequency
%freq = 1;       % 1 for quarterly, 2 for monthly %

%
% ndfy = 1948;             % First Year of data set %
% ndfm = 1;                % First quarter/month of data set %
% ndly = 2004;             % Last Year of data set %
% ndlm = 2;                % Last quarter/month of data set %
%***********************************************************************************************/
%*                                                                                             */
%*  Modified BBQ program       MATLAB version                                                  */
%*                                                                                             */
%*  Computes turning points and imposes restrictions in one step                               */
%*                                                                                             */
%*    NB: rule is peak greater than turnphase on either side,                                  */
%*             trough less than turnphase on either side.           (can be modified)          */
%*    Restriction: Phase quarters/months, Cycle  quarters/months Alternating Troughs and Peaks,*/
%*                  if two peaks(troughs) in a row chooses highest(lowest), Trough must be     */
%*                  lower then preceeding peak, and No truning point within phase length       */
%*                  of end points                                                              */
%*                                                                                             */
%*                                                                                              */
%*                                                                                             */
%*    Date: 15th November 2005                                                                  */
%*     Author: James Engel - Code modified from Adrian Pagan and Don Harding's BBQ code        */
%***********************************************************************************************/
clc;

clock1 = clock;

%if simulating will want to set no reps. If analysing data nrep=1

nrep=1;% set to one if analysing real data  - line 112 to enter data/model %
complete = 1;  % switch: 1-only use complete cycles,0-use incomplete cycles (excess still on complete cycle) %

%cycle characteristics%...

turnphase = 2;
phase = 2;          % censoring rules %
cycle=5;
thresh=10.4;         % bypasses phase and cycle restriction if peak to trough is > than thresh %




nsfy=ndfy;nsfm=ndfm;nsly=ndly;nslm=ndlm;


if freq==1;
    nd= 4*(ndly-1-ndfy) + (5-ndfm) + (ndlm);     % Number of data points %
elseif freq==2;
    nd= 12*(ndly-1-ndfy) + (13-ndfm) + (ndlm);
end;


nn=nd;


notentp=0;
durs=zeros(nd,1);bcp5=zeros(nd,1);bct5=zeros(nd,1);st=zeros(nd,1);
pds=zeros(1,1);tds=zeros(1,1);pdsa=zeros(1,1);tdsa=zeros(1,1);bvec=zeros(nd,1);
pdmmat=zeros(nrep+1,1);tdmmat=zeros(nrep+1,1);
pdmmata=zeros(nrep,1);tdmmata=zeros(nrep,1);


%set icensor%...
%icensor=0 if no censoring=1 if is%...

icensor=1;
iseed=45891546;
iseed=7654;


ibomb=1;
%ibomb is set at replication before bombed out%...




pdcv=0;tdcv=0;pacv=0;tacv=0;pdm=0;pdma=0;tdm=0;tdma=0;pdcm=0;tdcm=0;
pdem=0;tdem=0;tdema=0;pdema=0;epcv=0;etcv=0;


nbt=0;
nbp=0;
iter=1;




disp('no iters')
disp(nrep)

while iter<=nrep;
    
    
    
    nd=nn;
    
    
    %data - has to be in LN form unless log command is not commented out
    
    %x=load('ozr.txt');
    %x=load('ozgdp.txt');
    %x=load('unem1.txt');
    %x=load('usrfull.txt');
    if lg==0
        x=log(x);
    end
    %simulated data - place model m file in working directory
    %cosim and usar are in zip files usar simulates from ar process
    
    %x = cosim(nd);
    %[x,iseed] = usar(iseed,nd);
    
    
    y=zeros(nd,1);c=zeros(nd,1);in=zeros(nd,1);
    nbp=0;nbt=0;
    
    [bcp5,bct5,nbp,nbt]=rawall(x(1:nd),turnphase,nd,phase,cycle,thresh);   % calculates turning points with restrictions %
    
    
    if nbp+nbt<=2;
        notentp=notentp+1;
    else;
        
        ntr=nbt;npk=nbp;
        
        nr=[nbt;nbp];
        nv=max(nr);
        
        pdc=zeros(nv,1);tdc=zeros(nv,1);
        pda=zeros(nv,1);tda=zeros(nv,1);td=zeros(nv,1);pd=zeros(nv,1);
        pdc=zeros(nv,1);tdc=zeros(nv,1);pde=zeros(nv,1);tde=zeros(nv,1);
        tdea=zeros(nv,1);pdea=zeros(nv,1);
        % calculate peak to trough durations & amps
        % p st&s for peaks,t for troughs,p gives cntractions,t expansions
        % code is pd,td is durations; pda,tda is amps; pdc,tdc is cum move pde, tde exces
        % excess is measured differently to avoid case that amps is close to zero in part cycle so that denom can become neg.so use cum movements as denom
        % now vs triangle in early paper
        % calculate peak to trough durations & amps
        % p st&s for peaks,t for troughs,p gives cntractions,t expansions
        % code is pd,td is durations; pda,tda is amps; pdc,tdc is cum move pde,tde exces
        
        if bcp5(1,1) < bct5(1,1);       % Peaks are first
            
            nr=[nbt;nbp];
            r=nbt;
            pd=bct5(1:r,1)-bcp5(1:r,1);
            pda=x(bct5(1:r,1))-x(bcp5(1:r,1));
            
            k=1;
            while k<=r;
                pdc(k)=sumc(x(bcp5(k,1):bct5(k,1),1)-x(bcp5(k,1),1));
                k=k+1;
            end;
        else;                      % troughs are First
            
            r=nbt-1;
            pd=bct5(2:r+1,1)-bcp5(1:r,1);
            pda=x(bct5(2:r+1,1))-x(bcp5(1:r,1));
            
            k=1;
            while k<=r;
                pdc(k)=sumc(x(bcp5(k,1):bct5(k+1,1),1)-x(bcp5(k,1),1));
                k=k+1;
            end;
            
            r1=r;
            
        end;
        
        % calculate trough to peak durations & amplitudes
        
        if bct5(1,1) < bcp5(1,1);        %  Troughs are first
            r=nbp;
            td=bcp5(1:r,1)-bct5(1:r,1);
            tda=x(bcp5(1:r,1))-x(bct5(1:r,1));
            k=1;
            while k<=r;
                tdc(k)=sumc(x(bct5(k,1):bcp5(k,1),1)-x(bct5(k,1),1));
                
                k=k+1;
            end;
            
        else;                      % peaks are First
            r=nbp-1;
            td=bcp5(2:r+1,1)-bct5(1:r,1);
            tda=x(bcp5(2:r+1,1))-x(bct5(1:r,1));
            
            
            k=1;
            while k<=r;
                
                tdc(k)=sumc(x(bct5(k,1):bcp5(k+1,1),1)-x(bct5(k,1),1));
                
                k=k+1;
            end;
            
        end;
        pdc=pdc(1:rows(pd));
        tdc=tdc(1:rows(td));
        
        % compute excesses
        %excess is percentage of triangle area
        za=(pd.*pda)/2;
        pde=100*(pdc-za-.5*pda)./za;
        
        pdea=100*(pdc-((pd.*pda)/2))./za;
        
        za=(td.*tda)/2;
        tde=100*(tdc-za-.5*tda)./za;
        
        tdea=100*(tdc-((td.*tda)/2))./za;
        
        %***********************************************************************************************/
        
        
        
        if complete == 0;    % switch: 1-only use complete cycles,0-use incomplete cycles (excess still on complete cycle)%
            
            bct5u=bct5;
            bcp5u=bcp5;
            
            if bcp5(1,1) < bct5(1,1);     % modifies code to include incomplete cycles %
                
                bct5u  = [1;bct5];
                
            elseif 	bct5(1,1) < bcp5(1,1);
                
                bcp5u = [1;bcp5];
                
            end;
            
            
            nbtu = rows(bct5u);
            nbpu = rows(bcp5u);
            
            if bcp5u(nbpu,1) < bct5u(nbtu,1);    % modifies code to include incomplete cycles %
                
                bcp5u  = [bcp5u;nd];
                
            elseif 	bct5u(nbt,1) < bcp5u(nbp,1);
                
                bct5u = [bct5u;nd];
                
            end;
            
            
            
            
            ntr=rows(bct5u);npk=rows(bcp5u);
            nr=[ntr;npk];
            nv=max(nr);
            
            
            pdc=zeros(nv,1);tdc=zeros(nv,1);
            pda=zeros(nv,1);tda=zeros(nv,1);td=zeros(nv,1);pd=zeros(nv,1);
            pdc=zeros(nv,1);tdc=zeros(nv,1);
            
            
            if bcp5u(1,1) < bct5u(1,1);       % Peaks are first
                
                nr=[ntr;npk];
                r=ntr;
                pd=bct5u(1:r,1)-bcp5u(1:r,1);
                pda=x(bct5u(1:r,1))-x(bcp5u(1:r,1));
                k=1;
                while k<=r;
                    pdc(k)=sumc(x(bcp5u(k,1):bct5u(k,1),1)-x(bcp5u(k,1),1));
                    k=k+1;
                end;
            else;                      % troughs are First
                r=ntr-1;
                pd=bct5u(2:r+1,1)-bcp5u(1:r,1);
                pda=x(bct5u(2:r+1,1))-x(bcp5u(1:r,1));
                k=1;
                while k<=r;
                    pdc(k)=sumc(x(bcp5u(k,1):bct5u(k+1,1),1)-x(bcp5u(k,1),1));
                    k=k+1;
                end;
                
                r1=r;
                
            end;
            
            % calculate trough to peak durations & amplitudes
            
            if bct5u(1,1) < bcp5u(1,1);        %  Troughs are first
                
                r=npk;
                td=bcp5u(1:r,1)-bct5u(1:r,1);
                tda=x(bcp5u(1:r,1))-x(bct5u(1:r,1));
                k=1;
                while k<=r;
                    tdc(k)=sumc(x(bct5u(k,1):bcp5u(k,1),1)-x(bct5u(k,1),1));
                    
                    k=k+1;
                end;
                
                
            else;                      % peaks are First
                r=npk-1;
                td=bcp5u(2:r+1,1)-bct5u(1:r,1);
                tda=x(bcp5u(2:r+1,1))-x(bct5u(1:r,1));
                
                
                k=1;
                while k<=r;
                    
                    tdc(k)=sumc(x(bct5u(k,1):bcp5u(k+1,1),1)-x(bct5u(k,1),1));
                    
                    k=k+1;
                end;
                
            end;
            pdc=pdc(1:rows(pd));
            tdc=tdc(1:rows(td));
            
        end;
        
        
        
        %***********************************************************************************************/
        
        
        % cumulate.when nrep=1 then gives raw data,otherwsise sums over monte carlo amounts
        
        pdcv=pdcv+stdc(pd)/meanc(pd);       % compute cv's
        tdcv=tdcv+stdc(td)/meanc(td);
        pacv=pacv+stdc(pda)/meanc(pda);
        tacv=tacv+stdc(tda)/meanc(tda);
        epcv = epcv + stdc(pde)/meanc(pde);  %  cv of xss
        etcv=etcv+stdc(tde)/meanc(tde);
        
        pdm=pdm+meanc(pd);          % durations
        pdma=pdma+meanc(pda);
        
        
        tdm=tdm+meanc(td);          % durations
        tdma=tdma+meanc(tda);
        
        
        pdcm=pdcm+meanc(pdc);       % cumulative
        tdcm=tdcm+meanc(tdc);
        
        % pdem=pdem+meanc(pde);
        tdem=tdem+meanc(tde);       % xss
        tdema=tdema+meanc(tdea);
        pdem=pdem+meanc(pde);
        pdema=pdema+meanc(pdea);
        
        
    end;
    
    if nrep==1.&(nbp+nbt>2);
        
        %nbp;nbt;%  % number of peaks and number of troughs %
        
        
        disp('peaks at')
        disp(bcp5(1:nbp))
        
        disp('troughs at')
        disp(bct5(1:nbt))
        
        
    end;
    
    
    pdmmat(iter+1)=pdm-pdmmat(iter);
    tdmmat(iter+1)=tdm-tdmmat(iter);
    iter=iter+1;
    
    
end;






if nrep==1.&(nbp+nbt>2);
    %remove below if want states ...
    st=zeros(nd,1);
    
    
    [st]=states(bcp5,bct5,nbp,nbt,nd);
    
    %determine obs in which states have been completed%...
    na=min([bct5(1);bcp5(1)])';
    nb=max([bct5(nbt);bcp5(nbp)])';
    nb-na+1;
    z=[x(na:nb) st(na:nb)];
    z1= [x st];
    
    disp('dated series')
    disp(z1)
    
end;


nrep1=nrep-notentp;

disp('statistics on average cycle')
disp('contractions/expansions')
disp('   ')

disp('durations')
disp([pdm/nrep1 tdm/nrep1])

disp('amplitudes')
disp([pdma/nrep1 tdma/nrep1])

disp('cumulative movements')
disp([pdcm/nrep1 tdcm/nrep1])

disp('excess movements percent of triangle area')
disp([pdem/nrep1 tdem/nrep1])

disp('cv of dur')
disp([pdcv/nrep1 tdcv/nrep1])

disp('cv of amps')
disp([pacv/nrep1 tacv/nrep1])

disp('cv of xss')
disp([epcv/nrep1 etcv/nrep1])

disp('no of its skipped since no peaks+troughs<=2')
disp(notentp)
pa=sortrows(pdmmat,1);
pb=sortrows(tdmmat,1);

%disp(pb)

clock2=clock;

timerun = clock2-clock1;
disp('running time')
disp(etime(clock2,clock1))