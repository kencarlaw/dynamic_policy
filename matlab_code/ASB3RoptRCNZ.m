% ASBcon2RoptRCNZ.m 
% generates data for ASBconV2Fig51Fig72.m specifically fig 7.2 over
% variations in N and z
% Aug. 31, 2022, K.I. Carlaw

clear

%parameters 

N=100;      %population of agents
NN=N+1;
%M=20;
MM=200;

Block=50000;

%T=MM*Block;     %Time/iteration index for the sim

F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
mu=0.75;     %mean value of gi, individual benefit from ASA
sig=0.28;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
SWrho=5*rho;
mm=0;
Z=2;        %z-history length
Disc=0.995;

bine=44;
Rgb=38;
Rbb=68;

dim=1;
Rc1=dim;
Rc2=dim;
bin=dim;
RS1=0;
RS2=0;
critconv1=0.01;

%BinE=zeros(bin,1);
%R2=zeros(Rc2,1);
%R1=zeros(Rc1,1);
cost=zeros(dim,1);
%cost2=zeros(Rc2, Rc1,1);
Ev=zeros(dim,1);
%sw=zeros(dim,1);
Eq=zeros(dim,1);
%diagq=zeros(NN,1);
Ea=zeros(dim,1);
%diaga=zeros(NN,1);
ER=zeros(dim,1);
Egbar=zeros(dim,1);
scount=zeros(dim,1);
%Eswi=zeros(Rc2,Rc1,1);
RC=zeros(dim,1);
edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end

%opts = spreadsheetImportOptions("NumVariables", 5);
%opts.Sheet = "Sheet1";
%opts.DataRange = "A2:E28";
%opts.VariableNames = ["BE", "R1", "R2", "mu1", "sig"];
%opts.VariableTypes = ["double", "double", "double", "double", "double"];
%tbl = readtable("C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\RCRAT2.xlsx", opts, "UseExcel", false);
%BE = tbl.BE;
%R1 = tbl.R1;
%R2 = tbl.R2;
%mu1 = tbl.mu1;
%sig = tbl.sig;
%clear opts tbl

BinE=55;
R1=50;
R2=58;

for k=1:dim
    %BinE(k)=k; %k-1;
    %for r2=1:dim %Rc2 %1:Rc2
        %R2(r2)=r2-1;
        %for r1=1:dim %1:r2
            %R1(r1)=r1-1;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
                vz=zeros(Block,1);
                az=zeros(Block,1);
                eta=zeros(Block,1);
                R=zeros(Block,1);
                gbar=zeros(Block,N);
                gbar3=zeros(Block,1);
                A=zeros(Block,1);
                NCw=zeros(Block,1);NCw2=zeros(Block,1);
                vw=zeros(Block,1);aw=zeros(Block,1);
                b=b+1;
                if b<2
                    in=Z+1;
                else
                    in=1;
                end
                for t=1:Block
                    if t>1
                        if vw(t-1)<BinE
                            R(t)=R1;
                        else
                            R(t)=R2;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(R2-R1); 
                            swic(t)=1;
                        end                    
                    end                
                    if (b<2)
                        if (t<Z+1)
                            for z=1:Z
                                vz(z)=Z*unifrnd(0,N);
                                az(z)=Z*unifrnd(0,vz(z));
                            end
                        else
                            vz(t)=sum(vw(t-Z:t-1));
                            az(t)=sum(aw(t-Z:t-1));
                        end
                    else
                        if t<1+Z
                            if t<2
                                vz(t)=sum(v(t+(b-1)*Block-Z:t+(b-1)*Block-1));
                                az(t)=sum(a(t+(b-1)*Block-Z:t+(b-1)*Block-1));
                            else
                                vz(t)=sum(v(t+(b-1)*Block-Z-1:(b-1)*Block))+sum(vw(t-(Z-(Z-(t-1))):t-1));
                                az(t)=sum(a(t+(b-1)*Block-Z-1:(b-1)*Block))+sum(aw(t-(Z-(Z-(t-1))):t-1));  
                            end
                        else
                            vz(t)=sum(vw(t-Z:t-1));
                            az(t)=sum(aw(t-Z:t-1));
                        end
                    end
                    eta(t)=(aa+az(t))/(aa+bb+vz(t));
                    for n=1:N
                        if eta(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                    swicc=swic;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    swihold=cat(2,swicc',swic');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
                    sw=swihold';
                end

                vconv=v(1:b*Block);
                vconvlag=v(1:(b-1)*Block);
                freqvconv=histcounts(vconv(:),edges)/((b)*Block);
                freqvconvlag=histcounts(vconvlag(:),edges)/((b-1)*Block);
                TESTCONV1=zeros(1,NN);
                if b>1
                    for j=1:NN
                        TESTCONV1(j)=abs(freqvconv(j)-freqvconvlag(j));
                    end
                    TESTCONV(b)=sum(TESTCONV1);                
                end
                if b > 4
                    if (TESTCONV(b)<=critconv1) && (TESTCONV(b-1)<=critconv1) ...
                        && (TESTCONV(b-2)<=critconv1) && (TESTCONV(b-3)<=critconv1)...
                        && (TESTCONV(b-4)<=critconv1)
                        c=b;
                        cc=1;
                    end
                end
            end
            cost(k)=mean(NC);
            %cost2(r2,r1)=mean(NC2);
            Ev(k)=mean(v);
            Eq(k)=mean(eta);
            Ea(k)=mean(a);
            ER(k)=mean(R);
            Egbar(k)=mean(gb);
            RC(k)=max(0,(ER(k)-Ev(k))/ER(k));
            scount(k)=sum(sw)/(b*Block);       
            %end   
        %R2(k)
    %end
    %BinE(k)
    k
end

%costh=cost;
%costh2=cost2;
%MH=max(max(cost));
%MH2=max(max(cost2));
%for i=1:Rc2
%    for j=1:Rc1
%        if cost(i,j)==0
%            cost(i,j)=MH;
%        end
%        if cost2(i,j)==0
%            cost2(i,j)=MH2;
%            costh2(i,j)=NaN;
%        end
%    end
%end
%M=min(min(cost));
%[Ru,Rc]=find(cost==M);
%M2=min(min(cost2));
%[B2, Rbb2]=find(cost2==M2);

%figure 
%mesh(costh2)
%costb=costh2;


%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\costb.mat','costb')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Ev.mat','Ev')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Ea.mat','Ea')

%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\ERR.mat','ER')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Egg.mat','Egbar')

