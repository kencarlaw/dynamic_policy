% ASB3RV6opt.m 
%August, 2023, Kenneth I. Carlaw

clear

%parameters 

N=100;      %population of agents
NN=N+1;
%M=20;
MM=200;

Block=50000;

F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
SWrho=5*rho;
mm=0;
Z=2;        %z-history length
Disc=0.995;

bine=35;
Rgb=30;
Rbb=57;

dim=71;
Rc1=dim;
Rc2=dim;
bin=dim;
RS1=0;
RS2=0;
critconv1=0.01;

BinE=zeros(bin,1);
R2=zeros(Rc2,1);
R1=zeros(Rc1,1);
%F1=zeros(Rc1,1);
%F2=zeros(Rc2,1);
cost=zeros(Rc2, Rc1,1);
cost2=zeros(Rc2, Rc1,1);
Ev=zeros(Rc2,Rc1,1);
%diagv=zeros(NN,1);
sw=zeros(Rc2,Rc1,1);
Eq=zeros(Rc2,Rc1,1);
diagq=zeros(NN,1);
Ea=zeros(Rc2,Rc1,1);
diaga=zeros(NN,1);
ER=zeros(Rc2,Rc1,1);
Egbar=zeros(Rc2,Rc1,1);
Eswi=zeros(Rc2,Rc1,1);

edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end
for k=1:dim
    BinE(k)=k-1;
    for r2=1:dim %Rc2 %1:Rc2
        R2(r2)=r2-1;
        for r1=Rgb:Rgb %1:r2
            R1(r1)=r1;
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
                        if vw(t-1)<BinE(k)
                            R(t)=R1(r1);
                        else
                            R(t)=R2(r2);
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(R2(r2)-R1(r1)); 
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
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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



            cost(k,r2)=mean(NC);
            cost2(k,r2)=mean(NC2);
            Ev(k,r2)=mean(v);
            Eq(k,r2)=mean(eta);
            Ea(k,r2)=mean(a);
            ER(k,r2)=mean(R);
            Egbar(k,r2)=mean(gb);

        end   
    end
    BinE(k)
end

costh=cost;
costh2=cost2;
MH=max(max(cost));
MH2=max(max(cost2));
for i=1:Rc2
    for j=1:Rc1
        if cost(i,j)==0
            cost(i,j)=MH;
        end
        if cost2(i,j)==0
            cost2(i,j)=MH2;
            costh2(i,j)=NaN;
        end
    end
end
M=min(min(cost));
[Ru,Rc]=find(cost==M);
M2=min(min(cost2));
[B2, Rbb2]=find(cost2==M2);

figure 
mesh(costh2)
costRgb=costh2;

%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2021\costRgb.mat','costRgb')

