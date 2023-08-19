% ASBcon2RCDMC150.m 105 MC simulation for the optimal neighborhood of CD, Table 4.1 
% April, 2022, Kenneth I. Carlaw
% MC 150 optimal Rgb Rbb BinE

clear
%parameters 
N=100;      %population of agents
NN=N+1;
MM=200;

Block=20000;

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
%X=5;

MC=150;
count=0;
bd=10;
BinE1=zeros(bd,1);
BinE=30;
rd1=4;
RGB1=zeros(rd1,1);
RGB=29;
rd2=10;RBB1=zeros(rd2,1);
RBB=53;
BE=zeros(bd*rd1*rd2,1);RG=zeros(bd*rd1*rd2,1);RB=zeros(bd*rd1*rd2,1);
cost=zeros(MC,1);
EvMC=zeros(MC,1);SvMC=zeros(MC,1);
SDcost=zeros(bd*rd1*rd2,1);MinC=zeros(bd*rd1*rd2,1);MaxC=zeros(bd*rd1*rd2,1);Mcost=zeros(bd*rd1*rd2,1);
hcost=zeros(bd*rd1*rd2,MC,1);surfC=zeros(bd,rd1*rd2);

for i=1:bd
    BinE1(i)=i+(BinE-1);
    for rg=1:rd1
        RGB1(rg)=rg+(RGB-1);
        for rb=1:rd2
            RBB1(rb)=rb+(RBB-1);
            count=count+1;           
            for m=1:MC

                critconv1=0.01;

                edges=zeros(NN+1,1);
                for j=1:NN+1
                    edges(j)=j-1.5;
                end

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
                    q=zeros(Block,1);
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
                            if vw(t-1)<BinE1(i)
                                R(t)=RGB;
                            else
                                R(t)=RBB1(rb);
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
                        q(t)=(aa+az(t))/(aa+bb+vz(t));
                        for n=1:N
                            if q(t)*F<=g(t,n) %g(t+ind(r),n)
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
            cost(m)=mean(NC2);
            EvMC(m)=mean(v);
            SvMC(m)=std(v);
            hcost(count,m)=cost(m);
            %m
            end
            MinC(count)=min(cost);
            MaxC(count)=max(cost);
            Mcost(count)=mean(cost);
            surfC(i,rb)=mean(cost);
            SDcost(count)=std(cost);
            BE(count)=BinE1(i);
            RG(count)=RGB(rg);
            RB(count)=RBB1(rb);
            
%            MinC(2*(i-1)*rd1+(rg-1)*rd2+rb)=min(cost);
%            MaxC(2*(i-1)*rd1+(rg-1)*rd2+rb)=max(cost);
%            Mcost(2*(i-1)*rd1+(rg-1)*rd2+rb)=mean(cost);
%            SDcost((i-1)*rd1+(rg-1)*rd2+rb)=std(cost);
%            BE(2*(i-1)*rd1+(rg-1)*rd2+rb)=BinE1(i);
%            RG(2*(i-1)*rd1+(rg-1)*rd2+rb)=RGB1(rg);
%            RB(2*(i-1)*rd1+(rg-1)*rd2+rb)=RBB1(rb);

        end
    end
    i
end
[MC,IC]=min(Mcost);

%OPTCD=[BE(IC) RG(IC) RB(IC) MC];

OPTCD1=[BE RG RB Mcost SDcost];
            
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\OPTCDMC150.txt','OPTCD1','-ascii')