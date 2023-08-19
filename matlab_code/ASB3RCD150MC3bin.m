% ASBcon3R25RCDMC150.m MC 150 on best 25 RCD policies 
% April, 2022, Kenneth I. Carlaw
%optimal Rgb Rtb Rbb BinE BinE2

clear
%parameters 
N=100;      %population of agents
NN=N+1;
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
%X=5;

MC=150;
dim=25;
count=0;

BE1=zeros(dim,1);BE2=zeros(dim,1);RG=zeros(dim,1);RT=zeros(dim,1);RB=zeros(dim,1);
hcost2=zeros(dim,MC,1);cost2=zeros(MC,1);SDcost2=zeros(dim,1);
cost3=zeros(MC,1);MG=zeros(MC,1);MR=zeros(MC,1);Mcost3=zeros(dim,1);
MinC2=zeros(dim,1);MaxC2=zeros(dim,1);Mcost2=zeros(dim,1);
BinE=[30 30 31 32 32 32 32 32 33 33 33 33 34 34 34 34 34 34 34 34 35 35 35 36 36];
RGB=[29 29 29 29 30 30 30 30 30 30 30 31 30 30 30 30 30 30 31 31 30 31 31 31 31];
RTB=[37 44 47 40 42 46 47 50 45 42 47 47 47 47 46 72 84 57 44 43 50 50 55 52 50];
BinE2=[38 42 40 45 37 47 43 67 49 41 50 48 44 38 51 51 56 72 42 48 61 64 51 47 61];
RBB=[59 65 59 64 60 65 62 69 64 62 64 66 64 60 62 73 92 83 57 71 65 77 70 62 65];

for i=1:dim
    %BinE2(i)=BE2+(i-1);
    %for rb=1:dim
        %RBB(rb)=RB2+(rb-1);
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
            count1=0;
            count2=0;
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
                        if vw(t-1)<BinE(i)
                            R(t)=RGB(i);
                        elseif (BinE(i)<=vw(t-1))&&(vw(t-1)<BinE2(i))
                            R(t)=RTB(i);
                            count1=count1+1;
                        elseif BinE2(i)<=vw(t-1)
                            R(t)=RBB(i);
                            count2=count2+1;
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
                    RR=R;
                    GG=gbar3;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    RRhold=cat(2,RR',R');
                    GGhold=cat(2,GG',gbar3');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
                    RR=RRhold';
                    GG=GGhold';
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
            hcost2(count,m)=mean(NC2);
            cost2(m)=mean(NC2);
            MR(m)=mean(RR);
            MG(m)=mean(GG);
            cost3(m)=rho*MR(m)+(lam-1)*MG(m);
            %m
        end
        SDcost2(count,:)=std(cost2);
        MinC2(count)=min(cost2);
        MaxC2(count)=max(cost2);
        Mcost2(count)=mean(cost2);
        Mcost3(count)=mean(cost3);
        BE1(count)=BinE(i);
        BE2(count)=BinE2(i);
        RG(count)=RGB(i);
        RT(count)=RTB(i);
        RB(count)=RBB(i);            
    %end
    i
end
[MC,IC]=min(Mcost2);

OPTRCD1=[BE1 BE2 RG RT RB Mcost2 SDcost2];

OPTRCD=[BE1(IC) RG(IC) RT(IC) BE2(IC) RB(IC) MC];

%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2021\OPTRCD.txt','OPTRCD1','-ascii')