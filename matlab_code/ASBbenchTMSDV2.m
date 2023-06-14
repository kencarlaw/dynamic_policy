%ASBbenchTMSDV2.m 
%Benchmarking Table A.1 
% Updated Aug 26, 2022, K. I. Carlaw 
% Creates 51 stationary distributions and saves them in SD.mat from the ... 
% transition matrix for the N=50, z=1 benchmarking exercise


%clear
%parameters
N=100;      %population of agents
NN=N+1;
F=1;        %individual cost of apprehnsion given ASA


%apprehension probability APR(v)=gam*(1-a^(-(R/v)))
gam=0.8;
aa=1;
bb=0.25;

T=5000;

%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=0.6;    %mean value of gi, individual benefit from ASA
sig=0.2;   %varance of gi

% Calculating collapsed transition probability matrix (TP) and SD(R) for Markov process
vTP=zeros(NN,1);
aTP=zeros(NN,1);
SD=zeros(NN,NN);
SD2=zeros(NN,NN);
d=zeros(NN,NN);
R=zeros(NN,1);


for r=1:1
    %variables used inside loop structure
    simlength=T;
    TP=zeros(NN,NN);
    aaTP=zeros(NN,NN);
    gTP=zeros(NN,NN);
    etaTP=zeros(NN,NN);
    ATP=zeros(NN,1);
    check=zeros(NN,1);
    x0=zeros(NN,1);
    R(r)=40; %r-1;
    for j=1:NN
        vTP(j)=j-1;
        aTP(j)=j-1;
        x0(j)=0;
    end
x0(1)=1;

    for j=1:NN
        for k=1:NN
            if aTP(k)<=vTP(j)
                    etaTP(j,k)=(aa+aTP(k))/(aa+bb+vTP(j));
                gTP(j,k)=1-normcdf(F*(etaTP(j,k)),mu,sig);
            else
                gTP(j,k)=0;
            end
        end
        ATP(j)=gam*min(1,R(r)/vTP(j));
    end
    for j=1:NN
        for k=1:NN
            aaTP(j,k)=binopdf(aTP(k),vTP(j),ATP(j));
        end
    end
    for j=1:NN
        for k=1:NN
            for q=1:NN
                TP(j,k)=TP(j,k)+aaTP(j,q)*binopdf(vTP(k),N,gTP(j,q));
            end
        end
        check(j)=sum(TP(j,:));
        j
    end
    mc=dtmc(TP);
    TMv=mc.simulate(simlength-1,'X0',x0);
    d=TP^1000;
    SD2(r,:)=d(1,:);
    SD(r,:)=asymptotics(mc);
    %R(r)

end

save('SD.mat','SD');