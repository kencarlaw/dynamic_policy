%ASBbenchConvSDV2.m
%Aug 22, 2021, K. I. Carlaw
%Benchmarking of the sim to the N=50, z=1 transition matrix and stationary 
%distribution. Table A.1  
%Uses file SD.mat created by ASBbenchTMSDV2.m

clear
%parameters
N=50;      %population of agents
NN=N+1;
MM=200;
Block=50000;
%T=MM*Block;  %Time/iteration index for the sim and for simlegth in Markov chain
MC=1000;        %Benchmark MC iterations
%Rstar=0;      %resources devoted to apprehension
Rs=[5 17 45];
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;   % Apprehension prob parameter
aa=1;       % Bayesian parameters
bb=0.25;  % Bayesian parameters 
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
critconv1=0.01;
Z=1;

%Dimensions of Varialbles
%R=zeros(T,1);
%av=zeros(T,1);
%A=zeros(T,1); %Aprehnsion probability

%eta=zeros(T,1);   %subjective prob of apprehention

nn=zeros(NN,1);
edges=zeros(NN+1,1);
%TESTCONV=ones(B,MM,1);
%freqc=zeros(MC,NN,1);
%c=ones(B,1);

ATEST=zeros(MC,1);

%mgbar3=zeros(B,1);
%freqvconv=zeros(MM,NN);
%gbar=zeros(T,N);
%gbar3=zeros(T,1);
%b=zeros(B,1);
%cc=zeros(B,1);
BH=zeros(MC,1);

load SD

for j=1:NN+1
    edges(j)=j-1.5;
end
for i=1:NN
    nn(i)=i-1;
end

RR=Rs(1);
for i=1:MC
    %v=zeros(T,1); %number of violators per period
    %v(1)=round(unifrnd(0,N));
    %a=zeros(T,1);   %number of apprehended violators per period
    %a(1)=round(unifrnd(0,v(1)));

   %currentFile = sprintf('g%d.mat',ind);
   %load(currentFile,'g')

    % Sim loop
    c=1;
    b=0;
    cc=0;
    TESTCONV=zeros(MM,1);
    TESTCONV(1)=1;

    while ((cc<1) && (b<MM))
        g=normrnd(mu,sig,Block,N); % unifrnd(0,1,Block,N); %
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
        vw=zeros(Block,1);aw=zeros(Block,1);ARatw=zeros(Block,1);
        b=b+1;
        
        if b<2
            in=Z+1;
        else
            in=1;
        end
        for t=1:Block            
            R(t)=RR;
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
            NCw(t)=rho*R(t)+(lam-1)*gbar3(t);
            ARatw(t)=aw(t)/vw(t);
        end
        if b<2
            v=vw;
            a=aw;
            gb=gbar3;
            NC=NCw;
            Ap=A;
            qp=eta;
            ARat=ARatw;
        else
            vhold=cat(2,v',vw');
            ahold=cat(2,a',aw');
            gbhold=cat(2,gb',gbar3');
            NChold=cat(2,NC',NCw');
            Aphold=cat(2,Ap',A');
            qphold=cat(2,qp',eta');
            ARathold=cat(2,ARat',ARatw');
            v=vhold';
            a=ahold';
            gb=gbhold';
            NC=NChold';
            Ap=Aphold';
            qp=qphold';
            ARat=ARathold';
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
    BH(i)=b;
    %freqc(i,:)=freqvconv;
    ATEST(i)=sum(abs(freqvconv(:)-SD(RR+1,:)'));
    %clear (currentFile,'g')
    %mgbar3(i)=mean(gbar3);
i
end
%figure
%histogram(ATEST)
MATEST=mean(ATEST);
SATEST=std(ATEST);
MaxATEST=max(ATEST);
mnstar=mean(BH);
