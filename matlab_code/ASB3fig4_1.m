%ASBconV2fig31to9.m
%Baseline, convergence, two bin, perisitence, passive R 
%K. I. Carlaw Aug, 2022

clear

%parameters 

N=100;      %population of agents
NN=N+1;
MM=200;

Block=30000;

    %baseline parameterization
F=1;        %saction for ASB
gam=0.8;    % max objective apprehension prob.
aa=1;       %shape poarameter 1 for Beta distribution (Bayesian) 
bb=0.25;    %shape poarameter 2 for Beta distribution (Bayesian)
mu=0.6;     %mean value of g, individual benefit from ASA
sig=0.2;    %varance of g 
lam=5;      %socail cost parameter for ASB 
rho=2;      %unit cost of enforcement resources
Z=2;        %z-history length

BinE=50; % Bin edge passive policies

eps=8;      %parameter for exponenetial, objective apprehension function

critconv1=0.01;
RR=zeros(NN,1);
mv=zeros(NN,1);
mov=zeros(NN,1);
sdv=zeros(NN,1);
coefv=zeros(NN,1);
coefvw=zeros(NN,1);
ma=zeros(NN,1);
sda=zeros(NN,1);
gpiv=zeros(NN,1);
vcon=zeros(NN,Block);OBP=zeros(NN,Block);SBP=zeros(NN,Block);
cost=zeros(NN,1);
cost2=zeros(NN,1);
mgb=zeros(NN,1);mq=zeros(NN,1);mA=zeros(NN,1);
ARATE=zeros(NN,1);
ARATER=NaN(NN,1);
IRATE=zeros(NN,1);
IRATER=NaN(NN,1);
EXCAP=NaN(NN,1);
EI=zeros(NN,1);
elsc=zeros(NN,1);elsu=zeros(NN,1);pc=zeros(NN,1);pu=zeros(NN,1);Tc=zeros(NN,1);Tu=zeros(NN,1);
mvc=zeros(NN,1);mvu=zeros(NN,1);sdvc=zeros(NN,1);sdvu=zeros(NN,1);
SIGGc=zeros(NN,1);SIGGu=zeros(NN,1);
OPERRc=zeros(NN,1);OPERRu=zeros(NN,1);PrEvc=zeros(NN,1);PrEvu=zeros(NN,1);
BH=zeros(NN,1);RRC=zeros(NN,1);NRRC=zeros(NN,1);CMPL=zeros(NN,1);corqp=zeros(NN,1);
hfreq=zeros(NN,NN,1);

edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end
%load gfix.mat;
%g=gfix;
for r=1:NN
    RR(r)=r-1;
   
    c=1;
    b=0;
    cc=0;
    TESTCONV=zeros(MM,1);
    TESTCONV(1)=1;
    while ((cc<1) && (b<MM))
        g=normrnd(mu,sig,Block,N);
        %g=unifrnd(0.2,1.1,Block,N); %
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
        RCw=zeros(Block,1);NRCw=zeros(Block,1);CMPw=zeros(Block,1);
        b=b+1;
        BH(r)=b;
        if b<2
            in=Z+1;
        else
            in=1;
         end
        for t=1:Block            
            R(t)=RR(r);
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
                if eta(t)*F<=g(t,n) 
                    vw(t)=vw(t)+1;
                    gbar(t,n)=g(t,n);
                end
            end
            gbar3(t)=sum(gbar(t,:));
            A(t)=gam*min(1,R(t)/vw(t));
            %A(t)=gam*(1-1/(eps^(R(t)/vw(t))));
            aw(t)=binornd(vw(t),A(t));
            NCw(t)=rho*R(t)+(lam-1)*gbar3(t);
            ARatw(t)=aw(t)/vw(t);
            if R(t)-vw(t) > 0
                RCw(t)=(R(t)-vw(t));
            end
            if vw(t)-R(t)>0
                NRCw(t)=(vw(t)-R(t));
            end
            CMPw(t)=(N-vw(t))/N;
        end
        if b<2
            v=vw;
            a=aw;
            gb=gbar3;
            NC=NCw;
            Ap=A;
            qp=eta;
            ARat=ARatw;
            RC=RCw;
            NRC=NRCw;
            CMP=CMPw;
        else
            vhold=cat(2,v',vw');
            ahold=cat(2,a',aw');
            gbhold=cat(2,gb',gbar3');
            NChold=cat(2,NC',NCw');
            Aphold=cat(2,Ap',A');
            qphold=cat(2,qp',eta');
            ARathold=cat(2,ARat',ARatw');
            RChold=cat(2,RC',RCw');
            NRChold=cat(2,NRC',NRCw');
            CMPhold=cat(2,CMP',CMPw');
            v=vhold';
            a=ahold';
            gb=gbhold';
            NC=NChold';
            Ap=Aphold';
            qp=qphold';
            ARat=ARathold';
            RC=RChold';
            NRC=NRChold';
            CMP=CMPhold';
        end

        vconv=v(1:b*Block);
        if RR(r)==35
            vcon2=vconv;
        end
        vconvlag=v(1:(b-1)*Block);
        freqvconv=histcounts(vconv(:),edges)/((b)*Block);
        %if RR(r)==39
            hfreq(r,:)=freqvconv;
        %end
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
    SGGc=zeros(c*Block,1);SGGu=zeros(c*Block,1);
    vp=zeros(c*Block,1);vcc=zeros(c*Block,1);vu=zeros(c*Block,1);
    TMBC=zeros(2,2);
    TMB=zeros(2,2);
    per=zeros(2,2);        
    for t=1:c*Block
        if t>1000  
            if v(t)<BinE %(r)
                vp(t)=0;
                vcc(t)=v(t);
                vu(t)=-1;
                SGGc(t)=gb(t);                    
            else
                vp(t)=1;
                vcc(t)=-1;
                vu(t)=v(t);
                SGGu(t)=gb(t);                    
            end
            if (vp(t)==0) && (vp(t-1)==0)
                TMBC(1,1)=TMBC(1,1)+1;
            elseif (vp(t)==0) && (vp(t-1)==1)
                TMBC(1,2)=TMBC(1,2)+1;
            elseif (vp(t)==1) && (vp(t-1)==1)
                TMBC(2,2)=TMBC(2,2)+1;
            else
                TMBC(2,1)=TMBC(2,1)+1;
            end
        end
    end
    
    for j=1:2
        for i=1:2
            TMB(j,i)=TMBC(j,i)/sum(TMBC(j,:));
            per(j,i)=1/(1-TMB(j,i));
        end
    end
    elsc(r)=per(1,1);
    elsu(r)=per(2,2);
    pc(r)=TMB(1,1);
    if isnan(pc(r))
        pc(r)=0;
    end
    pu(r)=TMB(2,2);
    if isnan(pu(r))
        pu(r)=0;
    end    
    Tc(r)=(TMBC(1,1)+TMBC(1,2))/(c*Block-1000);
    Tu(r)=(TMBC(2,2)+TMBC(2,1))/(c*Block-1000);
    vcc(vcc==-1)=[];
    vtc=length(vcc);
    vu(vu==-1)=[];
    vtu=length(vu);
    SGGc(SGGc==0)=[];
    stc=length(SGGc);
    SGGu(SGGu==0)=[];
    stu=length(SGGu);
    SIGGc(r)=mean(SGGc);
    SIGGu(r)=mean(SGGu);    
    mvc(r)=mean(vcc(1000:vtc)); %mean(vcc);
    mvu(r)=mean(vu(1000:vtu)); %mean(vu);
    sdvc(r)=std(vcc(1000:vtc));
    sdvu(r)=std(vu(1000:vtu));
    OPERRc(r)=Tc(r)*pc(r)+(1-Tc(r))*(1-pu(r));
    OPERRu(r)=Tu(r)*pu(r)+(1-Tu(r))*(1-pc(r));
    PrEvc(r)=Tc(r)*mvc(r);
    PrEvu(r)=Tu(r)*mvu(r);
       
    vcon(r,:)=v((c-2)*Block+1:(c-1)*Block);
    mv(r)=mean(v);
    mov(r)=mode(v);
    sdv(r)=std(v);
    ma(r)=mean(a);
    sda(r)=std(a);
    coefv(r)=sdv(r)/mv(r);
    mgb(r)=mean(gb);
    mq(r)=mean(qp);
    mA(r)=mean(Ap);
    cost(r)=rho*RR(r)+(lam-1)*mgb(r); 
    cost2(r)=mean(NC);
    OBP(r,:)=Ap((c-2)*Block+1:(c-1)*Block);
    SBP(r,:)=qp((c-2)*Block+1:(c-1)*Block);
    xx=corrcoef(Ap,qp);
    corqp(r)=xx(1,2);

    RRC(r)=sum(RC)/(length(RC)*RR(r));
    NRRC(r)=sum(NRC)/(length(NRC)*mv(r));
    CMPL(r)=sum(CMP)/length(CMP);


    ARATE(r)=ma(r)/mv(r);
   
    IRATE(r)=min(RR(r),mv(r))/mv(r);
    
    if RR(r)>0
        EXCAP(r)=max(0,(RR(r)-mv(r))/RR(r));
        ARATER(r)=ARATE(r)/RR(r);
        IRATER(r)=IRATE(r)/RR(r);
    end
    if RR(r)==41
        vf=v;
        ARATEf=ARat;
        pf=Ap;
        qf=qp;
    end
    if RR(r)==48
        vvff=v;
    end
    
    gp=-1;
    eqn=0;
    while eqn < mv(r)
        gp=gp+0.0001;
        eqn=N*normcdf(gp,mu,sig);
    end
    gpiv(r)=normpdf(gp,mu,sig);
    RR(r)

end
TD=(mv(1)-mv(NN));

RDROP3=zeros(NN,1);RDROP4=zeros(NN,1);RDROP5=zeros(NN,1);RDROP8=zeros(NN,1);
TDROP3=zeros(NN,1);TDROP4=zeros(NN,1);TDROP5=zeros(NN,1);TDROP8=zeros(NN,1);

for r=1:NN
    if r>3
        RDROP3(r)=(mv(r-3)-mv(r))/TD;
        TDROP3(r)=(mv(r-3)-mv(r))/N;
    end
    if r>4
        RDROP4(r)=(mv(r-4)-mv(r))/TD;
        TDROP4(r)=(mv(r-4)-mv(r))/N;
    end
    if r>5
        RDROP5(r)=(mv(r-5)-mv(r))/TD;
        TDROP5(r)=(mv(r-5)-mv(r))/N;
    end
    if r>8
        RDROP8(r)=(mv(r-8)-mv(r))/TD;
        TDROP8(r)=(mv(r-8)-mv(r))/N;
    end
end
MRDROP3=max(RDROP3);
MRDROP4=max(RDROP4);
MRDROP5=max(RDROP5);
MRDROP8=max(RDROP8);

MTDROP3=max(TDROP3);
MTDROP4=max(TDROP4);
MTDROP5=max(TDROP5);
MTDROP8=max(TDROP8);

%MRD=[MRDROP3 MRDROP4 MRDROP5 MRDROP6]



%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\RRC.mat','RRC')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\NRRC.mat','NRRC')
%save('C:\Users\kcarlaw\Documents\MATLAB\ASB2023postRED\CMPL.mat','CMPL')

[M, I]=min(cost)
Rst=RR(I);
cost(1)
clifeqdata=table(RR(38:45),mvc(38:45),sdvc(38:45),mvu(38:45),sdvu(38:45));
[acf38,lag38]=autocorr(vcon(39,Block-10000:Block-5000));
[acf39,lag39]=autocorr(vcon(40,Block-10000:Block-5000));
[acf40,lag40]=autocorr(vcon(41,Block-10000:Block-5000));
[acf41,lag41]=autocorr(vcon(42,Block-10000:Block-5000));
[acf42,lag42]=autocorr(vcon(43,Block-10000:Block-5000));
[acf43,lag43]=autocorr(vcon(44,Block-10000:Block-5000));

figure %figure 4.1
tile=tiledlayout(5,5);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile([1 3])
plot(vcon(40,Block-10000:Block-5000),'Color',[0.5 0.5 0.5])
ylabel('R = 39','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(vcon(40,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 100])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lag39(2:4),acf39(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile([1 3])
plot(vcon(41,Block-10000:Block-5000),'Color',[0.5 0.5 0.5])
ylabel('R = 40','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(vcon(41,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 100])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lag40(2:4),acf40(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile([1 3])
plot(vcon(42,Block-10000:Block-5000),'Color',[0.5 0.5 0.5])
ylabel('R = 41','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(vcon(42,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 100])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lag41(2:4),acf41(2:4),'Filled','Color','k')
xlim([0 4]);
%title('Autocorrelation in violations')
%xlabel('Lags')
nexttile([1 3])
plot(vcon(43,Block-10000:Block-5000),'Color',[0.5 0.5 0.5])
ylabel('R = 42','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(vcon(43,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 100])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lag42(2:4),acf42(2:4),'Filled','Color','k')
xlim([0 4]);
nexttile([1 3])
plot(vcon(44,Block-10000:Block-5000),'Color',[0.5 0.5 0.5])
%title('Panel 6: R = 43')
ylim([0 100])
xlim([0 5000])
ylabel('R = 43','FontSize',20)
xlabel('Period (T = 0,...,5000)','FontSize',20)
nexttile
histogram(vcon(44,Block-10000:Block-5000),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 100])
ylim([0 0.15])
ylabel('Frequency of violations','FontSize',20)
view(90,-90)
nexttile
stem(lag43(2:4),acf43(2:4),'Filled','Color','k')
xlim([0 4]);
xlabel('AC for 3 lags','FontSize',20)

pmv=mv/N;
pRR=RR/N;
%output from baseline
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\movbl mov -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mvbl mv -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pmvbl pmv -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\pRR pRR -ASCII -DOUBLE;

%output for figure 5.1 
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Rcon103 RR -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\Gbar103 mgb -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\mabl103 ma -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TGB Tc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\TBB Tu -ASCII -DOUBLE;

%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDGB3 elsc -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\EDBB3 elsu -ASCII -DOUBLE;


