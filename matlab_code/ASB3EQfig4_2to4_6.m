% ASBNatEQoptcost.m
% Nat equilibrium optimal cost policy comparison codebase
% Updated Nov. 08, 2022
clear
%parameters
N=100;      %population of agents
NN=N+1;
F=1;        %individual cost of apprehnsion given ASA

%apprehension probability APR(v)=gam*(1-a^(-(R/v)))
gam=0.8;
aa=1;
bb=0.25;

%T=100000;

%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=0.6;    %mean value of gi, individual benefit from ASA
sig=0.2;   %varance of gi
rho=2;
lam=5;
Block=20000;

eps=8;

g=normrnd(mu,sig,Block,N);
%gfix=g;
Ec=zeros(Block,1);Ev=zeros(Block,1);ER=zeros(Block,1);I=zeros(Block,1);
MinC=zeros(Block,1);minCost=zeros(Block,1);minCost2=zeros(Block,1);

cost=zeros(Block,NN,1);cost2=zeros(Block,NN,1);
vh=zeros(NN);Rh=zeros(NN);

v=zeros(NN,1);P=zeros(NN,NN,1);P2=zeros(NN,NN,1);
RR=zeros(NN,1);
eqs=zeros(Block,NN,1);leqs=zeros(Block,NN,1);heqs=zeros(Block,NN,1);
ecn3=zeros(Block,NN,1);eqv3=zeros(Block,NN,1);eqvx=zeros(Block,NN,5);
RC3=zeros(Block,NN,1);RCPO=zeros(Block,NN,1);eqv4=zeros(Block,NN,1);
RC4=zeros(Block,NN,1);RCPO2=zeros(Block,NN,1);test=zeros(Block,NN,1);
CMP=zeros(Block,NN,1);CMP2=zeros(Block,NN,1);mCMPL=zeros(NN,1);mCMPL2=zeros(NN,1);
Ecost=zeros(NN,1);meqs=zeros(NN,1);lmeqs=zeros(NN,1);hmeqs=zeros(NN,1);mecn2=zeros(NN,NN,1);
Ecost2=zeros(NN,1);Ecost3=zeros(NN,1);Ev3=zeros(NN,1);Ev4=zeros(NN,1);
mRC3=zeros(NN,1);mRCPO=zeros(NN,1);
mRC4=zeros(NN,1);mRCPO2=zeros(NN,1);
PD=makedist('Normal','mu',mu,'sigma',sig);

edges=zeros(NN+1,1);
for j=1:NN+1
    edges(j)=j-1.5;
end

for r=1:NN
    RR(r)=r-1;
   for j=1:NN
        v(j)=j-1;
        P(r,j)=gam*min(1,(RR(r)/v(j)));
        P2(r,j)=gam*min(1,(RR(r)/(v(j)+1))); 
   end
end

gt=zeros(1,N,1);
for i=1:N
    gt(1,i)=i*.01;
end

for t=1:Block
    ec0=zeros(NN,1);ecN=zeros(NN,N,1);ecN2=zeros(NN,N,1);mecN=zeros(NN,1);
    vn=zeros(NN,N,1);vnn=zeros(NN,N,1);eqv=zeros(NN,N,1);leqv=zeros(NN,N,1);heqv=zeros(NN,N,1);
    ecn=zeros(NN,N,1);mecn=zeros(NN,1);ecn2=zeros(NN,N,1);
    testc=zeros(NN,1);HR=zeros(NN,N,1);In=zeros(NN,1);IN=zeros(NN,1);Inn=zeros(NN,1);
    HHR=zeros(NN,1);ecost=zeros(NN,N,1);ecost2=zeros(NN,N,1);

    sg=sort(g(t,:),'descend');%g(t,:)
    
    for r=1:NN
        %RR(r)=r-1;
            if sg(1)/(P(r,1))<F
                ec0(r)=rho*RR(r)+(lam-1)*sg(1);
            end

        for n=1:N
            if n==N
                if sg(N)>=P(r,N+1)*F
                    HR(r,n)=RR(r);
                    ecN(r,n)=rho*RR(r)+(lam-1)*sum(sg(1:N));
                    ecN2(r,n)=ecN(r,n);
                    eqv(r,n)=v(n+1);
                    heqv(r,n)=eqv(r,n);
                end
            end
            if n<N
                if (sg(n)>=P(r,n)*F) %&& (sg(n)<=gam*min(1,RR(r)/(n+1)))
                    vn(r,n)=v(n+1);
                end
                if (sg(n+1)<P2(r,n)*F)
                    vnn(r,n)=v(n+1);
                end
                if (vn(r,n)==vnn(r,n)) && (vn(r,n)~=0)
                    HR(r,n)=RR(r);
                    eqv(r,n)=vn(r,n);
                    if eqv(r,n)>=70
                        heqv(r,n)=eqv(r,n);
                    elseif eqv(r,n)<=30
                        leqv(r,n)=eqv(r,n);
                    end
                    ecn(r,n)=rho*RR(r)+(lam-1)*sum(sg(1:vn(r,n)));
                    ecn2(r,n)=ecn(r,n);
                end
            end
            if ecn(r,n)==0
                ecn(r,n)=1000;
            end
            if ecN(r,n)==0
                ecN(r,n)=1000;
            end
            ecost(r,n)=min(ecn(r,n),ecN(r,n));
            if ecn2(r,n)==0
                ecost2(r,n)=ecn2(r,n)+ecN2(r,n);
            else
                ecost2(r,n)=ecn2(r,n);
            end
        end
        eqv2=eqv(r,:);
        eqv2(eqv2==0)=[];
        for k=1:length(eqv2)
            eqvx(t,r,k)=eqv2(k);
        end
        test(t,r)=sum(eqvx(t,r,:)>0);
        ind=randi([1 length(eqv2)]);
        eqv3(t,r)=eqv2(ind);
        eqv4(t,r)=min(eqv2);
        freqv3(r,:)=histcounts(eqv3(:,r)',edges)/(Block);
        ecn3(t,r)=rho*RR(r)+(lam-1)*sum(sg(1:eqv3(t,r)));
        cost(t,r)=rho*RR(r)+(lam-1)*sum(sg(1:eqv4(t,r)));
        if RR(r)-eqv3(t,r)>0
            RC3(t,r)=(RR(r)-eqv3(t,r));
        end
        if RR(r)-eqv4(t,r)>0
            RCPO(t,r)=(RR(r)-eqv4(t,r));
        end
        if eqv3(t,r)-RR(r)>0
            RC4(t,r)=(eqv3(t,r)-RR(r));
        end
        if eqv4(t,r)-RR(r)>0
            RCPO2(t,r)=(eqv4(t,r)-RR(r));
        end
        CMP(t,r)=(N-eqv3(t,r))/N;
        CMP2(t,r)=(N-eqv4(t,r))/N;
        XX=ecost2(r,:);
        x=min(XX(XX>0));
        mecn2(t,r)=x;
        eqs(t,r)=sum(ecn2(r,:)>0)+sum(ecN2(r,:)>0);
        leqs(t,r)=sum(leqv(r,:)>0);
        heqs(t,r)=sum(heqv(r,:)>0);
        if eqs(t,r)>0
            [mecn(r),In(r)]=min(ecn(r,:));
            [mecN(r),IN(r)]=min(ecN(r,:));
        else
            mecn(r)=rho*RR(r);
            mecN(r)=rho*RR(N);
        end
        Inn=min(max(IN(r),In(r)));
        %HHR(r)=HR(r,max(IN(r),In(r)));

        if mecn(r)<=mecN(r)
            cost2(t,r)=mecn(r);
        else
            cost2(t,r)=mecN(r);
        end
        
    end
%    [Ec(t),I(t)]=min(cost(t,:));
    eec=min(min(ecost));
    [rr,nn]=find(ecost==eec);
    MinC(t)=eec;
    Ev(t)=v(nn);
    ER(t)=RR(rr);
    minCost(t)=min(cost(t,:));
    minCost2(t)=min(ecn3(t,:));
end
for r=1:NN
    Ecost(r)=mean(cost(:,r));
    meqs(r)=mean(eqs(:,r));
    lmeqs(r)=mean(leqs(:,r));
    hmeqs(r)=mean(heqs(:,r));
    Ecost2(r)=mean(mecn2(:,r));
    Ecost3(r)=mean(ecn3(:,r));
    Ev3(r)=mean(eqv3(:,r));
    Ev4(r)=mean(eqv4(:,r));
    mRC3(r)=mean(RC3(:,r))/RR(r);
    mRCPO(r)=mean(RCPO(:,r))/RR(r);
    mCMPL(r)=mean(CMP(:,r));
    mCMPL2(r)=mean(CMP2(:,r));
    mRC4(r)=mean(RC4(:,r))/Ev3(r);
    mRCPO2(r)=mean(RCPO2(:,r)/Ev4(r));
end
MEc=min(Ecost);
EEv=mean(Ev);
EER=mean(ER);
[acfR10,lagR10]=autocorr(eqv3(:,11));
[acfR20,lagR20]=autocorr(eqv3(:,21));
[acfR30,lagR30]=autocorr(eqv3(:,31));
[acfR40,lagR40]=autocorr(eqv3(:,41));
[acfR50,lagR50]=autocorr(eqv3(:,51));

[acfP10,lagP10]=autocorr(eqv4(:,11));
[acfP20,lagP20]=autocorr(eqv4(:,21));
[acfP30,lagP30]=autocorr(eqv4(:,31));
[acfP40,lagP40]=autocorr(eqv4(:,41));
[acfP50,lagP50]=autocorr(eqv4(:,51));

figure % Figure 4.2
tile=tiledlayout(2,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([1 1])
plot(RR,meqs,'Color','k')
title('Panel 1: Expected number of total equlibria as a function of R','FontSize',20)
ylim([0.75 2.5])
nexttile([1 1])
hold on
plot(RR,lmeqs,'Color','k','LineStyle','-')
plot(RR,hmeqs,'Color','k','LineStyle','--')
title('Panel 2: Expected number of good and bad equlimbria as a function of R','FontSize',20)
xlabel('Enforcement resources (R)','FontSize',20)
legend('# of equlibria {\leq} 30','# of equilibria {\geq} 70', 'Location','east','FontSize',15)
legend box off
ylim([0 1.25])
hold off

figure %figure 4.3
tile=tiledlayout(5,5);
tile.Padding='tight';
tile.TileSpacing='tight';
nexttile([1 3])
plot(eqv4(Block-10000:Block-5000,6),'Color',[0.5 0.5 0.5])
ylabel('R = 10','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv4(Block-10000:Block-5000,6),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagP10(2:4),acfP10(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4]);
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv4(Block-10000:Block-5000,13),'Color',[0.5 0.5 0.5])
ylim([0 100])
xlim([0 5000])
ylabel('R = 20','FontSize',20)
nexttile
histogram(eqv4(Block-10000:Block-5000,13),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagP20(2:4),acfP20(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv4(Block-10000:Block-5000,20),'Color',[0.5 0.5 0.5])
ylabel('R = 30','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv4(Block-10000:Block-5000,20),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagP30(2:4),acfP30(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv4(Block-10000:Block-5000,41),'Color',[0.5 0.5 0.5])
title('Panel 4: R = 26')
ylabel('R = 40','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv4(Block-10000:Block-5000,41),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagP40(2:4),acfP40(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv4(Block-10000:Block-5000,51),'Color',[0.5 0.5 0.5])
xlabel('Tic (T = 0,...,5000)','FontSize',20)
ylabel('R = 50','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv4(Block-10000:Block-5000,51),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
ylabel('Frequency of violations','FontSize',20)
view(90,-90)
nexttile
stem(lagP50(2:4),acfP50(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4]);
xticks([0 1 2 3 4]);
xlabel('AC for 3 lags','FontSize',20)

figure %figure 4.4
tile=tiledlayout(5,5);
tile.Padding='tight';
tile.TileSpacing='tight';
nexttile([1 3])
plot(eqv3(Block-10000:Block-5000,11),'Color',[0.5 0.5 0.5])
ylabel('R=10','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv3(Block-10000:Block-5000,11),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagR10(2:4),acfR10(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4]);
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv3(Block-10000:Block-5000,21),'Color',[0.5 0.5 0.5])
ylabel('R = 20','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv3(Block-10000:Block-5000,21),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagR20(2:4),acfR20(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv3(Block-10000:Block-5000,46),'Color',[0.5 0.5 0.5])
ylabel('R = 30','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv3(Block-10000:Block-5000,46),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagR30(2:4),acfR30(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv3(Block-10000:Block-5000,51),'Color',[0.5 0.5 0.5])
ylabel('R = 40','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv3(Block-10000:Block-5000,51),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
view(90,-90)
nexttile
stem(lagR40(2:4),acfR40(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4])
xticks([0 1 2 3 4]);
nexttile([1 3])
plot(eqv3(Block-10000:Block-5000,66),'Color',[0.5 0.5 0.5])
xlabel('Tic (T = 0,...,5000)','FontSize',20)
ylabel('R = 50','FontSize',20)
ylim([0 100])
xlim([0 5000])
nexttile
histogram(eqv3(Block-10000:Block-5000,66),'Normalization','probability','FaceColor',[0.2 0.2 0.2],'BinWidth',1)
xlim([0 102])
ylim([0 0.15])
ylabel('Frequency of violations','FontSize',20)
view(90,-90)
nexttile
stem(lagR50(2:4),acfR50(2:4),'Filled','Color','k')
ylim([-0.02 1]);
xlim([0 4]);
xticks([0 1 2 3 4]);
xlabel('AC for 3 lags','FontSize',20)

%figure % Figure
%tile=tiledlayout(3,1);
%tile.Padding='none';
%tile.TileSpacing='tight';
%nexttile ([2 1])
%plot(RR,Ecost3,'Color','k')
%xlabel('R')
%ylabel('E(cost|R')
%title('Equlimbirm lowest cost given R')
%nexttile([1 1])
%plot(RR,meqs,'Color','k')
%xlabel('R')
%ylabel('Number of equilibria')
%ylim([0.5 2.5])

load mvbl

figure %figure 4.5 
hold on
box on
plot(RR,mvbl,'Color','k')
plot(RR,Ev3,'Color','k','LineStyle','--','LineWidth',1.0)
plot(RR,Ev4,'Color','k','LineStyle',':','LineWidth',1.0)
legend('E(v|R) D-model','E(v|R) R-model','E(v|R) P-model','FontSize',15)
legend box off
xlabel('Enforcement resources (R)','FontSize',20)
ylabel('Number of violations','FontSize',20)
hold off

load RRC.mat
load NRRC.mat
load CMPL.mat

figure % Figure4.6
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([1 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005])
box off
plot(RR,CMPL,'Color','k');
plot(RR,mCMPL,'Color','k','LineStyle','-.')
plot(RR,mCMPL2,'Color','k','LineStyle','--')
%xlabel('R')
%ylabel('Reserve capacity')
title('Panel 1: Compliance rate for different levels of deterrence resources','FontSize',17,'FontWeight','normal')
legend('D-model','R-model','P-model','Location','best','FontSize',14)
legend box off
ylim([-0.05 1]);
legend('AutoUpdate','off')
xline(100)
yline(1)
hold off
nexttile([1 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005])
box off
plot(RR,RRC,'Color','k')
plot(RR,mRC3,'Color','k','LineStyle','-.')
plot(RR,mRCPO,'Color','k','LineStyle','--')
%xlabel('R')
%ylabel('Capacity','FontSize',20)
title('Panel 2: Reserve capacity rate for different levels of deterrence resources','FontSize',17,'FontWeight','normal')
legend('D-model','R-model','P-model','Location','best','FontSize',14)
legend box off
ylim([-0.05 1]);
legend('AutoUpdate','off')
xline(100)
yline(1)
hold off
nexttile([1 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005])
box off
plot(RR,NRRC,'Color','k')
plot(RR,mRC4,'Color','k','LineStyle','-.')
plot(RR,mRCPO2,'Color','k','LineStyle','--')
xlabel('Enforcement resources (R)','FontSize',14)
%ylabel('Reserve capacity')
title('Panel 3: Deficient capacity rate for different levels of deterrence resources','FontSize',17,'FontWeight','normal')
legend('D-model','R-model','P-model','Location','best','FontSize',14)
legend box off
ylim([-0.05 1]);
yticks([0 0.2 0.4 0.6 0.8 1])
yticklabels({'0','-0.2','-0.4','-0.6','-0.8','-1'})
legend('AutoUpdate','off')
xline(100)
yline(1)
hold off
%ylim([0.5 2.5])


