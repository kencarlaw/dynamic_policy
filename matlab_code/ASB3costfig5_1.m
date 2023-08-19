% ASBconV2Fig51Fig72 costs 
% Figure 5.1 and 7.2
% K.I. Carlaw Aug. 2022

clear
N=100;
NN=N+1;
X=5000;
F=1;
%lam=2;
PR=2;
PR2=0;
chi1=1;
dim=71;

load Rcon103;
load Gbar103;
load mabl103;
load POC2.mat;load POC5.mat;load POC50.mat;
load ROC2.mat;load ROC5.mat;load ROC50.mat;
load POC3.mat;

RR=Rcon103;
%ma=acn103;
mgbar2=Gbar103;
a=mabl103;
RR2=Rcon1032;
%ma=acn103;
mgbar22=Gbar1032;

cost=zeros(NN,1);
cost2=zeros(NN,1);
cost12=zeros(NN,1);
cost3=zeros(NN,1);
ccost3=zeros(X,NN,1);
ccost=zeros(X,NN,1);
ccost2=zeros(X,NN,1);
costT=zeros(NN,1);
llim=zeros(X,NN,1);
ulim=zeros(X,NN,1);
llim2=zeros(X,NN,1);
ulim2=zeros(X,NN,1);
llim3=zeros(X,NN,1);
ulim3=zeros(X,NN,1);
lam=zeros(X,1);
chi=zeros(X,1);
M=zeros(X,1);M2=zeros(X,1);MM2=zeros(X,1);
I=zeros(X,1);I2=zeros(X,1);II2=zeros(X,1);
Rstar=zeros(X,1);Rstar2=zeros(X,1);
RRcost=zeros(dim,dim,1);RRM=zeros(X,1);
RRcost2=zeros(dim,dim,1);RRM2=zeros(X,1);

for i=1:X
    lam(i)=0.05*i;

    for r=1:NN
        cost(r)=PR*RR(r)+(lam(i)-1)*mgbar2(r);
        cost2(r)=2*PR*RR(r)+2*(lam(i)-1)*mgbar2(r);
        cost12(r)=PR*RR2(2)+(lam(i)-1)*mgbar22(r);
        ccost(i,r)=PR*RR(r)+(lam(i)-1)*mgbar2(r); %+ma(r)*F;
        ccost2(i,r)=PR2*RR(r)+(lam(i)-1)*mgbar2(r);
        llim(i,r)=cost(1);
        ulim(i,r)=cost(NN);
        llim2(i,r)=cost2(1);
        ulim2(i,r)=cost2(NN);
        
        costT(r)=PR*RR(r)+(10000000000000-1)*mgbar2(r);
        
    end    
    [M(i), I(i)]=min(cost);
    [M2(i),I2(i)]=min(cost2);
    [MM2(i),II2(i)]=min(cost12); 

    Rstar(i)=RR(I(i));
    Rstar2(i)=RR(I2(i));
    for r1=1:dim
        for r2=1:dim
            RRcost(r1,r2)=PR*ER(r1,r2)+(lam(i)-1)*Egbar(r1,r2);
            RRcost2(r1,r2)=PR*ER2(r1,r2)+(lam(i)-1)*Egbar2(r1,r2);
            if RRcost(r1,r2)==0
                RRcost(r1,r2)=NaN;
            end
            if RRcost2(r1,r2)==0
                RRcost2(r1,r2)=NaN;
            end
        end
    end
    RRM(i)=min(min(RRcost));
    RRM2(i)=min(min(RRcost2));
end

[MM,II]=min(costT);
CS=zeros(X,1);
CS2=zeros(X,1);
for i=1:X
    CS(i)=(M(i)-RRM(i))/M(i);
    CS2(i)=(MM2(i)-RRM2(i))/MM2(i);
end


figure %figure 5.1
tile=tiledlayout(1,2);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([1 1])
hold on
box on
plot(RR,ccost(40,:)','Color','k')
plot(RR,POC2,'Color','k','LineStyle','--')
plot(RR,ROC2,'Color','k','LineStyle','-.')
title('Panel 1: Low cost of ASB (spitting on the sidewalk)','FontSize',20)
ylim([llim(40,1)-10 ulim(40,1)+10])
ylabel('E(cost|R)','FontSize',20)
xlabel('Deterrence resources (R)','FontSize',20)
legend('Dynamic model','Pareto Eqn','Random Eqn','Location','best','FontSize',15)
legend boxoff
hold off
nexttile ([1 1])
hold on
box on
plot(RR,ccost(100,:)','Color','k')
plot(RR,POC3,'Color','k','LineStyle','--')
plot(RR,ROC5,'Color','k','LineStyle','-.')
title('Panel 2: Intermediate cost of ASB (property crime)','FontSize',20)
xlabel('Deterrence resources (R)','FontSize',20)
%ylabel('E(cost|R)')
legend('Dynamic model','Pareto Eqn','Random Eqn','Location','best','FontSize',15)
legend boxoff
hold off



