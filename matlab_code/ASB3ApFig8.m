% ASB3ApFig8 
% appendix Figure 8
% Kenneth I. Carlaw Aug. 2023
% The data in lamrho.xlsx are genereted using ASB3mainconv.m,
% ASB3RoptRCNZ.M and ASB3RoptCLR.M

clear
N=100;
NN=N+1;
X=500000;
F=1;
%lam=2;
PR=2;
PR2=0;
chi1=1;
dim=71;

opts = spreadsheetImportOptions("NumVariables", 10);
opts.Sheet = "lamrho";
opts.DataRange = "A2:J9";
opts.VariableNames = ["lrat", "DR1", "DR2", "DR3", "eff1b", "eff2b", "eff3b","RC1","RC2","RC3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
tbl = readtable("C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\lamrho.xlsx", opts, "UseExcel", false);
lrat = tbl.lrat;
DR1 = tbl.DR1;
DR2 = tbl.DR2;
DR3 = tbl.DR3;
eff1b = tbl.eff1b;
eff2b = tbl.eff2b;
eff3b = tbl.eff3b;
RC1 = tbl.RC1;
RC2 = tbl.RC2;
RC3 = tbl.RC3;
clear opts tbl

%load acn103;
load Rcon103;
load Gbar103;
load mabl103;
load Rcon1032;
load Gbar1032;
load mabl1032;
%load Rcon106;
%load Gbar106;
%load mabl106;
load ERR;
load Egg;
load ERR2;
load Egg2;

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

figure
tile=tiledlayout(1,2);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile
hold on
box on
plot(lrat,RC1*100,'-*','Color','k')
plot(lrat,RC2*100,'-o','Color','k','LineStyle','--')
plot(lrat,RC3*100,'-s','Color','k','LineStyle',':')
%ylabel('% cost saving','FontSize',20)
%xlabel('RAT = $\frac{({\lambda}-1)}{\rho}$','Interpreter','latex','FontSize',20)
xlabel('RAT','FontSize', 20)
title('Panel 1: reserve capacity as a function of RAT','FontSize',20)
xlim([0 10])
ylim([-1.5 100])
ytickformat('percentage')
legend('DROP(0.6,0.2)=0.704','DROP(0.3,0.35)=0.517','DROP(0.75,0.3)=0.355', 'Location','NorthEast','FontSize',12)
legend boxoff
hold off
nexttile
hold on
box on
plot(lrat,DR1*100,'-*','Color','k')
plot(lrat,DR2*100,'-o','Color','k','LineStyle','--')
plot(lrat,DR3*100,'-s','Color','k','LineStyle',':')
%ylabel('% cost saving','FontSize',20)
%xlabel('RAT = $\frac{({\lambda}-1)}{\rho}$','Interpreter','latex','FontSize',20)
xlabel('RAT','FontSize', 20)
title('Panel 2: cost saving as a function of RAT','FontSize',20)
xlim([0 10])
ylim([-1.5 18])
ytickformat('percentage')
legend('DROP(0.6,0.2)=0.704','DROP(0.3,0.35)=0.517','DROP(0.75,0.3)=0.355', 'Location','NorthEast','FontSize',12)
legend boxoff
hold off
%nexttile
%hold on
%box on
%plot(lrat,eff1b*100,'-*','Color','k')
%plot(lrat,eff2b*100,'-o','Color','k','LineStyle','--')
%plot(lrat,eff3b*100,'-s','Color','k','LineStyle',':')
%ylabel('% efficiency','FontSize',20)
%xlabel('RAT = $\frac{({\lambda}-1)}{\rho}$','Interpreter','latex','FontSize',20)
%title('Panel 3: efficiency as a function of RAT','FontSize',20)
%xlim([0 10])
%ylim([-1.5 100])
%ytickformat('percentage')
%legend('DROP(0.6,0.2)=0.704','DROP(0.3,0.35)=0.517','DROP(0.75,0.3)=0.355', 'Location','NorthEast','FontSize',12)
%legend boxoff
%hold off





figure %figure 8 Appendix
tile=tiledlayout(3,4);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile ([2 1])
hold on
box on
plot(RR,ccost(40,:)','Color','k')
plot(RR,llim(40,:)','Color','k','LineStyle',':')
plot(RR,ulim(40,:)','Color','k','LineStyle','--')
title('Panel 1: very low cost of ASB (spitting on the sidewalk)')
ylim([llim(40,1)-10 ulim(40,1)+10])
%xlabel('Deterrence resources (R)')
ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=2','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off
nexttile ([2 2])
hold on
box on
plot(RR,ccost(100,:)','Color','k')
plot(RR,llim(100,:)','Color','k','LineStyle',':')
plot(RR,ulim(100,:)','Color','k','LineStyle','--')
title('Panel 2: intermediate cost of ASB (property crime)')
xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=5','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off
nexttile ([2 1])
hold on
box on
plot(RR,ccost(500000,:)','Color','k')
plot(RR,llim(500000,:)','Color','k','LineStyle',':')
plot(RR,ulim(500000,:)','Color','k','LineStyle','--')
%ylim([5000 25000])
title('Panel 3: very high cost of ASB (murder)')
%xlabel('Deterrence resources (R)')
%ylabel('E(cost|R)')
legend('Expected costs of crime {\rho}=2 and {\lambda}=400','E(cost|R=0)','E(cost|R=N)','Location','Best')
legend boxoff
hold off
nexttile ([1 4])
hold on
box on
plot((lam-1)/2,Rstar,'Color','k')
%plot((lam-1)/2,Rstar2,'*','markerIndices',1:50:length(lam),'LineStyle','--','Color','k')
ylabel('PR^{*}')
xlabel('RAT')
xlim([0 50])
ylim([-5 70])
legend('PR^{*} for baseline parameters','Location','East')
legend boxoff
hold off

for i=1:X
    lam(i)=5; %0.05*i;
    chi(i)=0.02*i;
    if chi(i)<1
        chi(i)=1;
    end

    for r=1:NN
        cost3(r)=PR*RR(r)+(lam(i)-1)*mgbar2(r)+a(r)*chi(i);
        ccost3(i,r)=PR*RR(r)+(lam(i)-1)*mgbar2(r)+a(r)*chi(i);
        llim3(i,r)=cost3(1);
        ulim3(i,r)=cost3(NN);
    end    
end


