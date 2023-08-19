clear

%parameters
T=20000;    %Time/iteration index for the sim
N=100;      %population of agents
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
PR=2;
BinE=53;
Rc=39;
Ru=60;
Rstar=Rc;
Z=2;

%Dimensions of Varialbles

%for r=1:RR
g=normrnd(mu,sig,T,N);
R=zeros(T,1);
R(1)=Rstar;
R(2)=Rstar;
R(3)=Rstar;
R(4)=Rstar;
v=zeros(T,1);   %number of violators per period
vz=zeros(T,1);
a=zeros(T,1);   %number of apprehended violators per period
az=zeros(T,1);
A=zeros(T,1); %Aprehnsion probability
q=zeros(T,1);   %subjective prob of apprehention
bin=zeros(T,1);ARATE=zeros(T,1);

% Sim loop
    for t=1:T
        if t>1
            if v(t-1)<BinE
                R(t)=Rc;
            else
                R(t)=Ru;
            end
        end
        if (t<=Z)
            for z=1:Z
                vz(z)=Z*unifrnd(0,N*2);
                az(z)=Z*unifrnd(0,vz(z));
            end
        else
            vz(t)=sum(v(t-Z:t-1));
            az(t)=sum(a(t-Z:t-1));
        end
        q(t)=(aa+az(t))/(aa+bb+vz(t));
        for n=1:N
            if q(t)*F<=g(t,n)
                v(t)=v(t)+1;
            end
        end
        A(t)=gam*min(1,R(t)/v(t));
        %A(t)=gam*(1-alpha^(-(R/v(t))));
        a(t)=binornd(v(t),A(t));
        bin(t)=BinE;
        ARATE(t)=a(t)/v(t);
        if isnan(ARATE(t))
            ARATE(t)=0;
        end
    end
%figure 
%histogram(v)
%figure 
%plot([v a])

figure 
tile=tiledlayout(3,1);
tile.Padding='none';
tile.TileSpacing='tight';
nexttile
hold on
box on
colororder({'k','k'})
yyaxis left
%plot(v,'+','MarkerSize',4,'MarkerIndices',1:5:T,'Color','k','LineStyle',"-",'LineWidth',0.5)
plot(v,'Color','k')
plot(bin,'Color','k','LineStyle','-.');
ylim([0 N])
xlim([1000 1500])
ylabel('Number of violations','Color','k','FontSize',20)
yyaxis right
%plot(q,'*','MarkerSize',4,'MarkerIndices',1:5:T,'Color','k','LineStyle',"-.",'LineWidth',0.5)
plot(ARATE,'Color','k','LineStyle',"--")
ylim([0 1])
xlim([1000 1500])
ylabel('Apprehension rate','Color','k')
xline(2169,'Color','k','LineStyle',':');
xline(2174,'Color','k','LineStyle',':');
%xline(2221,'Color','k','LineStyle',':');
%xline(2227,'Color','k','LineStyle',':');
legend('violations (v)','Bin edge = 39','apprehension rate (a/v)','FontSize',15)
title('Panel 1: Illustrative path of violation and apprehensions','FontSize',20)
hold off
nexttile
hold on
box on
plot(R,'Color','k','LineWidth',1)
ylim([30 65])
xlim([1000 1500])
xline(2169,'Color','k','LineStyle',':');
xline(2174,'Color','k','LineStyle',':');
%xline(2221,'Color','k','LineStyle',':');
%xline(2227,'Color','k','LineStyle',':');
ylabel('Deterrence resources (R)')
txt1='\leftarrow Deployment of crackdown resources';
txt2='Deployment of crackdown resources \rightarrow';
text(2172,40,txt2);
%text(2180,30,txt2);
title('Panel 2: Illustrative path of deterrence resource (R) deployment','FontSize',20)
hold off 
nexttile
hold on
box on
plot(A,'Color','k')
plot(q,'Color','k','LineStyle','--')
xlim([1000 1500])
legend('Objective probability of apprehension','Subjective probability of apprehension','FontSize',15)
xlabel('Period','FontSize',20)
title('Panel 3: Illustrative path of objective and subjective probabilities','FontSize',20)
hold off



%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61v v -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61bin bin -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61ARATE ARATE -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61R R -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61A A -ASCII -DOUBLE;
%save C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\fig61q q -ASCII -DOUBLE;
