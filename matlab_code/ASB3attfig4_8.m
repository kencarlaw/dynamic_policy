% ASBattractorsblFig34V2PostRED.m
% Uses input from ASBEquilibMDL.m
% Revised Apr 26, 2023, K. I. Carlaw 
clear

%parameters
%T=500000000;    %Time/iteration index for the sim
N=100;      %population of agents
NN=N+1;
%Rstar=41;      %fixed resources devoted to apprehension
F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25; %0.1765;
%critical distribution of ASA (anti-scoial act) RN(gi;mu,sig,0,inf)
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
PR=2;
R=[44]; %[44]; % vector of R values. Toggle
RR=1;
Z=2;
Block=10000000;  %number of sim periods in convergence loop

%BinE=50;
Q=0;QQ=zeros(NN,1);XX=zeros(NN,1);

Evv=zeros(RR,NN,1);Eqq=zeros(RR,NN,1);Eqq2=zeros(RR,NN,1);Eqq3=zeros(RR,NN,1);
U=zeros(RR,NN,1);
V=zeros(RR,NN,1);
L=zeros(RR,NN,1);
count=zeros(NN,1);
vp=zeros(NN,1);qp=ones(NN,1);
ct=0;

for i=1:NN
    vp(i)=i-1;
end
%load vvb.mat
for r=1:RR
    vv=zeros(NN,NN,1);
    qq=zeros(NN,NN,1);qm=zeros(NN,1);
    Pvv=zeros(NN,NN,1);Pqq=zeros(NN,NN,1);
    Q=0;QQ=zeros(NN,1);XX=zeros(NN,1);
    edges=zeros(NN+1,1);
    for j=1:NN+1
        edges(j)=(j-1.5)/N;
    end

    ct=0;
    %if r==1
    %    Block=100000000;
    %elseif r==2
    %    Block=10000000;
    %else
    %    Block=50000000;
    %end
    
%    for i=1:NN
%        vp(i)=i-1;
    while Q < 1
        q=zeros(Block,1);A=zeros(Block,1);
        vz=zeros(Block,1);az=zeros(Block,1);
        v=zeros(Block,1);a=zeros(Block,1);
        for t=1:Block
           if (t<=Z)
                for z=1:Z
                    vz(z)=Z*unifrnd(0,N);
                    az(z)=Z*unifrnd(0,vz(z));
                end
            else
                vz(t)=sum(v(t-Z:t-1));
                az(t)=sum(a(t-Z:t-1));
           end

            q(t)=(aa+az(t))/(aa+bb+vz(t));

            g=normrnd(mu,sig,N);
            for n=1:N      
                if q(t)*F<=g(n)
                    v(t)=v(t)+1;
                end
            end
            %vv(r,t)=v(t);
            A(t)=gam*min(1,R(r)/v(t));
            %A(t)=gam*(1-alpha^(-(R/v(t))));
            a(t)=binornd(v(t),A(t));
        end
        if ct<2
            vw=v;
        else
            vhold=cat(2,vw',v');
            vw=vhold';
        end
        for t=1:Block    
            for i=1:NN
                for j=1:NN
                    if t>3
                        if v(t-1)==i-1
                            if v(t)==j-1
                                vv(i,j)=vv(i,j)+1;
                            end
                        end
                    end
                end
            end
            for i=1:NN
                qm(i)=(i-1)/(N);
                for j=1:NN
                    if t>1
                        %if j<=i
                            if (edges(i) < q(t-1)) && (q(t-1) <= edges(i+1))
                                if (edges(j) < q(t)) && (q(t) <= edges(j+1))
                                        qq(i,j)=qq(i,j)+1;
                                end
                            end  
                        %end
                    end
                end
            end
            
            for i=1:NN
                QQ(i)=sum(vv(:,i))>850;
                XX(i)=sum(vv(:,i));
            end
        end
        ct=ct+1
        if ct>1
            Q=1;%all(QQ,'all');
        end        
    end
    %vp(i)
    %end
        %TPvv=(vv.*vp)/(Block-1);
    for i=1:NN
        V(r,i)=0;
        for j=1:NN
            Pvv(i,j)=vv(i,j)/sum(vv(i,:));
            Pqq(i,j)=qq(i,j)/sum(qq(i,:));
            %Pqm(i,j)
            if isnan(Pvv(i,j))
                Pvv(i,j)=0;
            end
        end
    end
%    for i=1:NN
%        for j=1:NN
%            Pqq(i,j)=qq(i,j)/sum(qq(i,:));
%            if isnan(Pqq(i,j))
%                Pqq(i,j)=0;
%            end
%        end
%    end
    Evv(r,:)=Pvv*vp;
    %PPqq=reshape(Pqq,[],1);
    %PPqq(PPqq==0)=[];
    %PPqm=reshape(qm,[],1);
    %PPqm(PPqm==0)=[];
    Eqq(r,:)=Pqq*qm;
    %for i=1:NN
    %    Eqq3(r,i)=Eqq2(r,i)/sum(Eqq2(r,:));
    %    if i>1
    %        Eqq(r,i)=Eqq(r,i-1)+Eqq3(r,i);
    %    end
    %end
    U(r,:)=Evv(r,:)-vp';
    for i=1:NN
        if U(r,i)>=0
            L(r,i)=1;
        else
            L(r,i)=-1;
        end
        if i>1            
            if L(r,i)-L(r,i-1)~=0
                count(i)=i-1;
            end
        end
    end
    count(count==0)=[];
    
end 
xvec=zeros(NN,2);
yvec=zeros(NN,2);
for i=1:count(1)
    for j=1:2
        xvec(i,j)=0.039+(i-1)*0.0095 + (j-1)*0.0094;
        yvec(i,j)=0.165;
    end
end
for i=count(1)+1:NN %count(2)
    for j=1:2
        xvec(i,j)=0.045+NN*0.0095 -((i-1-count(1))*0.0093 + (j-1)*0.0093);
        yvec(i,j)=0.165;
    end
end
%for i=count(2)+1:count(3)
%    for j=1:2
%        xvec(i,j)=0.045+count(2)*0.0095 + ((i-1-count(2))*0.0094 + (j-1)*0.0094);
%        yvec(i,j)=0.165;
%    end
%end
%for i=count(3)+1:NN
%    for j=1:2
%        xvec(i,j)=0.040+NN*0.0094 -((i-1-count(3))*0.0091 + (j-1)*0.0091);
%        yvec(i,j)=0.165;
%    end
%end

[acf,lag]=autocorr(vw(100:Block),10);
%[pacf,plag]=parcorr(vw(100:(ct-1)*Block),10);
xx=[count(1),count(1)];
yy=[0,count(1)];
%xx1=[count(2),count(2)];
%yy1=[0,count(2)];
%xx2=[count(3),count(3)];
%yy2=[0,count(3)];

figure %figure 4.7
tile=tiledlayout(3,1);
tile.Padding='tight';
tile.TileSpacing='tight';
nexttile ([3 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',20)
plot(vp, Evv(1,:),'Color','k','LineStyle','--')
plot(vp, vp,'Color','k','LineStyle','-')
plot(xx,yy,'Color','k','Linestyle',':')
%plot(xx1,yy1,'Color','k','LineStyle',':');
%plot(xx2,yy2,'Color','k','Linestyle',':')
for i=1:N/2
    annotation('textarrow',xvec(i*2,:),yvec(i*2,:))
end 
xlabel('Violations in current period (v^{t})','FontSize',20)
ylabel('E(v^{t+1}|v^{t})','FontSize',20)
legend('E(v^{t+1}|v^{t})','Location','east','FontSize',14) 
legend boxoff
txt=['v^{t}=',num2str(count(1)-1,'%2.f')];
text(count(1)-1,-3,txt);
%txt1=['v^{t}=',num2str(count(2)-1,'%2.f')];
%text(count(2)-1,-3,txt1);
%txt2=['v^{t}=',num2str(count(3)-1,'%2.f')];
%text(count(3)-1,-3,txt2);
ylim([-6 N])
legend('AutoUpdate','off')
xline(N)
yline(0)
yline(N)
hold off
handaxes2=axes('position',[0.1 0.65 0.20 0.20]);
stem(lag(2:4),acf(2:4),'Filled','Color','k')
xlim([0 4]);
xticks([0 1 2 3 4]);
ylim([-0.5 1]);
title('Autocorrelation in violations','FontSize',14)
xlabel('Lags','FontSize',12)
handaxes3=axes('position',[0.7 0.23 0.20 0.20]);
histogram(vw(100:(ct-1)*Block),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('Frequency of violations','FontSize',14)
xlim([0 100])
xlabel('v')
