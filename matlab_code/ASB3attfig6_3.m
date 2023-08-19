
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
mu=.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %social cost conversion of individual ASA for social damage function
PR=2;

Ri=[23,48]; % vectore of R values
%BinE=35;
Binq=0.745;
RR=1;

Z=2;
Block=100000;  %number of sim periods in convergence loop

Evv=zeros(RR,NN,1);
U=zeros(RR,NN,1);
V=zeros(RR,NN,1);
L=zeros(RR,NN,1);
count=zeros(NN,1);
Evvb=zeros(RR,NN,1);
Ub=zeros(RR,NN,1);
Vb=zeros(RR,NN,1);
Lb=zeros(RR,NN,1);
countb=zeros(NN,1);
vp=zeros(NN,1);
for i=1:NN
    vp(i)=i-1;
end
load vvr2q.mat
%vv=vvr2;
for r=1:RR
    %vv=zeros(NN,NN,1);
    qq=zeros(NN,NN,1);qm=zeros(NN,1);
    Pvv=zeros(NN,NN,1);Pqq=zeros(NN,NN,1);
    Q=0;QQ=zeros(NN,1);XX=zeros(NN,1);
    edges=zeros(NN+1,1);
    for j=1:NN+1
        edges(j)=(j-1.5)/N;
    end

    ct=0;
    while Q < 1
        q=zeros(Block,1);A=zeros(Block,1);qb=zeros(Block,1);
        vz=zeros(Block,1);az=zeros(Block,1);
        v=zeros(Block,1);a=zeros(Block,1);
        R=zeros(Block,1);
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
            %qb(t)=(0.5+(az(t))/(vz(t)));
            g=normrnd(mu,sig,N);
            if t>1
                if q(t)>Binq
                    R(t)=Ri(1);
                else
                    R(t)=Ri(2);
                end
            end
            
            for n=1:N      
                if q(t)*F<=g(n)
                    v(t)=v(t)+1;
                end
            end
            A(t)=gam*min(1,R(t)/v(t));
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
                        if vw(t-1)==i-1
                            if vw(t)==j-1
                                vv(i,j)=vv(i,j)+1;
                            end
                        end
                    end
                end
            end
            %for i=1:NN
            %    qm(i)=(i-1)/(N);
            %    for j=1:NN
            %        if t>1
            %            %if j<=i
            %                if (edges(i) < q(t-1)) && (q(t-1) <= edges(i+1))
            %                    if (edges(j) < q(t)) && (q(t) <= edges(j+1))
            %                            qq(i,j)=qq(i,j)+1;
            %                    end
            %                end  
            %            %end
            %        end
            %    end
            %end
            
            for i=1:NN
                QQ(i)=sum(vv(:,i))>1050;
                XX(i)=sum(vv(:,i));
            end
        end
        Q=all(QQ,'all');
        ct=ct+1
    end
    for i=1:NN
        V(r,i)=0;
        for j=1:NN
            Pvv(i,j)=vv(i,j)/sum(vv(i,:));
            if isnan(Pvv(i,j))
                Pvv(i,j)=0;
            end
        end
    end
    Evv(r,:)=Pvv*vp;
    U(r,:)=Evv(r,:)-vp';
    for i=1:NN
        if U(r,i)>=0
            L(r,i)=1;
        else
            L(r,i)=-1;
        end
        if i>2            
            if L(r,i)-L(r,i-1)~=0
                count(i)=i-1;
            end
        end
    end
    count(count==0)=[];
    
end 
load vv23.mat
load v23.mat
Pvv23=zeros(NN,NN,1);Evv23=zeros(NN,NN,1);
U23=zeros(NN,NN,1);L23=zeros(NN,NN,1);
count23=zeros(NN,1);
for r=1:RR
    for i=1:NN
        V(r,i)=0;
        for j=1:NN
            Pvv23(i,j)=vv23(i,j)/sum(vv23(i,:));
            if isnan(Pvv23(i,j))
                Pvv23(i,j)=0;
            end
        end
    end
    Evv23(r,:)=Pvv23*vp;
    U23(r,:)=Evv23(r,:)-vp';
    for i=1:NN
        if U23(r,i)>=0
            L23(r,i)=1;
        else
            L23(r,i)=-1;
        end
        if i>1            
            if L23(r,i)-L23(r,i-1)~=0
                count23(i)=i-1;
            end
        end
    end
    count23(count23==0)=[];
    
end 
xvec=zeros(NN,2);
yvec=zeros(NN,2);
for i=1:count(1)
    for j=1:2
        xvec(i,j)=0.04+(i-1)*0.0095 + (j-1)*0.0095;
        yvec(i,j)=0.5935;
    end
end
for i=count(1)+1:NN %count(2)
    for j=1:2
        xvec(i,j)=0.04+NN*0.0095 -((i-1-count(1))*0.0095 + (j-1)*0.0095);
        yvec(i,j)=0.5935;
    end
end
%for i=count(2)+1:count(3)
%    for j=1:2
%        xvec(i,j)=0.04+count(2)*0.0095 + ((i-1-count(2))*0.0095 + (j-1)*0.0095);
%        yvec(i,j)=0.5617;
%    end
%end
%for i=count(3)+1:NN
%    for j=1:2
%        xvec(i,j)=0.037+NN*0.0095 -((i-1-count(3))*0.0094 + (j-1)*0.0094);
%        yvec(i,j)=0.5617;
%    end
%end
xvec23=zeros(NN,2);
yvec23=zeros(NN,2);
for i=1:count23(1)
    for j=1:2
        xvec23(i,j)=0.04+(i-1)*0.0095 + (j-1)*0.0095;
        yvec23(i,j)=0.124;
    end
end
for i=count23(1)+1:count23(2)
    for j=1:2
        xvec23(i,j)=0.044+count23(2)*0.0093 -((i-1-count23(1))*0.0093 + (j-1)*0.0093);
        yvec23(i,j)=0.124;
    end
end
for i=count23(2)+1:count23(3)
    for j=1:2
        xvec23(i,j)=0.052+count23(2)*0.0094 + ((i-1-count23(2))*0.0093 + (j-1)*0.0093);
        yvec23(i,j)=0.124;
    end
end
for i=count23(3)+1:NN
    for j=1:2
        xvec23(i,j)=0.055+NN*0.0093 -((i-1-count23(3))*0.0093 + (j-1)*0.0093);
        yvec23(i,j)=0.124;
    end
end


[acf,lag]=autocorr(vw,10);
[acfb,lagb]=autocorr(v23,10);
xx=[count(1),count(1)];xxb=[count23(1),count23(1)];
yy=[0,count(1)];yyb=[0,count23(1)];
xx1=[54,54];xx1b=[count23(2),count23(2)];
yy1=[0,N];yy1b=[0,count23(2)];
xx2=[88,88];xx2b=[count23(3),count23(3)];
yy2=[0,88];yy2b=[0,count23(3)];



figure % Figure6.2
tile=tiledlayout(4,1);
tile.Padding='tight';
tile.TileSpacing='tight';
nexttile ([2 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',20)
%box on
plot(vp, Evv(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xx,yy,'Color','k','Linestyle',':')
%plot(xx1,yy1,'Color','k','LineStyle',':');
%plot(xx2,yy2,'Color','k','Linestyle',':')
%quiver(vp(1:2:NN),V(1,(1:2:NN)),U(1,(1:2:NN)),V(1,(1:2:NN)),'Color',[0.2 0.2 0.2],'AutoScaleFactor',2.5,'LineWidth',0.9,'MaxHeadSize',0.25)
for i=1:N/2
    annotation('textarrow',xvec(i*2,:),yvec(i*2,:))
end 
%handaxes=axes('position',[0.15 0.52 0.3 0.3]);
%stem(lag,acf,'Filled','Color','k')
%xlabel('Violations in current period (v^{t})','FontSize',20)
ylabel('E(v^{t+1}|v^{t})','FontSize',20)
legend('E(v^{t+1}|v^{t})','FontSize',12,'Location','north','FontSize',14) 
title('Panel 1: E(v^{t+1}|v^{t}), (BEq^{*}, R^{*}_{1}, R^{*}_{2})=(0.745, 23,47)','FontSize',25,'FontWeight','normal')
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
nexttile ([2 1])
hold on
set(gca,'TickDir','out','TickLength',[0.005,0.005],'FontSize',20)
%box on
plot(vp, Evv23(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xxb,yyb,'Color','k','Linestyle',':')
plot(xx1b,yy1b,'Color','k','LineStyle',':');
plot(xx2b,yy2b,'Color','k','Linestyle',':')
for i=1:N/2
    annotation('textarrow',xvec23(i*2,:),yvec23(i*2,:))
end 
xlabel('Violations in current period (v^{t})','FontSize',20)
ylabel('E(v^{t+1}|v^{t})','FontSize',20)
legend('E(v^{t+1}|v^{t})','FontSize',12,'Location','north') 
title('Panel 2: E(v^{t+1}|v^{t}), R = 23','FontSize',25,'FontWeight','normal')
legend boxoff
txt=['v^{t}=',num2str(count23(1)-1,'%2.f')];
text(count23(1)-1,-3,txt);
txt1=['v^{t}=',num2str(count23(2)-1,'%2.f')];
text(count23(2)-1,-3,txt1);
txt2=['v^{t}=',num2str(count23(3)-1,'%2.f')];
text(count23(3)-1,-3,txt2);
ylim([-6 N])
legend('AutoUpdate','off')
xline(N)
yline(0)
yline(N)
hold off

handaxes1=axes('position',[0.07 0.70 0.17 0.18]);
box off
set(gca,'TickDir','out','TickLength',[0.005,0.005])
stem(lag(2:4),acf(2:4),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 4])
xticks([0 1 2 3 4])
ylim([0 1])

handaxes2=axes('position',[0.07 0.225 0.17 0.18]);
box off
set(gca,'TickDir','out','TickLength',[0.005,0.005])
stem(lagb(2:4),acfb(2:4),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 4])
xticks([0 1 2 3 4])
ylim([0 1])

handaxes3=axes('position',[0.78 0.65 0.17 0.18]);
box off
set(gca,'TickDir','out','TickLength',[0.005,0.005])
histogram(vw(1000:Block),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequency of violations')
xlabel('v')

handaxes4=axes('position',[0.78 0.17 0.17 0.18]);
box off
set(gca,'TickDir','out','TickLength',[0.005,0.005])
histogram(v23,'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequency of violations')
xlabel('v')

