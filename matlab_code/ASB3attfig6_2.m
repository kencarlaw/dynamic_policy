
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

Ri=[30,57]; % vectore of R values
BinE=35;
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
load vvr2.mat
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
        q=zeros(Block,1);A=zeros(Block,1);
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
            g=normrnd(mu,sig,N);
            if t>1
                if v(t-1)<BinE
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
load vv30.mat
load v30.mat
Pvv30=zeros(NN,NN,1);Evv30=zeros(NN,NN,1);
U30=zeros(NN,NN,1);L30=zeros(NN,NN,1);
count30=zeros(NN,1);
for r=1:RR
    for i=1:NN
        V(r,i)=0;
        for j=1:NN
            Pvv30(i,j)=vv30(i,j)/sum(vv30(i,:));
            if isnan(Pvv30(i,j))
                Pvv30(i,j)=0;
            end
        end
    end
    Evv30(r,:)=Pvv30*vp;
    U30(r,:)=Evv30(r,:)-vp';
    for i=1:NN
        if U30(r,i)>=0
            L30(r,i)=1;
        else
            L30(r,i)=-1;
        end
        if i>1            
            if L30(r,i)-L30(r,i-1)~=0
                count30(i)=i-1;
            end
        end
    end
    count30(count30==0)=[];
    
end 
xvec=zeros(NN,2);
yvec=zeros(NN,2);
for i=1:count(1)
    for j=1:2
        xvec(i,j)=0.04+(i-1)*0.0095 + (j-1)*0.0095;
        yvec(i,j)=0.594;
    end
end
for i=count(1)+1:NN %count(2)
    for j=1:2
        xvec(i,j)=0.04+NN*0.0095 -((i-1-count(1))*0.0095 + (j-1)*0.0095);
        yvec(i,j)=0.594;
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
xvec30=zeros(NN,2);
yvec30=zeros(NN,2);
for i=1:count30(1)
    for j=1:2
        xvec30(i,j)=0.04+(i-1)*0.0095 + (j-1)*0.0095;
        yvec30(i,j)=0.124;
    end
end
for i=count30(1)+1:count30(2)
    for j=1:2
        xvec30(i,j)=0.048+count30(2)*0.0095 -((i-1-count30(1))*0.0095 + (j-1)*0.0095);
        yvec30(i,j)=0.124;
    end
end
for i=count30(2)+1:count30(3)
    for j=1:2
        xvec30(i,j)=0.043+count30(2)*0.0095 + ((i-1-count30(2))*0.0094 + (j-1)*0.0094);
        yvec30(i,j)=0.124;
    end
end
for i=count30(3)+1:NN
    for j=1:2
        xvec30(i,j)=0.040+NN*0.0095 -((i-1-count30(3))*0.0094 + (j-1)*0.0094);
        yvec30(i,j)=0.124;
    end
end


[acf,lag]=autocorr(vw,10);
[acfb,lagb]=autocorr(v30,10);
xx=[count(1),count(1)];xxb=[count30(1),count30(1)];
yy=[0,count(1)];yyb=[0,count30(1)];
xx1=[54,54];xx1b=[count30(2),count30(2)];
yy1=[0,N];yy1b=[0,count30(2)];
xx2=[88,88];xx2b=[count30(3),count30(3)];
yy2=[0,88];yy2b=[0,count30(3)];



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
%xlabel('Violations in current period (v^{t})','FontSize',14)
ylabel('E(v^{t+1}|v^{t})','FontSize',20)
legend('E(v^{t+1}|v^{t})','FontSize',14,'Location','north') 
title('Panel 1: E(v^{t+1}|v^{t}), (BE^{*}, R^{*}_{1}, R^{*}_{2})=(35, 30, 57)','FontSize',25,'FontWeight','normal')
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
plot(vp, Evv30(1,:),'Color','k','LineStyle','--')
%PP=get(gca,'Position');
plot(vp, vp,'Color','k','LineStyle','-')
plot(xxb,yyb,'Color','k','Linestyle',':')
plot(xx1b,yy1b,'Color','k','LineStyle',':');
plot(xx2b,yy2b,'Color','k','Linestyle',':')
for i=1:N/2
    annotation('textarrow',xvec30(i*2,:),yvec30(i*2,:))
end 
xlabel('Violations in current period (v^{t})','FontSize',20)
ylabel('E(v^{t+1}|v^{t})','FontSize',20)
legend('E(v^{t+1}|v^{t})','FontSize',14,'Location','north') 
title('Panel 2: E(v^{t+1}|v^{t}), R = 30','FontSize',25,'FontWeight','normal')
legend boxoff
txt=['v^{t}=',num2str(count30(1)-1,'%2.f')];
text(count30(1)-1,-3,txt);
txt1=['v^{t}=',num2str(count30(2)-1,'%2.f')];
text(count30(2)-1,-3,txt1);
txt2=['v^{t}=',num2str(count30(3)-1,'%2.f')];
text(count30(3)-1,-3,txt2);
ylim([-6 N])
legend('AutoUpdate','off')
xline(N)
yline(0)
yline(N)
hold off

handaxes1=axes('position',[0.07 0.70 0.17 0.18]);
stem(lag(2:4),acf(2:4),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 4])
xticks([0 1 2 3 4])
ylim([0 1])

handaxes2=axes('position',[0.07 0.225 0.17 0.18]);
stem(lagb(2:4),acfb(2:4),'Filled','Color','k')
title('Autocorrelation in violations')
xlabel('Lags')
xlim([0 4])
xticks([0 1 2 3 4])
ylim([0 1])

handaxes3=axes('position',[0.78 0.65 0.17 0.18]);
histogram(vw(1000:(ct-1)*Block),'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequency of violations')
xlabel('v')

handaxes4=axes('position',[0.78 0.17 0.17 0.18]);
histogram(v30,'Normalization','probability','FaceColor',[0.17 0.17 0.17],'BinWidth',1)
title('frequency of violations')
xlabel('v')

