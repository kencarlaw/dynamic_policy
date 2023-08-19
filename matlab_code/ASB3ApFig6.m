%ASB3ApFig6.m
%August 14, 2023, Kenneth I. Carlaw   

clear

%parameters 

N=100;      %population of agents
NN=N+1;
%M=20;
MM=200;

Block=50000;

%T=MM*Block;     %Time/iteration index for the sim

Rstar=0;      %resources devoted to apprehension
%F=1;        %individual cost of apprehnsion given ASA
gam=0.80;
aa=1;
bb=0.25;
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=5;     %socail cost conversion of individual ASA for social damage function
rho=2;
init=10;
mm=0;
Z=2;        %z-history length
X=76;

eps=4;

RR=zeros(NN,1);
mvF=zeros(NN,X,1);
%movF=zeros(NN,X,1);
%sdvF=zeros(NN,X,1);
%coefvF=zeros(NN,X,1);
XX=zeros(X,1);
XXr=zeros(X,1);
YY=zeros(NN,1);
YYf=zeros(NN,1);

critconv1=0.01;
F=zeros(X,1);

for r=1:NN
    for x=1:X
        F(x)=(x-1)*0.02;
        XX(x)=x;
        XXr(x)=45;
        YY(r)=r-1;
        YYf(r)=50;        
        RR(r)=r-1;

        edges=zeros(NN+1,1);
        for j=1:NN+1
            edges(j)=j-1.5;
        end
        c=1;
        b=0;
        cc=0;
        TESTCONV=zeros(MM,1);
        TESTCONV(1)=1;
        while ((cc<1) && (b<MM))
            g=normrnd(mu,sig,Block,N);
            vz=zeros(Block,1);
            az=zeros(Block,1);
            q=zeros(Block,1);
            R=zeros(Block,1);
            gbar=zeros(Block,N);
            gbar3=zeros(Block,1);
            A=zeros(Block,1);
            NCw=zeros(Block,1);NCw2=zeros(Block,1);
            vw=zeros(Block,1);aw=zeros(Block,1);
            b=b+1;
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
                q(t)=(aa+az(t))/(aa+bb+vz(t));
                for n=1:N
                    if q(t)*F(x)<=g(t,n) %g(t+ind(r),n)
                        vw(t)=vw(t)+1;
                        gbar(t,n)=g(t,n); %g(t,n);
                    end
                end
                gbar3(t)=sum(gbar(t,:));
                A(t)=gam*min(1,R(t)/vw(t));
                %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                aw(t)=binornd(vw(t),A(t));
                NCw(t)=rho*R(t)+(lam-1)*gbar3(t);
            end
            if b<2
                v=vw;
                a=aw;
                gb=gbar3;
                NC=NCw;
            else
                vhold=cat(2,v',vw');
                ahold=cat(2,a',aw');
                gbhold=cat(2,gb',gbar3');
                NChold=cat(2,NC',NCw');
                v=vhold';
                a=ahold';
                gb=gbhold';
                NC=NChold';
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
        
        mvF(r,x)=mean(v);
%        movF(r,x)=mode(v);
%        sdvF(r,x)=std(v);
%        coefvF(r,x)=sdvF(r,x)/mvF(r,x);.

    end
    
    RR(r)
end


%[M, I]=min(cost);
%Rstar=RR(I);
figure 
mesh(mvF,'EdgeColor',[0.7 0.7 0.7],'FaceColor','none')
hold on
plot3(XX,XXr,mvF(45,XX),'Color','k','LineStyle','-','LineWidth',2)
plot3(YYf,YY,mvF(YY+1,50),'Color','k','LineStyle','--','LineWidth',2)
ylabel('Enforcement resources, R = {0,...,100}', 'fontsize', 12,'color', [0 0 0])
xlabel('Level of sanction, F = {0,...,1.5}', 'fontsize', 12,'color', [0 0 0])
zlabel('E(v|F,R)', 'fontsize', 12,'color', [0 0 0])
legend('E(v|F,R)','R=44','F=1','Location','Best')
xlim([0 76])
ylim([0 100])
xticks([0 13 26 38 51 63 76])
xticklabels({'0.00','0.25','0.50','0.75','1.00','1.25','1.50',})
hold off

