% ASBcon3R1000sim.m 3 Bin RCD lowest cost search 
%April, 2022
%optimal BE BE2 Rgb Rtb Rbb

clear
%parameters 
N=100;      %population of agents
NN=N+1;
MM=200;

Block=20000;

F=1;        %individual cost of apprehnsion given ASA
gam=0.8;
aa=1;
bb=0.25;        %0.1765;
mu=0.6;     %mean value of gi, individual benefit from ASA
sig=0.2;    %varance of gi 
lam=21;     %social cost conversion of individual ASA for social damage function
rho=2;
SWrho=5*rho;
mm=0;
Z=2;        %z-history length
X=5;

MC=50; 
BnE1=zeros(MC,1);BnE2=zeros(MC,1);RRRgb=zeros(MC,1);RRRtb=zeros(MC,1);RRRbb=zeros(MC,1);CM=zeros(MC,1);
count=zeros(MC,1);

for m=1:MC
    %from ASBconvperR2V3.m
    Rgb=round(unifrnd(0,NN));
    Rtb=round(unifrnd(Rgb,NN));
    Rbb=round(unifrnd(Rtb,NN));
    BE=round(unifrnd(0,NN));
    BE2=round(unifrnd(BE,NN));

    crit=1;
    crit2=138.3939;
    RbbM=138.3939+1;

    critconv1=0.01;

    count3=0;
    RR=zeros(NN,1);
    for r=1:NN
        RR(r)=r-1;
    end

    edges=zeros(NN+1,1);
    for j=1:NN+1
        edges(j)=j-1.5;
    end

    while crit>0.02 %RbbM>=crtit2
        costold=RbbM;

        %Search for Bin edge 1
        BinE=zeros(N,1);
        Bcost=zeros(N,1);
        if BE-X<1
            BEd=1;
        else
            BEd=BE-X;
        end
        if BE+X>=N
            BEu=N;
        else
            BEu=BE+X;
        end      
        if BE2<=BEu
            BEu=BE2-1;
        end
        if BEu<BEd
            BEu=BEd;
        end
        for i=BEd:BEu
            BinE(i)=i;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
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
                    if t>1
                        if vw(t-1)<BinE(i)
                            R(t)=Rgb;
                        elseif (BinE(i)<=vw(t-1))&&(vw(t-1)<BE2)
                            R(t)=Rtb;
                        elseif BE2<=vw(t-1)
                            R(t)=Rbb;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(Rbb-Rgb); 
                            swic(t)=1;
                        end                    
                    end                
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
                        if q(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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
            Bcost(i)=mean(NC2);
            %BinE(i)
        end
        [BM,Ib]=min(Bcost(BEd:BEu));
        BE=Ib+BEd-1;
        if BE>N
            BE=N;
        end

        %search for Bin edge 2
        BinE2=zeros(N,1);
        Bcost2=zeros(N,1);
        if BE2<BE
            BE2=BE;
        end
        if BE2-X<1
            BE2d=1;
        else
            BE2d=BE2-X;
        end
        if BE2+X>=N
            BE2u=N;
        else
            BE2u=BE2+X;
        end            
        if BE2d<BE+1
            BE2d=BE+1;
        end
        if BE2u<BE2d
            BE2u=BE2d;
        end
        for i=BE2d:BE2u
            BinE2(i)=i;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
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
                    if t>1
                        if vw(t-1)<BE
                            R(t)=Rgb;
                        elseif (BE<=vw(t-1))&&(vw(t-1)<BinE2(i))
                            R(t)=Rtb;
                        elseif BinE2(i)<=vw(t-1)
                            R(t)=Rbb;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(Rbb-Rgb); 
                            swic(t)=1;
                        end                    
                    end                
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
                        if q(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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
            Bcost2(i)=mean(NC2);
            %BinE2(i)
        end
        [BM2,Ib2]=min(Bcost2(BE2d:BE2u));
        BE2=Ib2+BE2d-1;
        if BE2<BE
            BE2=BE;
        end
        if BE2>N
            BE2=N;
        end

        %Search for Rgb
        RRgb=zeros(NN,1);
        Rgbcost=zeros(NN,1);
        if Rgb-X<1
            Rgbd=1;
        else
            Rgbd=Rgb-X;
        end
        if Rgb+X>=NN
            Rgbu=NN;
        else
            Rgbu=Rgb+X;
        end
        if Rtb<=Rgbu
            Rgbu=Rtb;
        end
        if Rgbu<Rgbd
            Rgbu=Rgbd;
        end
        for i=Rgbd:Rgbu
            RRgb(i)=i-1;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
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
                    if t>1
                        if vw(t-1)<BE
                            R(t)=RRgb(i);
                        elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                            R(t)=Rtb;
                        elseif BE2<=vw(t-1)
                            R(t)=Rbb;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(Rbb-RRgb(i)); 
                            swic(t)=1;
                        end                    
                    end                
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
                        if q(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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
            Rgbcost(i)=mean(NC2);
            %RRgb(i)
        end
        [RgbM, Igb]=min(Rgbcost(Rgbd:Rgbu));
        Rgb=RR(Igb-1+Rgbd);

        %Search for Rtb
        RRtb=zeros(NN,1);
        Rtbcost=zeros(NN,1);
        if Rtb<Rgb
            Rtb=Rgb;
        end
        if Rtb-X<1
            Rtbd=1;
        else
            Rtbd=Rtb-X;
        end
        if Rtb+X>=NN
            Rtbu=NN;
        else
            Rtbu=Rtb+X;
        end
        if Rtbd<Rgb
            Rtbd=Rgb;
        end
        if Rtbu>Rbb
            Rtbu=Rbb;
        end
        if Rtbu<Rtbd
            Rtbu=Rtbd;
        end
        for i=Rtbd:Rtbu
            RRtb(i)=i-1;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
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
                    if t>1
                        if vw(t-1)<BE
                            R(t)=Rgb;
                        elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                            R(t)=RRtb(i);
                        elseif BE2<=vw(t-1)
                            R(t)=Rbb;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(Rbb-RRgb(i)); 
                            swic(t)=1;
                        end                    
                    end                
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
                        if q(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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
            Rtbcost(i)=mean(NC2);
            %RRtb(i)
        end
        [RtbM, Itb]=min(Rtbcost(Rtbd:Rtbu));
        Rtb=RR(Itb-1+Rtbd);
        if Rtb<Rgb
            Rtb=Rgb;
        end

        %search for Rbb
        RRbb=zeros(NN,1);
        Rbbcost=zeros(NN,1);
        if Rbb<Rtb
            Rbb=Rtb;
        end
        if Rbb-X<1
            Rbbd=1;
        else
            Rbbd=Rbb-X;
        end
        if Rbb+X>=NN
            Rbbu=NN;
        else
            Rbbu=Rbb+X;
        end                
        if Rbbd<Rtb
            Rbbd=Rtb;
        end
        if Rbbu<Rbbd
            Rbbu=Rbbd;
        end
        for i=Rbbd:Rbbu
            RRbb(i)=i-1;
            c=1;
            b=0;
            cc=0;
            TESTCONV=zeros(MM,1);
            TESTCONV(1)=1;
            count1=0;
            count2=0;
            while ((cc<1) && (b<MM))
                g=normrnd(mu,sig,Block,N);
                swi=zeros(Block,1);
                swic=zeros(Block,1);
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
                    if t>1
                        if vw(t-1)<BE
                            R(t)=Rgb;
                        elseif (BE<=vw(t-1))&&(vw(t-1)<BE2)
                            R(t)=Rtb;
                            count1=count1+1;
                        elseif BE2<=vw(t-1)
                            R(t)=RRbb(i);
                            count2=count2+1;
                        end
                        if R(t)~=R(t-1)
                            swi(t)=abs(Rbb-RRgb(i)); 
                            swic(t)=1;
                        end                    
                    end                
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
                        if q(t)*F<=g(t,n) %g(t+ind(r),n)
                            vw(t)=vw(t)+1;
                            gbar(t,n)=g(t,n); %g(t,n);
                        end
                    end
                    gbar3(t)=sum(gbar(t,:));
                    A(t)=gam*min(1,R(t)/vw(t));
                    %A(t)=gam*(1-1/(eps^(R(t)/v(t))));
                    aw(t)=binornd(vw(t),A(t));
                    NCw(t)=rho*R(t)+(lam-1)*gbar3(t)+SWrho*swi(t);
                    NCw2(t)=rho*R(t)+(lam-1)*gbar3(t);
                end
                if b<2
                    v=vw;
                    a=aw;
                    gb=gbar3;
                    NC=NCw;
                    NC2=NCw2;
                else
                    vhold=cat(2,v',vw');
                    ahold=cat(2,a',aw');
                    gbhold=cat(2,gb',gbar3');
                    NChold=cat(2,NC',NCw');
                    NC2hold=cat(2,NC2',NCw2');
                    v=vhold';
                    a=ahold';
                    gb=gbhold';
                    NC=NChold';
                    NC2=NC2hold';
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
            Rbbcost(i)=mean(NC2);
            %RRbb(i)
        end
        [RbbM, Ibb]=min(Rbbcost(Rbbd:Rbbu));
        Rbb=RR(Ibb-1+Rbbd);
        if Rbb>N
            Rbb=N;
        end
        if Rbb<Rtb
            Rbb=Rtb;
        end
        crit=abs(costold-RbbM)/RbbM;
        count3=count3+1;
    end
    count(m)=count3;
    BnE1(m)=BE;
    BnE2(m)=BE2;
    RRRgb(m)=Rgb;
    RRRtb(m)=Rtb;
    RRRbb(m)=Rbb;
    CM(m)=RbbM;
    m
end
opt=[BE BE2 Rgb Rtb Rbb];

OPTBB=[BnE1 BnE2 RRRgb RRRtb RRRbb CM];

save('C:\Users\kcarlaw\Documents\MATLAB\ASB2022final\OPTBB.txt','OPTBB','-ascii')