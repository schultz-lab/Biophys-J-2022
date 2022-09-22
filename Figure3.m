weights1 = [0 10 2 1 0 1 1];
weights2 = [0 10 1 2 0 1 1];
weights3 = [1 10 1 2 1 1 1];
D=50;
Ki=1;

N1=300;
N2=300;

Rr=[1:N2]/10;
Hr=[1:N1]/10;
Ha=[1:N1]/6;

HA = 25;

M1=zeros(N1,N2);
M2=zeros(N1,N2);
M3=zeros(N1,N2);
M4=zeros(N1,N2);
M5=zeros(N1,N2);
M6=zeros(N1,N2);

for i=1:N1
    for j=1:N2
        p1=[Rr(j) HA Hr(i)];
        M1(i,j)=cost(weights1,D,p1);
        M3(i,j)=cost(weights2,D,p1);
        M5(i,j)=cost(weights3,D,p1);
    end
end

[bla ii1]=min(M1(:,end));
HR1=Hr(ii1);
i2=find(Rr>1,1);
[bla ii2]=min(M3(:,i2));
HR2=Hr(ii2);
[p3,fval] = patternsearch(@(p) cost(weights3,D,p),[1 HA HR2]);
HA3=p3(2);
HR3=p3(3);

d1 = 1 + Ki*D/(HA+Ki+1);
d2 = 1 + Ki*D/(Ki+1);
d3 = 1 + Ki*D/(HA3+Ki+1);

ae1=Ki*D/(HR1-1)-Ki-1;
ae2=Ki*D/(HR2-1)-Ki-1;
ae3=Ki*D/(HR3-1)-Ki-1;

for i=1:N1
    for j=1:N2
        p2=[Rr(j)/2 Ha(i) HR1];
        M2(i,j)=cost(weights1,D,p2);
        p2=[Rr(j)/2 Ha(i) HR2];
        M4(i,j)=cost(weights2,D,p2);
        p2=[Rr(j)/2 Ha(i) HR3];
        M6(i,j)=cost(weights3,D,p2);
    end
end

opt=zeros(51,2);
lb=[1 d1];
ub=[Inf d2];
A=[1 -1];
b=0;

for i=1:51      
    p0=[10 20];
    w=[0 10 1+(i-1)/50 2-(i-1)/50 0 1 1];
    [pf,fval] = fmincon(@(p) cost(w,D,[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);
    opt(i,:)=pf; 
end

%% Figure 3A
subplot(3,3,1)
mmin=min(min(M1));
MM1=M1-mmin+1;
image(Rr,Hr,log10(MM1),'CDataMapping','scaled')
hold on
r=max(Rr);
h=max(Hr);
a=max(Ha);
x=[0 0 r d2 1 1];
y=[d1 h h d2 d2 d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
x=[0 0 r r d1];
y=[d1 0 0 h d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
plot([0 r],[1 1],'k')
plot([0 r],[0 r],'k')
plot([1 1],[0 h],'k')
plot([0 r],[d1 d1],'k')
plot([0 r],[d2 d2],'k')
plot([0 r],[HR1 HR1],'--k')
plot(HR1,HR1,'ko','MarkerFaceColor',[1 0.5 0],'MarkerSize',8)
set(gca,'YDir','normal','XTick',[0 10 20 30],'YTick',[0 10 20 30])
hh=colorbar;
caxis([0 1.5])
set(hh,'Ticks',[0 log10([50 60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'50','60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)

%% Figure 3B
subplot(3,3,2)
mmin=min(min(min(M1)),min(min(M2)));
MM1=M1-mmin+1;
MM2=M2-mmin+1;

for i=1:300
    HHa(:,i)=Ha;
end

surf(Rr/2,HR1*ones(N1,1),HHa,log10(MM2))
hold on
surf(Rr,Hr,HA*ones(N1,N1),log10(MM1))

plot3([0 r/2],[HR1 HR1],[1 1],'k')
plot3([1 1],[HR1 HR1],[0 a],'k')
plot3([HR1 HR1],[HR1 HR1],[0 a],'k')
plot3([0 r/2],[HR1 HR1],[ae1 ae1],'k')
plot3([0 r],[HR1 HR1],[HA HA],'--k')

[X,Y,Z] = sphere;
rad=0.35;
C=zeros(21,21,3);
C(:,:,1)=1;
C(:,:,2)=0.5;
C(:,:,3)=0;
surf(rad*X+HR1,rad*Y+HR1,2*rad*Z+HA,C)

plot3([0 r],[1 1],[HA HA],'k')
plot3([0 r],[0 r],[HA HA],'k')
plot3([1 1],[0 h],[HA HA],'k')
plot3([0 r],[d3 d3],[HA HA],'k')
plot3([0 r],[d2 d2],[HA HA],'k')
plot3([0 r],[HR1 HR1],[HA HA],'--k')

grid off
shading flat
hh=colorbar;
caxis([0 1.5])
set(hh,'Ticks',[0 log10([50 60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'50','60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)
zlabel('h_z','FontSize',15)
axis([0 15 0 15 20 50])
campos([150.1604 -150.9585  270.4557])

%% Figure 3D
subplot(3,3,4)
mmin=min(min(M3));
MM3=M3-mmin+1;
image(Rr,Hr,log10(MM3),'CDataMapping','scaled')
hold on
r=max(Rr);
h=max(Hr);
a=max(Ha);
x=[0 0 r d2 1 1];
y=[d1 h h d2 d2 d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
x=[0 0 r r d1];
y=[d1 0 0 h d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
plot([0 r],[1 1],'k')
plot([0 r],[0 r],'k')
plot([1 1],[0 h],'k')
plot([0 r],[d1 d1],'k')
plot([0 r],[d2 d2],'k')
plot([0 r],[HR2 HR2],'--k')

for i=1:51
        c=[0.85 0.85 0.85];
        plot(opt(i,1),opt(i,2),'.','MarkerFaceColor',c,'MarkerEdgeColor',c)
end

plot(HR1,HR1,'ko','MarkerFaceColor',[1 0.5 0],'MarkerSize',8)
plot(1,HR2,'ko','MarkerFaceColor',[0.85 0.5 1],'MarkerSize',8)
set(gca,'YDir','normal','XTick',[0 10 20 30],'YTick',[0 10 20 30])
hh=colorbar;
caxis([0 1.5])
set(hh,'Ticks',[0 log10([50 60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'50','60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)

%% Figure 3E
subplot(3,3,5)
mmin=min(min(min(M3)),min(min(M4)));
MM3=M3-mmin+1;
MM4=M4-mmin+1;
for i=1:300
HHa(:,i)=Ha;
end
surf(Rr/2,HR2*ones(N1,1),HHa,log10(MM4))
hold on
surf(Rr,Hr,HA*ones(N1,N1),log10(MM3))

plot3([0 r/2],[HR2 HR2],[1 1],'k')
plot3([1 1],[HR2 HR2],[0 a],'k')
plot3([HR2 HR2],[HR2 HR2],[0 a],'k')
plot3([0 r/2],[HR2 HR2],[ae2 ae2],'k')
plot3([0 r],[HR2 HR2],[HA HA],'--k')

[X,Y,Z] = sphere;
rad=0.35;
C=zeros(21,21,3);
C(:,:,1)=0.85;
C(:,:,2)=0.5;
C(:,:,3)=1;
surf(rad*X+1,rad*Y+HR2,2*rad*Z+HA,C)

plot3([0 r],[1 1],[HA HA],'k')
plot3([0 r],[0 r],[HA HA],'k')
plot3([1 1],[0 h],[HA HA],'k')
plot3([0 r],[d3 d3],[HA HA],'k')
plot3([0 r],[d2 d2],[HA HA],'k')
plot3([0 r],[HR2 HR2],[HA HA],'--k')

grid off
shading flat
hh=colorbar;
caxis([0 1.5])
set(hh,'Ticks',[0 log10([50 60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'50','60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)
zlabel('h_z','FontSize',15)
axis([0 15 0 15 20 50])
campos([150.1604 -150.9585  270.4557])

%% Figure 3G
subplot(3,3,7)
mmin=min(min(M5));
MM5=M5-mmin+1;
image(Rr,Hr,log10(MM5),'CDataMapping','scaled')
hold on
r=max(Rr);
h=max(Hr);
a=max(Ha);
x=[0 0 r d2 1 1];
y=[d1 h h d2 d2 d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
x=[0 0 r r d1];
y=[d1 0 0 h d1];
patch(x,y,[1 1 1],'EdgeColor','none','FaceAlpha',0.5)
plot([0 r],[1 1],'k')
plot([0 r],[0 r],'k')
plot([1 1],[0 h],'k')
plot([0 r],[d3 d3],'k')
plot([0 r],[d2 d2],'k')
plot([0 r],[HR3 HR3],'--k')
plot(1,HR3,'ko','MarkerFaceColor',[0.85 0.5 1],'MarkerSize',8)
set(gca,'YDir','normal','XTick',[0 10 20 30],'YTick',[0 10 20 30])
hh=colorbar;
caxis([0 1.5])
set(hh,'Ticks',[0 log10([60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)

%% Figure 3H
subplot(3,3,8)
mmin=min(min(min(M5)),min(min(M6)));
MM5=M5-mmin+1;
MM6=M6-mmin+1;

for i=1:300
HHa(:,i)=Ha;
end
surf(Rr/2,HR3*ones(N1,1),HHa,log10(MM6))
hold on
surf(Rr,Hr,HA3*ones(N1,N1),log10(MM5))

plot3([0 r/2],[HR3 HR3],[1 1],'k')
plot3([1 1],[HR3 HR3],[0 a],'k')
plot3([HR3 HR3],[HR3 HR3],[0 a],'k')
plot3([0 r/2],[HR3 HR3],[ae3 ae3],'k')
plot3([0 r],[HR3 HR3],[HA3 HA3],'--k')

[X,Y,Z] = sphere;
rad=0.35;
C=zeros(21,21,3);
C(:,:,1)=0.85;
C(:,:,2)=0.5;
C(:,:,3)=1;
surf(rad*X+1,rad*Y+HR3,2*rad*Z+HA3,C)

plot3([0 r],[1 1],[HA3 HA3],'k')
plot3([0 r],[0 r],[HA3 HA3],'k')
plot3([1 1],[0 h],[HA3 HA3],'k')
plot3([0 r],[d3 d3],[HA3 HA3],'k')
plot3([0 r],[d2 d2],[HA3 HA3],'k')
plot3([0 r],[HR3 HR3],[HA3 HA3],'--k')

grid off
shading flat
hh=colorbar;
caxis([0 1.2])
set(hh,'Ticks',[0 log10([60 70 80]-mmin+1)],'TickLabels',{num2str(round(mmin)),'60','70','80'},'FontSize',12)
ylabel(hh,'Cost','FontSize',15)
set(gca,'FontSize',12)
xlabel('R_r','FontSize',15)
ylabel('h_r','FontSize',15)
zlabel('h_z','FontSize',15)
axis([0 15 0 15 20 50])
campos([150.1604 -150.9585  270.4557])

p1=[11.3 25 11.3];
p2=[1.1 25 10.3];
p3=[1.1 28.7869 8.0889];

%% Figure 3C
subplot(3,3,3)
objective2(50,p1)

%% Figure 3F
subplot(3,3,6)
objective2(50,p2)

%% Figure 3I
subplot(3,3,9)
objective2(50,p3)