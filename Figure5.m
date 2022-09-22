N=300;
D=150;
Ki=1;

Rr=[1:N]*60/N;
Ha=25;
Hr=[1:N]*60/N;

[Ri, Ai, Am, Dm, Re, Ae, De, Reg4]=phenotypes(D , Rr, Ha, Hr);

load P2

i1=find(Rr>=1,1);
Reg=false(N,1,N);

for i=i1:N
    Reg(i1:i,1,i)=1;
end

d1 = 1 + Ki*D/(Ha+Ki+1);
d2 = 1 + Ki*D/(Ki+1);
    
% find optimal values for each phenotype
Rimin = min(Ri(Reg4==1));
Aimin = min(Ai(Reg4==1));
Ammin = min(Am(Reg4==1));
Dmmin = min(Dm(Reg4==1));
Remin = min(Re(Reg4==1));
Aemin = min(Ae(Reg4==1));
Demin = min(De(Reg4==1));

% colors: blue optimizes Dm and De, red optimizes Am and Ae 
f=@(x,xmin) min(1,1./(1+(x-xmin).^3));

rmDm=f(squeeze(Dm),Dmmin);
rmAm=f(squeeze(Am),Ammin);
rmDe=f(squeeze(De),Demin);
rmAe=f(squeeze(Ae),Aemin);

CR=1-0.25*(0.518*rmDm+rmDe);
CG=1-0.25*(0.765*rmAm+0.267*rmAe+rmDm);
CB=1-0.25*(rmAm+rmAe+0.031*rmDe);

CR(Reg==0)=1;
CB(Reg==0)=1;
CG(Reg==0)=1;

P=P';
CC=zeros([size(P) 3]);
CC(:,:,1)=CR';
CC(:,:,2)=CG';
CC(:,:,3)=CB';

% plots
image(Rr,Hr,CC)
hold on
set(gca,'YDir','normal')
r=max(Rr);
h=max(Hr);
CL=0.85;
plot([1 60],[1 60],'Color',CL*[1 1 1],'LineWidth',3)
plot([1 1],[1 60],'Color',CL*[1 1 1],'LineWidth',3)
axis([0 60 0 60])
ylabel('Repressor expression h_r','FontSize',15)
xlabel('Repressor self-repression threshold R_r','FontSize',15)

rr=[0.85 0.5 1];
gg=[1 0.85 0];
plot(54,57,'ko','MarkerSize',10,'MarkerFaceColor',gg,'LineWidth',1.5) % lacI [1 0.75 0]
plot(2,10,'ko','MarkerSize',10,'MarkerFaceColor',rr,'LineWidth',1.5) % tetr [0.5 0.8 1]
plot(1,20,'ko','MarkerSize',10,'MarkerFaceColor',rr,'LineWidth',1.5) % marRAB
plot(1,8,'ko','MarkerSize',10,'MarkerFaceColor',gg,'LineWidth',1.5) % iclR
plot(1,14.5,'ko','MarkerSize',10,'MarkerFaceColor',gg,'LineWidth',1.5) % mtlR
plot(4.5,13.4,'ko','MarkerSize',10,'MarkerFaceColor',gg,'LineWidth',1.5) % galS
plot(44,48,'ko','MarkerSize',10,'MarkerFaceColor',gg,'LineWidth',1.5) % galR