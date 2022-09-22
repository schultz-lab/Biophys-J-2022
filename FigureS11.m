tic;
N=80;
D=50;
Ki=1;

Rr=[1:N]*30/N;
Ha=[1:N]*25/N;
Hr=[1:N]*30/N;

[Ri, Ai, Am, Dm, Re, Ae, De, Reg4]=phenotypes(D , Rr, Ha, Hr);

try
    load P3
catch

% find pareto fronts
P=zeros(size(squeeze(Ri)));

for i=1:length(Rr)
  for j=1:length(Ha)
    for k=i:length(Hr)

        Qam = Am(i,j,k)-Am;
        Qdm = Dm(i,j,k)-Dm;
        Qae = Ae(i,j,k)-Ae;
        Qde = De(i,j,k)-De;
        Q= (Qam>=0) & (Qdm>=0) & (Qae>=0) & (Qde>=0);
        Q(i,j,k)=0;

        if find(Q==1)
            P(i,j,k)=0;
        else
            P(i,j,k)=1;
        end    
        
    end
  end
end

save P3 P

end
toc;

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

CR=1-0.518*rmDm-rmDe;
CG=1-0.765*rmAm-0.267*rmAe-rmDm;
CB=1-rmAm-rmAe-0.031*rmDe;

CR(P==0)=1;
CB(P==0)=1;
CG(P==0)=1;

CC=zeros([size(P) 3]);
CC(:,:,:,1)=CR;
CC(:,:,:,2)=CG;
CC(:,:,:,3)=CB;

% plots
figure(2)
hold off
[HA, RR, HR] = meshgrid(Ha,Rr,Hr);
scatter3(RR(P==1),HR(P==1),HA(P==1),10,[CR(P==1) CG(P==1) CB(P==1)],'filled')
axis([0 max(Rr) 0 max(Hr) 0 max(Ha)]);
ylabel('h_r','FontSize',20)
xlabel('R_r','FontSize',20)
zlabel('h_z','FontSize',20)
hold on;
plot3([1 1],[30 30],[0 25],'k--','LineWidth',2);
plot3([1 1],[0 0],[0 25],'k--','LineWidth',2);
plot3([30 30],[30 30],[0 25],'k--','LineWidth',2);
plot3([0 0],[0 0],[0 25],'k--','LineWidth',2);
plot3([0 30],[0 30],[0 0],'k--','LineWidth',2);
plot3([0 30],[0 30],[25 25],'k--','LineWidth',2);
plot3([1 1],[0 30],[0 0],'k--','LineWidth',2);
plot3([1 1],[0 30],[25 25],'k--','LineWidth',2);
plot3([1 1],[1 1],[0 25],'k--','LineWidth',2);