[Ri, Ai, Am, Dm, Re, Ae, De, Reg4]=phenotypes(D , Rr, Ha, Hr);

try
    load P2
catch

% find pareto fronts
P=zeros(size(squeeze(Ri)));

for i=1:length(Rr)
    for k=1:length(Hr)
        
        Qam = Am(i,1,k)-Am;
        Qdm = Dm(i,1,k)-Dm;
        Qae = Ae(i,1,k)-Ae;
        Qde = De(i,1,k)-De;
        Q= (Qam>=0) & (Qdm>=0) & (Qae>=0) & (Qde>=0);
        Q(i,1,k)=0;

        if find(Q==1)
            P(i,k)=0;
        else
            P(i,k)=1;
        end    
        
    end
end

save P2 P

end

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
d1 = 1 + Ki*D/(Ha+Ki+1);
d2 = 1 + Ki*D/(Ki+1);
plot([0 r],[1 1],'k')
plot([0 r],[0 r],'k')
plot([1 1],[0 h],'k')
plot([0 r],[d1 d1],'k')
plot([0 r],[d2 d2],'k')
plot([1 d2],[d2 d2],'k','LineWidth',3)
plot([d2 d1],[d2 d1],'k','LineWidth',3)
plot([1 d1],[d1 d1],'k','LineWidth',3)
plot([1 1],[d1 d2],'k','LineWidth',3)

ylabel('Repressor expression h_r','FontSize',15)
xlabel('Repressor self-repression threshold R_r','FontSize',15)