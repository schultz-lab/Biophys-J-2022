Ki=1;
w=[0 10 1 10 0 1 1];

N1=300;
N2=300;

Rr=[1:N2]/10;
Hr=[1:N1]/10;
Ha=[1:N1]/6;

HA = 25;
M=zeros(N1,N2);
loops = 201;

for i=1:loops      
    p0=[10 20];
    D(i)=0.5*10^(1+2*(i-1)/(loops-1));
    d1 = 1 + Ki*D(i)/(HA+Ki+1);
    d2 = 1 + Ki*D(i)/(Ki+1);
    lb=[1 d1];
    ub=[Inf d2];
    A=[1 -1];
    b=0;
    [pf,fval] = fmincon(@(p) cost(w,D(i),[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);
    opt(i,:)=pf;
    
end

rr=[1 0.5 0];
gg=[0 0.5 1];
p1 = plot(D,opt(:,1),'Color',rr,'LineWidth',3);
hold on
p2 = plot(D,opt(:,2),'Color',gg,'LineWidth',3);
set(gca,'XScale','log','FontSize',12)
ylabel('Optimal parameters','FontSize',15)
xlabel('Extracellular drug D_{out}','FontSize',15)
axis([5 500 0 30])
p3 = plot([D(61) D(61)],[0 opt(61,2)],'--k','LineWidth',3); %vertical
p4 = plot([D(101) D(101)],[0 opt(101,2)],'--r','LineWidth',3);
p5 = plot([D(161) D(163)],[0 opt(163,2)],'--k','LineWidth',3);
p6 = plot([D(1) D(61)],[opt(61,2) opt(61,2)],'--k','LineWidth',3);
p7 = plot([D(1) D(101)],[opt(101,2) opt(101,2)],'--r','LineWidth',3);
p8 = plot([D(1) D(163)],[opt(163,2) opt(163,2)],'--k','LineWidth',3);
legend([p1,p2],'R_r','h_r','Location','northwest')