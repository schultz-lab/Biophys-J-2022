function tracking(n) 
D=50;
Ki=1;
w=[0 10 1 1 0 1 1];

N1=300;
N2=300;

Rr=[1:N2]/10;
Hr=[1:N1]/10;
Ha=[1:N1]/6;

HA = 25;
d1 = 1 + Ki*D/(HA+Ki+1);
d2 = 1 + Ki*D/(Ki+1);

M=zeros(N1,N2);

lb=[1 d1];
ub=[Inf d2];
A=[1 -1];
b=0;

pareto2D_single;
loops = 201;

for i=1:loops      
    p0=[10 20];
    dd=10^(-1+2*(i-1)/(loops-1));
    %dd=10*(i-1)/(loops-1);
    %w=[0 10 2-dd 1+dd 0 1 1];
    w(n)=dd;
    [pf,fval] = fmincon(@(p) cost(w,D,[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);
    opt(i,:)=pf;
    
end

[pf0,fval] = fmincon(@(p) cost([0 10 1 1 0 1 1],D,[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);

rr=[1 0.5 0];
gg=[0 0.5 1];
gr=[0.7 0.7 0.7];

plot(opt(:,1),opt(:,2),'.','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',7)
plot(opt(1,1),opt(1,2),'ko','MarkerFaceColor',gg,'MarkerSize',7)
plot(opt(end,1),opt(end,2),'ko','MarkerFaceColor',rr,'MarkerSize',7)
plot(pf0(1),pf0(2),'ko','MarkerFaceColor',gr,'MarkerSize',7)
end