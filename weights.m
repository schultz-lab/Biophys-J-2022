function weights(n)
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

loops = 201;

for i=1:loops      
    p0=[10 20];
    dd(i)=10^(-1+2*(i-1)/(loops-1));
    w(n)=dd(i);
    [pf,fval] = fmincon(@(p) cost(w,D,[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);
    opt(i,:)=pf;
    
end

[pf0,fval] = fmincon(@(p) cost([0 10 1 1 0 1 1],D,[p(1) HA p(2)]),p0,A,b,[],[],lb,ub);

rr=[1 0.5 0];
gg=[0 0.5 1];

if n==3
    opt(167,:) = [];
    dd(167) = [];
end

plot(dd,opt(:,1),'Color',rr,'LineWidth',3)
hold on
plot(dd,opt(:,2),'Color',gg,'LineWidth',3)
legend('R_r','h_r')
set(gca,'XScale','log')
ylabel('Optimal parameters','FontSize',15)
xlabel('Weight','FontSize',15)
axis([0.1 10 0 30])
end