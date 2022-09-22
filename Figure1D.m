D=20;
p=[0.5 5 1.5];

Ki=1;
Ka=3000;
K1=0.0001;

plots=0;
lims= [1000,0.03];

regA=@(rf,p) p(2)*heaviside(p(1)-rf);
regR=@(rf,p) p(3)*heaviside(1-rf);
DD=@(t) D*heaviside(t);

df=@(d,r) (d-r-1+sqrt((d-r-1).^2+4*d))/2;
rf=@(d,r) r.*(1-d./(1+d));

r0=min(1,p(3));

if r0>p(1)
    a0=0;
else
    a0=p(2);
end

opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) equations(t,x,p),[0 10],[r0,a0,0],opts);
    
plot(tt,0.5*xx(:,3)/max(xx(:,3)),'b','LineWidth',5)
hold on
plot([-0.5 0],[0 0],'b','LineWidth',5)
plot([-0.5 tt(end)],0.3*p(1)*[1 1]/p(3),':k','LineWidth',2)
plot(tt,0.3*xx(:,1)/p(3),'Color',[0 0.7 0],'LineWidth',5)
plot([-0.5 0],0.3*[r0 r0]/p(3),'Color',[0 0.7 0],'LineWidth',5)
plot(tt,0.3*rf(xx(:,3),xx(:,1))/p(3),':','Color',[0 0.7 0],'LineWidth',5)
plot(tt,0.5*xx(:,2)/max(xx(:,2)),'r','LineWidth',5)
axis([-0.5 max(tt)/3 0 0.5])
set(gca,'Box','off','XTick',[],'YTick',[],'LineWidth',2)

function dx=equations(t,x,p)
    Ki=1;
    D=10;
    rf=@(d,r) r.*(1-d./(1+d));
    regA=@(rf,p) p(2)*heaviside(p(1)-rf);
    regR=@(rf,p) p(3)*heaviside(1-rf);
    DD=@(t) D*heaviside(t);
    K1=0.0001;
    dx=zeros(3,1);
    dx(1)= regR(rf(x(3),x(1)),p) - x(1);
    dx(2)= regA(rf(x(3),x(1)),p) - x(2);
    dx(3)= Ki*(DD(t)-x(3)) - x(2)*x(3) / (K1*x(3) + 1) - x(3);
end

function [value,isterminal,direction]=events(t,x,p)
    dx=equations(t,x,p);
    value=norm(dx./x)-1e-3;
    isterminal=1;
    direction=0;
end