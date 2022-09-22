function objective2(D,p)

Ki=1;

regA=@(rf,p) p(2)*heaviside(1-rf);
regR=@(rf,p) p(3)*heaviside(p(1)-rf);
DD=@(t) D*heaviside(t);
rf=@(d,r) r.*(1-d./(1+d));

function dx=equations(t,x,p)
    dx=zeros(3,1);
    dx(1)= regR(rf(x(3),x(1)),p) - x(1);
    dx(2)= regA(rf(x(3),x(1)),p) - x(2);
    dx(3)= Ki*(DD(t)-x(3)) - x(2)*x(3) - x(3);
end

function [value,isterminal,direction]=events(t,x,p)
    dx=equations(t,x,p);
    value=norm(dx./x)-1e-3;
    isterminal=1;
    direction=0;
end

        
r0=min(p(1),p(3));
if r0>1
    a0=0;
else
    a0=p(2);
end

opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) equations(t,x,p),[0 10],[r0,a0,0],opts);

plot([0 3],[1 1],':k','LineWidth',1)
hold on
plot(tt,xx(:,3),'b','LineWidth',3)
plot(tt,xx(:,1),'Color',[0 0.7 0],'LineWidth',3)
plot(tt,rf(xx(:,3),xx(:,1)),'--','Color',[0 0.7 0],'LineWidth',3)
plot(tt,xx(:,2),'r','LineWidth',3)
axis([0 3 0 15])
box off
set(gca,'XTick',[0 1 2 3],'YTick',[0 5 10 15],'FontSize',12)
xlabel('Time','FontSize',15)
ylabel('Concentrations','FontSize',15)

end