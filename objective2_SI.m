function objective2_SI(D,p)

Ki=1;
Ka=3000;
K1=0.0001;

plots=0;
lims= [1000,0.03];

regA=@(rf,p) p(2)*heaviside(1-rf);
regR=@(rf,p) p(3)*heaviside(p(1)-rf);
DD=@(t) D*heaviside(t);

df=@(d,r) (d-r-1+sqrt((d-r-1).^2+4*d))/2;
rf=@(d,r) r.*(1./(1+d));

function dx=equations(t,x,p)
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

        
r0=min(p(1),p(3));

if r0>1
    a0=0;
else
    a0=p(2);
end

opts = odeset('Events',@(t,x) events(t,x,p));
[tt,xx]=ode113(@(t,x) equations(t,x,p),[0 10],[r0,a0,0],opts);

    subplot(3,1,1)
    hold off;
    plot(tt,xx(:,2),'r','LineWidth',3)
    hold on
    axis([0 max(tt)/2 0 1.2*p(2)]);
    ylabel('Z','FontSize',16)
    titl=sprintf('K_i=1 D_{out}=%g R_r=%g h_z=%g h_r=%g',D,p(1),p(2),p(3));
    title(titl);
    set(gca,'FontSize',16);
    
    subplot(3,1,2)
    hold off;
    plot(tt,xx(:,1),'g','LineWidth',3)
    hold on
    plot(tt,rf(xx(:,3),xx(:,1)),'--g','LineWidth',2)
    plot([0 tt(end)],[1 1],'-.k','LineWidth',2)
    axis([0 max(tt)/2 0 1.2*p(3)])
    ylabel('R','FontSize',16)
    set(gca,'FontSize',16);
    
    subplot(3,1,3)
    hold off;
    plot(tt,xx(:,3),'b','LineWidth',3);
    hold on
    axis([0 max(tt)/2 0 1.2*max(xx(:,3))])
    set(gca,'FontSize',16);
    ylabel('D','FontSize',16)
    xlabel('time','FontSize',16)

end