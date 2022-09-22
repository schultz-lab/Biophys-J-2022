function [Ri, Ai, Am, Dm, Re, Ae, De, Reg4]=phenotypes(D , Rr, Ha, Hr)

Ki=1;

[HA, RR, HR] = meshgrid(Ha,Rr,Hr);

d1 = Ki*D./(HA+Ki+1);
d2 = Ki*D/(Ki+1);
de = HR - 1;
ae = Ki*D./(HR-1);

Reg0 = RR<1 | HR<1 ;
Reg3 = HR>=d1+1 & HR<=d2+1 & HR<RR ;
Reg4 = RR>=1 & HR>=d1+1 & HR<=d2+1 & HR>=RR ;
Reg5 = HR>d2+1 & HR<RR ;
Reg6 = RR>=1 & HR>d2+1 & HR>=RR ;

t1= (min(RR,HR)-1)/(Ki*D);
F=@(v) (0.61+0.46*v+0.23*v.^2-0.038*v.^3)./(1+0.02*v.^5)+(v+0.05)*0.02.*v.^5./(1+0.02*v.^5);
dm = D*Ki*sqrt(pi./(2*HA)).*F(t1./sqrt(pi./(2*HA)));

X = (RR-1).*HA./(2*(HR-RR)); 
Y = Ki*D*HA./(HR-RR); 
am = -X+sqrt(X.^2+Y);

Ri = min(RR,HR);
Ai = zeros(size(RR));
Ai(Reg0) = HA(Reg0);
Am = HA;
Dm = dm;
Re = HR;
Ae = HA;
De = d1;

%%%%%%%% region 0
Dm(Reg0) = d1(Reg0);

%%%%%%%% region 3
Am(Reg3) = ae(Reg3);
Ae(Reg3) = ae(Reg3);
De(Reg3) = de(Reg3);

%%%%%%%% region 4
Am(Reg4) = am(Reg4);
Ae(Reg4) = ae(Reg4);
De(Reg4) = de(Reg4);

%%%%%%%% region 5
Am(Reg5) = 0;
Dm(Reg5) = d2;
Ae(Reg5) = 0;
De(Reg5) = d2;

%%%%%%%% region 6
Am(Reg6) = 0;
Dm(Reg6) = d2;
Ae(Reg6) = 0;
De(Reg6) = d2;

end