function c=cost(weights, D , p)

Ki=1;

%%% concentrations of interest, par = [Rr Ha Hr]
d1 = Ki*D/(p(2)+Ki+1);
d2 = Ki*D/(Ki+1);
de = abs(p(3) - 1);
ae = abs(Ki*D/(p(3)-1));

t1= (min(p(1),p(3))-1)/(Ki*D);
F=@(v) (0.61+0.46*v+0.23*v.^2-0.038*v.^3)./(1+0.02*v.^5)+(v+0.05)*0.02.*v.^5./(1+0.02*v.^5);
dm = D*Ki*sqrt(pi./(2*p(2))).*F(t1./sqrt(pi./(2*p(2))));

X = (p(1)-1)*p(2)/(2*(p(3)-p(1))); 
Y = Ki*D*p(2)/(p(3)-p(1));
am = -X+sqrt(X^2+Y);

%%%%%%%% region 0
if p(3)<1 || p(1)<1
    c = weights*[ min(p(1),p(3)) p(2) p(2) d1 p(3) p(2) d1 ]';

%%%%%%%% region 1
elseif p(3)<d1+1
    if p(3)<=p(1)
        c = weights*[ p(3) 0 p(2) dm p(3) p(2) d1 ]';

%%%%%%%% region 2
    else
        c = weights*[ p(1) 0 p(2) dm p(3) p(2) d1 ]';
    end

%%%%%%%% region 3
elseif p(3)<d2+1
    if p(3)<=p(1)
        c = weights*[ p(3) 0 ae dm p(3) ae de ]';

%%%%%%%% region 4
    else
        c = weights*[ p(1) 0 am dm p(3) ae de ]';
    end

%%%%%%%% region 5
else
    if p(3)<=p(1)
        c = weights*[ p(3) 0 0 d2 p(3) 0 d2 ]';

%%%%%%%% region 6
    else
        c =  weights*[ p(1) 0 0 d2 p(3) 0 d2 ]';
    end
    
end

end
