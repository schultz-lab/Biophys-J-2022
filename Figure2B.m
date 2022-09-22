D=50;

Rr=[5 1.65 1.1 4 4/1.5 30 20]; 
Ha=[25 25 25 25 25 25 25];
Hr=[0.5 1.65 1.65 4 4 30 30];

for i=1:7
  subplot(2,4,i+floor(i/4))
  p=[Rr(i) Ha(i) Hr(i)];
  objective2(D,p)
end