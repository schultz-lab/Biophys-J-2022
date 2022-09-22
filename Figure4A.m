N=300;
D=50;
Ki=1;

Rr=[1:N]*30/N;
Ha=25;
Hr=[1:N]*30/N;

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

ylabel('h_r','FontSize',20)
xlabel('R_r','FontSize',20)

set(gca,'Color','none')
set(gca,'Box','on')
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])