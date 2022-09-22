D=50;

Rr=[5 5 1.2 5 3 34 27];
Ha=[25 25 25 25 25 25 25];
Hr=[0.5 1.4 1.4 2.5 7 30 34];

Hrmax=35;
Rrmax=35;
PRINT=1;

%% Figures S2-S8
for i=0:6
figure(i+1);
p=[Rr(i+1) Ha(i+1) Hr(i+1)];
objective2_SI(D,p)
fil=sprintf('example%d.pdf',i);
if PRINT
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[0.9*pos(3), pos(4)])
    print(gcf,'-dpdf',fil, '-opengl')
end
end

%% Figure S1
figure(8)
hold off;
Ha=25;
Ki=1;
a1=Ki*D/(Ki+Ha+1);
a2=Ki*D/(Ki+1);
plot([1 Rrmax],[1 Hrmax],'k','LineWidth',2);
hold on;
plot([1 Rrmax],[1 1],'k','LineWidth',2);
plot([1 Rrmax],[a1 a1],'k','LineWidth',2);
plot([1 Rrmax],[a2 a2],'k','LineWidth',2);
plot([1 1],[1 Hrmax],'k','LineWidth',2);
xlabel('R_r','FontSize',16);
ylabel('h_r','FontSize',16);
scatter(Rr,Hr,80,'filled', 'r');
text(Rr-1.,Hr,{'0','1','2','3','4','5','6'},'FontSize',16, 'Color', 'r');
text(Rrmax+0.3,a1,'\alpha_1','FontSize',20);
text(Rrmax+0.3,a2,'\alpha_2','FontSize',20);
set(gca,'FontSize',16);
if PRINT
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[0.9*pos(3), pos(4)])
    print(gcf,'-dpdf','HrRr.pdf', '-opengl')
end
