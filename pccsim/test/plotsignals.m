function plotsignals(t, sp,si,sq )
%PLOTSIGNALS Summary of this function goes here
%   Detailed explanation goes here

figure;
ax1 = subplot(3,1,1);plot(t,abs(sp));
xlabel('Time [\mu s]')
ylabel('Level');
grid;
axis([0 60e-6 -0.5 0.5]);
set(gca,'XTick',0:5e-6:60e-6);
set(gca,'YTick',-0.5:0.1:0.5);
set(gca,'XTickLabel',0:5:60);
%title(['Line ',num2str(line),'- ENV']);
title(['ENV']);
ax2 = subplot(3,1,2);plot(t,si);
xlabel('Time [\mu s]')
ylabel('Level');
grid;
axis([0 60e-6 -0.5 0.5]);
set(gca,'XTick',0:5e-6:60e-6);
set(gca,'YTick',-0.5:0.1:0.5);
set(gca,'XTickLabel',0:5:60);
title(['I']);
ax3 = subplot(3,1,3);plot(t,sq);
xlabel('Time [\mu s]')
ylabel('Level');
grid;
axis([0 60e-6 -0.5 0.5]);
set(gca,'XTick',0:5e-6:60e-6);
set(gca,'YTick',-0.5:0.1:0.5);
set(gca,'XTickLabel',0:5:60);
title(['Q']);
linkaxes([ax1,ax2,ax3],'xy');

end

