%Plotting PDOS after GetPdosRotated.m function has run, another way to plot
%orbitals neatly.

hold on

%Plot trigonal after rotation. Assumes the rotation varialbe is not set to 0
plot(yorbu(:,1)+yorbu(:,5),xout,'Color','#4819F0',LineWidth=3, DisplayName='Bi Lone Pair')
plot(yorbu(:,2)+yorbu(:,6),xout,'Color','#FF2EE3',LineWidth=3, DisplayName='Bi pz')
plot(yorbu(:,3)+yorbu(:,4)+yorbu(:,7)+yorbu(:,8),xout,'Color','#AD29F1',LineWidth=3, DisplayName='Bi px/py')

plot(yorbu(:,9)+yorbu(:,12),xout,'Color','#77AC30',LineWidth=3,DisplayName='Cl pz')
plot(yorbu(:,10)+yorbu(:,11)+yorbu(:,13)+yorbu(:,14),xout,'Color','#73F300',LineWidth=3,DisplayName='Cl px/py')

plot(yorbu(:,15)+yorbu(:,16)+yorbu(:,18)+yorbu(:,19),xout,'Color','#F31818',LineWidth=3,DisplayName='O px/py')
plot(yorbu(:,17)+yorbu(:,20),xout,'Color','#A2142F',LineWidth=3,DisplayName='O pz')

%Plot PDOS as initially read from Wannier90. Assumes rotation is set to 0.
%plot(yorbu(:,1),xout,'b',yorbu(:,2),xout,'r',yorbu(:,3),xout,'r',yorbu(:,4),xout,'r',yorbu(:,5),xout,'b','LineWidth',10)

%legend('Bi Lone Pair','Bi pz','Bi px / py','Cl pz','Cl px / py','O px / py', 'O pz')
xlabel('Total PDOS')
ylabel('$E_{nk}-E_F$ (eV)','interpreter','latex')
ylim([-12 6])
set(gca,'ytick',[-12:2:6])
%axis([0 max(ytot) -3 1])
set(gca,'fontsize',32)
axis([0 10 -12 6])
grid on
