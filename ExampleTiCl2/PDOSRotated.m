clear all

%if rotation=0, there's no basis rotation. Else the basis is rotated to
%diagonalize the density matrix. Note that this is specific to TiCl2 where
%all the d orbitals of one atom are 
rotation=1

EF= 8.3283;


kgrid=[20 20 20];
hamfilename='ticl2_hr.dat'
orbitals=[1:5];
xout=[-5:0.01:5];
smearing=0.02
soc=0;
%R = vector R at which we want a part of the local Wannier hamiltonian.
%R=[0 0 0] is the local R.
R=[0 0 0];
[ norbu,yorbu,ntotu,ytotu ,rho,rhoe,rhov, Hrask,Hrrot] = GetPdosAndNFromWan(rotation, hamfilename, kgrid, orbitals , EF, smearing, xout,soc,1,R );
% note that after the rotation, small (<10^-19) imaginary occupations may
% be found. 
% norb = occupation/orbital
% yorbu = DOS
% ntotu = total occupation
% ytotu = PDOS
% Hrask, Hrrot = Hamiltonians at R, before and after rotation

hold on

%Plot trigonal after rotation. Assumes the rotation varialbe is not set to 0
plot(yorbu(:,1),xout,'b',yorbu(:,2),xout,'b',yorbu(:,3),xout,'g',yorbu(:,4),xout,'r',yorbu(:,5),xout,'r','LineWidth',10)

%Plot PDOS as initially read from Wannier90. Assumes rotation is set to 0.
%plot(yorbu(:,1),xout,'b',yorbu(:,2),xout,'r',yorbu(:,3),xout,'r',yorbu(:,4),xout,'r',yorbu(:,5),xout,'b','LineWidth',10)

legend('PDOS')
xlabel('Total PDOS')
ylabel('E-E_F')
%axis([0 max(ytot) -3 1])
set(gca,'fontsize',25)

grid on



