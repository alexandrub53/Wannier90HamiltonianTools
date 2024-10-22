%To graph rotated PDOS and to obtain Hamiltonian, rotated orbitals 

clear all

%if rotation=0, there's no basis rotation. Else the basis is rotated to
%diagonalize the density matrix. Note that this is specific to TiCl2 where
%all the d orbitals of one atom are 
rotation=1;

EF= 8.0917;                         %Fermi Energy
kgrid=[20 20 20];                   %k point grid
hamfilename='BiOCl_hr.dat'          %DFT Bandstructure

orbitals=[1:20];                    %Orbitals to be considered from Wannier90
xout=[-15:0.01:10];                 %Energy limits
smearing=0.1                        %Peak smearing/broadening level
soc=0;                              %spin orbit coupling
on_bitals = 1-20;                   %orbitals to be left "on" or "off"
off_bitals =0;


%R = vector R at which we want a part of the local Wannier hamiltonian.
%R=[0 0 0] is the local R.
R=[0 0 0];
[ norbu,yorbu,eorb,ntotu,ytotu ,rho,rhoe,rhov, enkr,enkd,pnkr,pnkd,Hkm,Hkdm] = GetPdosAndNFromWan(rotation, hamfilename, kgrid, orbitals , EF, smearing, xout,soc,1,R,on_bitals,off_bitals );
% note that after the rotation, small (<10^-19) imaginary occupations may
% be found. 
% norb = occupation/orbital
% yorbu = DOS
% ntotu = total occupation
% ytotu = PDOS
% Hrask, Hrzeros = Hamiltonians at R, before and after zero'd entries

hold on

%Plot each orbital after rotation. Assumes the rotation varialbe is not set to 0
plot(yorbu(:,1),xout,'b','LineWidth',5, 'DisplayName', 'Bi1')
plot(yorbu(:,5),xout,'Color','#0072BD',LineWidth=2.5,DisplayName='Bi2')
plot(yorbu(:,9),xout,'g',LineWidth=5,DisplayName='Cl1')
plot(yorbu(:,12),xout,'Color','#77AC30',LineWidth=2.5,DisplayName='Cl2')
plot(yorbu(:,15),xout,'r',LineWidth=5,DisplayName='O1')
plot(yorbu(:,18),xout,'Color','#A2142F',LineWidth=2.5,DisplayName='O2')

plot(yorbu(:,2),xout,'b',LineWidth=5)
plot(yorbu(:,3),xout,'b',LineWidth=5)
plot(yorbu(:,4),xout,'b',LineWidth=5)

plot(yorbu(:,6),xout,'Color','#0072BD',LineWidth=2.5)
plot(yorbu(:,7),xout,'Color','#0072BD',LineWidth=2.5)
plot(yorbu(:,8),xout,'Color','#0072BD',LineWidth=2.5)

plot(yorbu(:,10),xout,'g',LineWidth=5)
plot(yorbu(:,11),xout,'g',LineWidth=5)

plot(yorbu(:,13),xout,'Color','#77AC30',LineWidth=2.5)
plot(yorbu(:,14),xout,'Color','#77AC30',LineWidth=2.5)

plot(yorbu(:,16),xout,'r',LineWidth=5)
plot(yorbu(:,17),xout,'r',LineWidth=5)
plot(yorbu(:,19),xout,'Color','#A2142F',LineWidth=2.5)
plot(yorbu(:,20),xout,'Color','#A2142F',LineWidth=2.5)


%Plot PDOS as initially read from Wannier90. Assumes rotation is set to 0.
%plot(yorbu(:,1),xout,'b',yorbu(:,2),xout,'r',yorbu(:,3),xout,'r',yorbu(:,4),xout,'r',yorbu(:,5),xout,'b','LineWidth',10)

legend('Bi1','Bi2','Cl','Cl2','O1','O2')
xlabel('Total PDOS')
ylabel('E-E_F')
ylim([-12 10])
%axis([0 max(ytot) -3 1])
set(gca,'fontsize',25)

grid on


%To find the band gap
[minenkd,maxenkd] = bounds(enkd,1);
for i = 1:length(maxenkd)
    if maxenkd(i) < EF
        Ms = maxenkd(i);
    else
        Ms;
    end
end
for j = 1:length(minenkd)
    if minenkd(j) > EF
        Mb = minenkd(j);
        break
    else 
        continue
    end
end
Egap = Mb - Ms
