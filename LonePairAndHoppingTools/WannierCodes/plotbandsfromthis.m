clear all

%Simple plotting scheme to compare the Wannier hamiltonian (hr_dat) with
%the output of a band structure from quantum espresso, as extracted from the output file using
%getbands.py. Can easily be modified for additional functionality (PDOS,
%Fermi Surfaces etc).

emin=-12
emax=6
EF=  8.0917;

%load the data file generated from the QE output.
load bands.dat
bands = bands(:,2:end);

grid

set(gca,'XTick',[1 21 41 61 81 101])
set(gca,'XTickLabel',{'\Gamma'; 'K'; 'M'; 'L'; 'A' ;'\Gamma'})
ylabel('$E_{nk}-E_F$ (eV)','interpreter','latex')
set(gca,'fontsize',28)

%plot the bands from QE
plot(bands-EF,'k-','linewidth',1)
set(gca,'ytick',[emin:2:emax])
hold on

axis([0 101 emin emax])


%load some k-points we use for the Wannier tight-binding hamiltonian
pointsforbandsh

%normalize them correctly.
kaplot = kaplot*2*pi;
klist=kaplot

on_bitals = 1:20;
off_bitals = 0;

%Rotate the given data for each atom
[rot,rote,rotv] = RotationMatrixThatIWantToUse('BiOCl_hr.dat',0,0,0);
[V,D]=eig(rot(1:4,1:4));
[V2,D2]=eig(rot(5:8,5:8));
[V3,D3]=eig(rot(9:11,9:11));
[V4,D4]=eig(rot(12:14,12:14));
[V5,D5]=eig(rot(15:17,15:17));
[V6,D6]=eig(rot(18:20,18:20));
rotv=zeros(20);
rotv(1:4,1:4)=rotv(1:4,1:4)+V;
rotv(5:8,5:8)=rotv(5:8,5:8)+V2;
rotv(9:11,9:11)=rotv(9:11,9:11)+V3;
rotv(12:14,12:14)=rotv(12:14,12:14)+V4;
rotv(15:17,15:17)=rotv(15:17,15:17)+V5;
rotv(18:20,18:20)=rotv(18:20,18:20)+V6;


%load the Wannier hamiltonian. Generate the E_nk and the projectors onto
%each orbital (p_nk). See inside the function description for more detail.
[enk,pnk,enkd,pnkd,Hkm,Hkdm] = wanbandsrotorbs('BiOCl_hr.dat',kaplot,rotv,on_bitals,off_bitals);
nk = size(kaplot,1);
set(gca,'XTick',[1 21 41 61 81 101])
set(gca,'XTickLabel',{'\Gamma'; 'K'; 'M'; 'L'; 'A'; '\Gamma'})
ylabel('$E_{nk}-E_F$ (eV)','interpreter','latex')
set(gca,'fontsize',28)
grid on

%Plot the Wannier bands from the Hamiltonian. They should match the bands
%from QE.

r = 32 / 255;
g = 85 / 255;
b = 107 / 255;

wan = plot(real(enkd)-EF, 'o', 'linewidth',1);


%Just for fun, some scatter plots for two orbitals. In this case these are
%the so-called 'projected bands'. In this example we have 20 orbitals, so the
%bands have 20 possible energies. We also have 20 possible orbitals.


hold on
orb1=1;
orb2=2;
orb3=3;
orb4=4;

for eloc=[1:20]
s1 = scatter([1:101]', enkd(1:101,eloc)-EF, (pnkd(1:101,orb4,eloc))*400+0.02, DisplayName='Bi pz');
s2 = scatter([1:101]', enkd(1:101,eloc)-EF, (pnkd(1:101,orb2,eloc)+pnkd(1:101,orb3,eloc))*400+0.02, DisplayName='Bi px / py');
s3 = scatter([1:101]', enkd(1:101,eloc)-EF, (pnkd(1:101,orb1,eloc))*400+0.03, DisplayName = 'Bi Lone Pair');

%Color in rgb
s1.MarkerEdgeColor =[(115/255) (243/255) (0/255)];
s2.MarkerEdgeColor =[(162/255) (20/255) (47/255)];
s3.MarkerEdgeColor =[(243/255) (24/255) (24/255)];
s4.MarkerEdgeColor =[(119/255) (172/255) (48/255)];
legend([s1,s2,s3])

hold on
end

axis([1 101 emin emax])


