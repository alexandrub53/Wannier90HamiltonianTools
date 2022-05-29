clear all


%Simple plotting scheme to compare the Wannier hamiltonian (hr_dat) with
%the output of a band structure from quantum espresso, as extracted from the output file using
%getbands.py. Can easily be modified for additional functionality (PDOS,
%Fermi Surfaces etc).

emin=-7
emax=3.0
EF=  8.3283;

%load the data file generated from the QE output.
load bands.dat
bands = bands(:,2:end);

grid

set(gca,'XTick',[1 21 41 61 81 101 ])
set(gca,'XTickLabel',{'\Gamma'; 'K'; 'M' ;'L'; 'A'; '\Gamma'})
ylabel('$E_{nk}-E_F$ (eV)','interpreter','latex')
set(gca,'fontsize',18)

%plot the bands from QE
plot(bands-EF,'k-','linewidth',1)
set(gca,'ytick',[emin:1.0:emax])
hold on

axis([0 101 emin emax])

%load some k-points we use for the Wannier tight-binding hamiltonian
pointsforbandsh

%normalize them correctly.
kaplot = kaplot*2*pi;
klist=kaplot


%load the Wannier hamiltonian. Generate the E_nk and the projectors onto
%each orbital (p_nk). See inside the function description for more detail.
[enk,pnk] = wanbands('ticl2_hr.dat',kaplot);
nk = size(kaplot,1);
set(gca,'XTick',[1 21 41 61 81 101 ])
set(gca,'XTickLabel',{'\Gamma'; 'K'; 'M' ;'L'; 'A'; '\Gamma'})
ylabel('$E_{nk}-E_F$ (eV)','interpreter','latex')
set(gca,'fontsize',60)
grid on

%Plot the Wannier bands from the Hamiltonian. They should match the bands
%from QE.
plot(enk-EF,'ok','linewidth',1)

%Just for fun, some scatter plots for two orbitals. In this case these are
%the so-called 'Fat bands'. In this example we have 5 orbitals, so the
%bands have 5 possible energies. We also have 5 possible orbitals.
hold on
orb1=1
orb2=4
for eloc=[1:5]
scatter([1:101]', enk(1:101,eloc)-EF,  (pnk(1:101,orb1,eloc))*1000+0.01,'b') 
scatter([1:101]', enk(1:101,eloc)-EF, (pnk(1:101,orb2,eloc))*1000,'g') 
hold on
end



axis([1 102 emin emax])



