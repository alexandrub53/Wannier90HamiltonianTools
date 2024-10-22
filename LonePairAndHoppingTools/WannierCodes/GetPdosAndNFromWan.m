%Code to calculate the necessary information for PDOS and bands functions,
%with additional calculation for the electronic density in each orbital
%before/after rotation.

function [ norb,yorb,eorb,ntot,ytot,rot,rote,rotv,enk,enkd,pnk,pnkd,Hkm,Hkdm] = GetPdosAndNFromWan(rotation, hamfilename, kgrid, orbitals , EF, smearing, xout,soc,spin,R,on_bitals,off_bitals)
btype=2;  %Type of peak broadening (2 = Gaussian)
griddensity=kgrid(1); %Grid density refers to the kgrid mesh to be used
griddensityz=kgrid(3);
kFplot=[]

for kz=1:griddensityz
for kx=1:griddensity
    for ky=1:griddensity
     kFplot=[kFplot; 2*pi*(kx-griddensity/2)/griddensity  2*pi*(ky-griddensity/2)/griddensity  2*pi*(kz-griddensityz/2)/griddensityz];
    end
end
end

[rot,rote,rotv] = RotationMatrixThatIWantToUse(hamfilename,0,0,0); %Call the rotation matrix from the script, if you are using Hamiltonian at R=[0,0,0]. If you are using Rho rotation, switch script callings.
% [rot,rote,rotv] = GetRho(hamfilename,kFplot,EF,spin);  

%Fixing up to be able to rotatable. Rotating each atom on its own.
[V,~]=eig(rot(1:4,1:4));
[V2,~]=eig(rot(5:8,5:8));
[V3,~]=eig(rot(9:11,9:11));
[V4,~]=eig(rot(12:14,12:14));
[V5,~]=eig(rot(15:17,15:17));
[V6,~]=eig(rot(18:20,18:20));

rotv=zeros(20);
rotv(1:4,1:4)=rotv(1:4,1:4)+V;
rotv(5:8,5:8)=rotv(5:8,5:8)+V2;
rotv(9:11,9:11)=rotv(9:11,9:11)+V3;
rotv(12:14,12:14)=rotv(12:14,12:14)+V4;
rotv(15:17,15:17)=rotv(15:17,15:17)+V5;
rotv(18:20,18:20)=rotv(18:20,18:20)+V6;

numE=size(rote,1);  
Rx=R(1,1); %Re-call the R you want
Ry=R(1,2);
Rz=R(1,3);

if rotation==0 %Just in case you don't want to rotate anything.
   rotv = eye(size(rhoe,1));
end

[enk,pnk,enkd,pnkd,Hkm,Hkdm] = wanbandsrotorbs(hamfilename,kFplot,rotv,on_bitals,off_bitals); %Retrieve the rotated/diagonalized Hamiltonians. 

n=zeros(size(orbitals,2));
for l=1:size(orbitals,2)
  yout=zeros(size(xout));
  for eloc=1:numE
    yout =yout+ broaden(enkd(:,eloc)-EF,pnkd(:,orbitals(l),eloc),xout,smearing,btype); %Creates the pdos peaks with Gaussian broadening using the 'Broaden.m' script. Uses the diagonalized eigen values.
  end
norb(l)=trapz(xout,(yout.*(1-sign(xout))/2))/(trapz(xout,yout)); %norb = number of electrons in orbital, using trapezoidal integrals
eorb(l)=trapz(xout,(yout.*xout))/(trapz(xout,yout)); %eorb = average energy in orbital, using trapezoidal integrals
yorb(:,l)=yout/(trapz(xout,yout)); %yorb = pdos of orbital, also uses trapezoidal integrals
end

%if there's no spin-orbit coupling, and the calculation is non-magnetic, we add the density over the spins.
if soc==0
    norb=norb*2;
    yorb=yorb*2;
end
ntot=sum(norb); %total number of electrons
ytot=sum(yorb,2); %total pdos



end
