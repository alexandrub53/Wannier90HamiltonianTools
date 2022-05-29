function [ norb,yorb,ntot,ytot ,rho,rhoe,rhov,Hrask,Hrrot] = GetPdosAndNFromWan(rotation, hamfilename, kgrid, orbitals , EF, smearing, xout,soc,spin,R)
btype=2;
griddensity=kgrid(1);
griddensityz=kgrid(3);
kFplot=[]

for kz=1:griddensityz
for kx=1:griddensity
    for ky=1:griddensity
     kFplot=[kFplot; 2*pi*(kx-griddensity/2)/griddensity  2*pi*(ky-griddensity/2)/griddensity  2*pi*(kz-griddensityz/2)/griddensityz];
    end
end
end


[rho,rhoe,rhov] = GetRho(hamfilename,kFplot,EF,spin);  


numE=size(rhoe,1);  
Rx=R(1,1);
Ry=R(1,2);
Rz=R(1,3);

if rotation==0
   rhov = eye(size(rhoe,1));
end


[enkr,pnkr,Hrask,Hrrot] = wanbandsrot(hamfilename,kFplot,rhov,Rx,Ry,Rz);

n=zeros(size(orbitals,2));
for l=1:size(orbitals,2)
  yout=zeros(size(xout));
  for eloc=1:numE
    yout =yout+ broaden(enkr(:,eloc)-EF,pnkr(:,orbitals(l),eloc),xout,smearing,btype);
  end
norb(l)=trapz(xout,(yout.*(1-sign(xout))/2))/(trapz(xout,yout));
yorb(:,l)=yout/(trapz(xout,yout));
end

%if there's no spin-orbit coupling, and the calculation is non-magnetic, we add the density over the spins.
if soc==0
    norb=norb*2;
    yorb=yorb*2;
end
ntot=sum(norb);
ytot=sum(yorb,2);



end

