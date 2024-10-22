function [enk,pnk,enkd,pnkd,Hkm,Hkdm] = wanbandsrotorbs(hrdat,klist,rot,on_bitals,off_bitals)

%This function takes the data from the Hamiltonian hrdat file, rotates this
%original matrix with the Hamiltonian rotation matrix (from RotationMatrix
%ThatIWantToUse.m script), then turns orbitals 'on' or 'off' depending on the set variables. 

% Compuate LDOS for Wannier functions based on _hr.dat file
% Header of file is read and skipped 
% let nw = # of Wannier functions
%
% Input vars:
%   hrdat = file name of _hr.dat file to read
%   klist = 3 x nk list of k vectors in lattice units
% Output vars:
%   enk = nk x nw eigenvalues
%   enkd = nk x nw eigenvalues after diagonalization
%   pnk = nk x nw x nw |psink|^2 probabilities 
%   pnkd = nk x nw x nw |psink|^2 probabilities after diagonalization
%   Hkm == Hrrask, originaly Hamiltonian, which is the hamiltonian of each
%   k point averaged over all kpoints-- related to enk and pnk
%   Hkdm == Hr but only the diagonal values (Without Ef change), which
%   takes the Hkm but only the diagonals, turning off the hopping--related
%   to enkd + pnkd

% open file and get the data
fprintf('Opening file %s\n',hrdat);
fid = fopen(hrdat);
txt = fgets(fid);
fprintf('hrdat: %s',txt);
tmp = textscan(fid,'%d',2);
nskip = tmp{1}(2);
fprintf('hrdat: %d %d  ---  skipping %d integers\n',tmp{1}(1),nskip,nskip);
nskip = tmp{1}(2);
tmp = textscan(fid,'%d',nskip);
dcell = textscan(fid,'%d %d %d %d %d %f %f');
fclose(fid);
clear fid tmp nskip txt

% parse the data: 
% # wannier -> nw , # lines -> nl , all data -> d
% R vectors -> R , H matrix elemnts -> HRij
nw = max(dcell{4});
nl = size(dcell{4},1);
d = zeros(nl,7);
for j = 1:7
    d(:,j) = dcell{j};
end
clear dcell j
R = d(:,1:3);
i = sqrt(-1);
HRij = d(:,6)+i*d(:,7);
invrot = inv(rot);

% # of kpoints and do the loop over k points
nk = size(klist,1);
enk = zeros(nk,nw);
pnk = zeros(nk,nw,nw);
enkd = zeros(nk,nw);
pnkd = zeros(nk,nw,nw);
Hkm = zeros(nw,nw);
Hkdm = zeros(nw,nw);

for j=1:nk
  ka = klist(j,:)';
  kaR = d(:,1:3)*ka;
  HkRij = HRij.*exp(i*kaR); %Fourier Transfroms from R space into k space
  HkR = reshape(HkRij,nw,nw,nl/nw^2); 
  Hk = squeeze(sum(HkR,3)); %Reshapes data entries to make Hamiltonian at a kpoint in the kgrid
  % Hk = invrot*(Hk + Hk')/2*rot; %Rotates Hamiltonian as though it were a Jacobian
  Hk=((invrot*Hk*rot)+(invrot*Hk*rot)')/2;
  [vk,ek] = eig(Hk); %Generates eigenvectors/eigenvalues of the Hamiltonian
  ek = diag(ek)';
  enk(j,:) = ek; %Generates eigenvalue matrix (nk x nw)
  pnk(j,:,:) = abs(vk).^2; %Generates probability matrix (nk x nw x nw)
  Hkm = Hkm + Hk; %builds final Hrrask hamiltonian (now with self rotation)

  Hkdiag = diag(diag(Hk)); %Just the on-site energies right now
  Nodiags = Hk - Hkdiag; %everything except the onsite energies
  numkeep = length(on_bitals);
  KeepPile = zeros(nw,nw);
  if off_bitals ==0 %If all orbitals are on, then don't change anything from Hamiltonian.
      Hkdiag = Hk;
  elseif on_bitals ==0 %If all orbitals are off, then only have the diagonals.
      Hkdiag = Hkdiag;

%Key Step! Depending on if you'r purpose is to turn an orbital off, or to
%turn an orbital on, Change the Final Else Statement to reflect intentions!

  % else %For each orbital left on, keep the total row/column of the orbital
  %     for b = 1:(numkeep)
  %         K = on_bitals(b);
  %         KeepPile(K,:) = Nodiags(K,:); 
  %         KeepPile(:,K) = Nodiags(:,K);
  %     end
  %     Hkdiag = Hkdiag + KeepPile;
  % 
  else             %Turn this on if you are more focused in turning something off as opposed to keeping orbitals on
      for b = 1:length(off_bitals)
          K = off_bitals(b);
          Nodiags(K,:) = zeros(1,nw);
          Nodiags(:,K) = zeros(nw,1);
      end
      Hkdiag = Hkdiag + Nodiags;

  end
  [vkd, ekd] = eig(Hkdiag);
  ekd = diag(ekd)';
  enkd(j,:) = ekd; %Rotated eigenvalues matrix
  pnkd(j,:,:) = abs(vkd).^2; %Rotates probability matrix
  Hkdm = Hkdm + Hkdiag;
end
Hkm = Hkm / nk;
Hkm = real(Hkm)
Hkdm = Hkdm / nk;
Hkdm = real(Hkdm)
end
