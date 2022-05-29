function [enk,pnk,Hrask,Hrrot] = wanbandsrot(hrdat,klist,rot,Rx,Ry,Rz)
% [enk,psink] = WANBANDS(hrdat,klist)
%
% Compuate LDOS for Wannier functions based on _hr.dat file, rotated by the
% rot matrix. 

% Header of file is read and skipped 
% let nw = # of Wannier functions
%
% Input vars:
%   hrdat = file name of _hr.dat file to read
%   klist = 3 x nk list of k vectors in lattice units
% Output vars:
%   enk = nk x nw eigenvalues
%   pnk = nk x nw x nw |psink|^2 probabilities 

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
invrot=inv(rot);
% # of kpoints and do the loop over k points
nk = size(klist,1);
enk = zeros(nk,nw);
pnk = zeros(nk,nw,nw);
for j=1:nk
  ka = klist(j,:)';
  kaR = d(:,1:3)*ka;
  HkRij = HRij.*exp(i*kaR);
  HkR = reshape(HkRij,nw,nw,nl/nw^2);
  Hk = squeeze(sum(HkR,3));
  Hk = invrot*(Hk + Hk')/2*rot;
  [vk,ek] = eig(Hk);
  vk=vk;
  ek = diag(ek)';
  enk(j,:) = ek;
  pnk(j,:,:) = (abs(vk).^2);
end
Hrask=[];
for scan=1:length(d(:,1))
  if (R(scan,:)==[Rx Ry Rz])
   Hrask=[Hrask, HRij(scan) ];
  end
end
Hrask=reshape(Hrask,nw,nw);
Hrask=real(Hrask);
Hrrot=invrot*Hrask*rot;
Hrrot=real(Hrrot);
end

