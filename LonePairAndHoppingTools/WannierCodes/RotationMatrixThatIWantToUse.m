function [Hrask,rote,rotv] = RotationMatrixThatIWantToUse(hrdat,Rx,Ry,Rz)

%This function creates the Rotation Hamiltonian Matrix, which for us is the
%Hamiltonian at R =[0,0,0], from the original hrdat file.

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

Hrask=[];
for scan=1:length(d(:,1))
  if (R(scan,:)==[Rx Ry Rz])
   Hrask=[Hrask, HRij(scan) ];
  end
end
Hrask=reshape(Hrask,nw,nw);
Hrask=real(Hrask);

[Hrv,HrE] = eig(Hrask);
HrE = diag(HrE);

rotv = Hrv;
rote = HrE;

end
