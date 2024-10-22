%To generate a k-point path from a n x n x n grid.

clear all
n=20
klist=0
for kz=0:n
    klist=klist+1;
    k(klist,1:3)=[0.5*(kz/n) 0.0 0.0];
    k(klist,4)=1/(3*n+1);
end
for ky=1:n
    klist=klist+1;
    k(klist,1:3)=[0.5  0.5*ky/n 0.0];
    k(klist,4)=1/(3*n+1);
end
for kx=1:n
    klist=klist+1;
    k(klist,1:3)=[0.5*(1-kx/n) 0.5*(1-kx/n) 0.0];
    k(klist,4)=1/(3*n+1);
end