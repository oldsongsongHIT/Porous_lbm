lx=6;
ly=1;
n=40;
nx=lx*n+1;
ny=ly*n+1;
obst=zeros(nx,ny);
cx=n*lx/6+1;
cy=n*ly/2+1;
r=10;
for i=1:nx
    for j=1:ny
        if (i-cx)^2+(j-cy)^2<r^2
            obst(i,j)=1;
        end
    end
end
filename=['geometry.dat'];
fid=fopen(filename,'wt');
for i=1:nx
    fprintf(fid, '%g ', obst(i,:));
end
fclose(fid);
            
