function [R,p,t,bdrynodes,np,nt,nt_aux,nf,rawData] = import_mesh(filename)
%filename = 'UniformCircle1.txt'
rawData = importdata(filename);

row1 = rawData(1,:);
np   = row1(1); % #points
nt   = row1(2); % #triangles

p1 = rawData(2:np+1,1:3);  % points with indices
t1 = [rawData(np+2:2:np+1+2*nt,:),rawData(np+3:2:np+2+2*nt,1)]; % triangle (3 indices per row)
b = rawData(np+2+2*nt:end,1:3); % boundary edges (2 indices per row, cycles)


p = p1(:,1:2)';
t = t1(:,1:3);%[t1(t1(:,4)==1,1:3);t1(t1(:,4)==0,1:3)];
nt_aux = 0;%size(t1(t1(:,4)==0),1);
bdrynodes = b(b(:,3)==1,1);%unique(b(b(:,3) ==1,1:2));
nf = setdiff(reshape(t1(:,1:3),[],1),reshape(bdrynodes,[],1));

R = sqrt(max(p1(:,1).^2 + p1(:,2).^2));
