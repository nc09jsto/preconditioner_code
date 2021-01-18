function [BCentre,ECentre,E,VA,EA,VE,TE,BT,BP,E2,Bbdrynodes,E_dual,N_dual,T_dual,EA2,VA2,TR2] = Barycentric_mesh(file)
[R,p,t,bdrynodes,nn,nt,nt_aux,nf,rawData] = import_mesh(file);

BCentre = 1/3*[p(1,t(:,1))+p(1,t(:,2))+p(1,t(:,3));p(2,t(:,1))+p(2,t(:,2))+p(2,t(:,3))];

TR = triangulation(t,p');

E = edges(TR);

ECentre = [p(1,E(:,1))+p(1,E(:,2));p(2,E(:,1))+p(2,E(:,2))]./2;

VA = vertexAttachments(TR);
EA = edgeAttachments(TR,E);

VE = cell(nn,1);
for i = 1:nn
[VE{i},~] = find(E(:,1)==i | E(:,2)==i);
end
VE = cellfun(@(x) x',VE,'UniformOutput', false);
TE = cell(size(E,1),1);
for i =1:size(E,1)
TE{i} = combvec([i],EA{i});
end
TE = cell2mat(cellfun(@(x) x',TE,'UniformOutput', false));
TE = sortrows(TE(:,2:-1:1))+nn;
TE(:,2) = TE(:,2)+nt;
BT = cell(nn,1);

temp1 = cellfun(@(x,y) combvec(x+nn,y+nn+nt)',VA,VE,'UniformOutput', false);
temp2 = cellfun(@(x) intersect(x,TE,'rows'),temp1,'UniformOutput', false);
for i =1:nn
    BT{i} = [i*ones(size(temp2{i},1),1),temp2{i}];
end
BT = cell2mat(BT);
BP = [p,BCentre,ECentre];

TR2 = triangulation(BT,BP');
E2 = edges(TR2);
Bbdrynodes = freeBoundary(TR2); 

E_dual = E2(sum(E2<=nn,2)==0,:);
N_dual = unique(E_dual);

T_dual = cellfun(@unique, temp1,'UniformOutput',false);

EA2 = edgeAttachments(TR2,E2);
VA2 = vertexAttachments(TR2);

%%
% trimesh(BT,BP(1,:),BP(2,:),zeros(size(BP,2),1),'FaceColor','none');hold on;
% trimesh(t,p(1,:),p(2,:),zeros(size(p,2),1),'FaceColor','none','EdgeColor','r'); hold on;plot(BCentre(1,:),BCentre(2,:),'o');view([0,90]); pbaspect([1 1 1]);plot(ECentre(1,:),ECentre(2,:),'ro');
% plot(p(1,:),p(2,:),'go');

% figure;
% for i = 1:size(E_dual,1)
%    plot(BP(1,E_dual(i,:)),BP(2,E_dual(i,:)),'b');hold on;
% end
% plot(BP(1,N_dual),BP(2,N_dual),'ro');
% plot(p(1,:),p(2,:),'go');

end



