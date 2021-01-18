function [h] = hgeom2f1_eval(a,b,c,z)
h = zeros(size(z));
indLess = find(abs(z)<1);
ind1 = find(abs(z)==1);
indMore = find(abs(z)>1);

h(indLess) = hyper2f1(a,b,c,z(indLess),1e-20);
h(ind1) = hypergeom([a,b],c,z(ind1));
h(indMore) = (1-z(indMore)).^(-a).*hyper2f1(a,c-b,c,z(indMore)./(z(indMore)-1),1e-20);
end