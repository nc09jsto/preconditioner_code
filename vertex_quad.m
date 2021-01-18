function ML = vertex_quad(nodl,nodm,sh_nod,p,s,areal,aream,p_c,w_4D)
xm = p(1, nodm);
ym = p(2, nodm);
xl = p(1, nodl);
yl = p(2, nodl);

x = p_c(:,1);
y = p_c(:,2);
z = p_c(:,3);
w = p_c(:,4);

local_l = find(nodl==sh_nod);
nsh_l = find(nodl~=sh_nod);
nsh_m = find(nodm~=sh_nod);
p_c = [xl(local_l), yl(local_l)];
Bl = [xl(nsh_l(1))-p_c(1) xl(nsh_l(2))-xl(nsh_l(1)); yl(nsh_l(1))-p_c(2) yl(nsh_l(2))-yl(nsh_l(1))];
Bm = [xm(nsh_m(1))-p_c(1) xm(nsh_m(2))-xm(nsh_m(1)); ym(nsh_m(1))-p_c(2) ym(nsh_m(2))-ym(nsh_m(1))];

d1 = ((x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); y']-Bm*[z';(z.*w)']).^2));
d2 = ((x.^2)'.*sum((Bm*[ones(1,size(w_4D,1)); y']-Bl*[z';(z.*w)']).^2));


r_1 = max(0,1-(x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); y']+repmat(p_c',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[z';(z.*w)']+repmat(p_c',1,size(w_4D,1))).^2))./d1;
r_2 = max(0,1-(x.^2)'.*sum((Bm*[ones(1,size(w_4D,1)); y']+repmat(p_c',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bl*[z';(z.*w)']+repmat(p_c',1,size(w_4D,1))).^2))./d2;

% f = @(t) t.^s./(t+1);
% 
% temp1 = arrayfun(@(x)[0:1/64:x],r_1,'UniformOutput',false);
% temp2 = arrayfun(@(x)[0:1/64:x],r_2,'UniformOutput',false);
% 
% f_1 = cellfun(f,temp1,'UniformOutput',false);
% f_2 = cellfun(f,temp2,'UniformOutput',false);
% 
% psi1 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_1,[1],ones(1,81)),f_1,'UniformOutput',false));
% psi2 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_2,[1],ones(1,81)),f_2,'UniformOutput',false));

% psi1 = (r_1.^s)/s.*hypergeom([1,s],s+1,-r_1);
% psi2 = (r_2.^s)/s.*hypergeom([1,s],s+1,-r_2);
if s==1/2
    psi1 = 2*atan(sqrt(r_1));
    psi2 = 2*atan(sqrt(r_2));
elseif s==3/4
%     psi1 = (s^(-1)*r_1.^s).*hgeom2f1_eval(1,s,s+1,-r_1);
%     psi2 = (s^(-1)*r_2.^s).*hgeom2f1_eval(1,s,s+1,-r_2);
    
%       psi1 = log((1-sqrt(2)*r_1.^(1/4)+r_1.^(1/2))./(1+sqrt(2)*r_1.^(1/4)+r_1.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1+(r_1.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1-(r_1.^(1/4)/sqrt(2)));
%       psi2 = log((1-sqrt(2)*r_2.^(1/4)+r_2.^(1/2))./(1+sqrt(2)*r_2.^(1/4)+r_2.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1+(r_2.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1-(r_2.^(1/4)/sqrt(2)));

     %psi1 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1)))))))));
     %psi2 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2)))))))));

    psi1 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_1)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_1)))));
    psi2 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_2)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_2)))));

%       psi1 = 2*real(r_1.^(7/4)./(-r_1).^(7/4).*(atan((-r_1).^(1/4))-atanh((-r_1).^(1/4))));
%       psi2 = 2*real(r_2.^(7/4)./(-r_2).^(7/4).*(atan((-r_2).^(1/4))-atanh((-r_2).^(1/4))));
elseif s==1/4
    psi1 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_1)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_1)))));
    psi2 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_2)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_2)))));   
elseif s == 7/10
    psi1 = -2*real(atan(r_1.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_1).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_1).^(1/10))+atanh((-1)^(7/10)*(r_1).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_1.^(1/10))));
    psi2 = -2*real(atan(r_2.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_2).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_2).^(1/10))+atanh((-1)^(7/10)*(r_2).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_2.^(1/10))));
else
    psi1 = r_1.^(s).*hgeom2f1_eval(1,s,s+1,-r_1)/s;
    psi2 = r_2.^(s).*hgeom2f1_eval(1,s,s+1,-r_2)/s;
end
J = (x.^3.*z)';

% if ~isempty(find([J.*d1.^(s-1).*psi1,J.*d2.^(s-1).*psi2]>=10))
%     display('>100')
% end

ML = (4*areal*aream)*sum(w_4D'.*J.*(d1.^(s-1).*psi1+d2.^(s-1).*psi2));
end