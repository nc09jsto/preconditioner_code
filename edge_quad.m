function ML = edge_quad(nodl,nodm,nod_diff,p,s,areal,aream,p_c,w_4D)
xm = p(1, nodm); 
ym = p(2, nodm);

xl = p(1, nodl); 
yl = p(2, nodl);

x = p_c(:,1); 
y = p_c(:,2); 
z = p_c(:,3);
w = p_c(:,4);

local_l = find(nodl~=nod_diff(1)); 
nsh_l = find(nodl==nod_diff(1)); 
nsh_m = find(nodm==nod_diff(2)); 

P1 = [xl(local_l(1)), yl(local_l(1))]; 
P2 = [xl(local_l(2)), yl(local_l(2))]; 

Bl = [P2(1)-P1(1) -P2(1)+xl(nsh_l); P2(2)-P1(2) -P2(2)+yl(nsh_l)];
Bm = [P2(1)-P1(1) -P2(1)+xm(nsh_m); P2(2)-P1(2) -P2(2)+ym(nsh_m)];

d1 = ((x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); (y.*w)']-Bm*[(1-y.*z)';(y.*(1-z))']).^2));
d2 = ((x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); y']-Bm*[(1-y.*z.*w)';(y.*z.*(1-w))']).^2));
d3 = ((x.^2)'.*sum((Bl*[(1-y.*z)'; (y.*(1-z))']-Bm*[ones(1,size(w_4D,1));(y.*z.*w)']).^2));
d4 = ((x.^2)'.*sum((Bl*[(1-y.*z.*w)'; (y.*z.*(1-w))']-Bm*[ones(1,size(w_4D,1));y']).^2));
d5 = ((x.^2)'.*sum((Bl*[(1-y.*z.*w)'; (y.*(1-z.*w))']-Bm*[ones(1,size(w_4D,1));(y.*z)']).^2));


r_1 = max(0,1-(x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); (y.*w)']+repmat(P1',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[(1-y.*z)';(y.*(1-z))']+repmat(P1',1,size(w_4D,1))).^2))./d1;
r_2 = max(0,1-(x.^2)'.*sum((Bl*[ones(1,size(w_4D,1)); y']+repmat(P1',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[(1-y.*z.*w)';(y.*z.*(1-w))']+repmat(P1',1,size(w_4D,1))).^2))./d2;
r_3 = max(0,1-(x.^2)'.*sum((Bl*[(1-y.*z)'; (y.*(1-z))']+repmat(P1',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[ones(1,size(w_4D,1));(y.*z.*w)']+repmat(P1',1,size(w_4D,1))).^2))./d3;
r_4 = max(0,1-(x.^2)'.*sum((Bl*[(1-y.*z.*w)'; (y.*z.*(1-w))']+repmat(P1',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[ones(1,size(w_4D,1));y']+repmat(P1',1,size(w_4D,1))).^2))./d4;
r_5 = max(0,1-(x.^2)'.*sum((Bl*[(1-y.*z.*w)'; (y.*(1-z.*w))']+repmat(P1',1,size(w_4D,1))).^2)).*max(0,1-(x.^2)'.*sum((Bm*[ones(1,size(w_4D,1));(y.*z)']+repmat(P1',1,size(w_4D,1))).^2))./d5;


% f = @(t) t.^s./(t+1);
% 
% psi1 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_1,[1],ones(1,81)),cellfun(f,arrayfun(@(x)[0:1/64:x],r_1,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
% psi2 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_2,[1],ones(1,81)),cellfun(f,arrayfun(@(x)[0:1/64:x],r_2,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
% psi3 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_3,[1],ones(1,81)),cellfun(f,arrayfun(@(x)[0:1/64:x],r_3,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
% psi4 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_4,[1],ones(1,81)),cellfun(f,arrayfun(@(x)[0:1/64:x],r_4,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));
% psi5 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_5,[1],ones(1,81)),cellfun(f,arrayfun(@(x)[0:1/64:x],r_5,'UniformOutput',false),'UniformOutput',false),'UniformOutput',false));

% psi1 = (r_1.^s)/s.*hypergeom([1,s],s+1,-r_1);
% psi2 = (r_2.^s)/s.*hypergeom([1,s],s+1,-r_2);
% psi3 = (r_3.^s)/s.*hypergeom([1,s],s+1,-r_3);
% psi4 = (r_4.^s)/s.*hypergeom([1,s],s+1,-r_4);
% psi5 = (r_5.^s)/s.*hypergeom([1,s],s+1,-r_5);
if s==1/2
    psi1 = 2*atan(sqrt(r_1));
    psi2 = 2*atan(sqrt(r_2));
    psi3 = 2*atan(sqrt(r_3));
    psi4 = 2*atan(sqrt(r_4));
    psi5 = 2*atan(sqrt(r_5));
elseif s==3/4
%     psi1 = (s^(-1)*r_1.^s).*hgeom2f1_eval(1,s,s+1,-r_1);
%     psi2 = (s^(-1)*r_2.^s).*hgeom2f1_eval(1,s,s+1,-r_2);
%     psi3 = (s^(-1)*r_3.^s).*hgeom2f1_eval(1,s,s+1,-r_3);
%     psi4 = (s^(-1)*r_4.^s).*hgeom2f1_eval(1,s,s+1,-r_4);
%     psi5 = (s^(-1)*r_5.^s).*hgeom2f1_eval(1,s,s+1,-r_5);

%      psi1 = log((1-sqrt(2)*r_1.^(1/4)+r_1.^(1/2))./(1+sqrt(2)*r_1.^(1/4)+r_1.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1+(r_1.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1-(r_1.^(1/4)/sqrt(2)));
%      psi2 = log((1-sqrt(2)*r_2.^(1/4)+r_2.^(1/2))./(1+sqrt(2)*r_2.^(1/4)+r_2.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1+(r_2.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1-(r_2.^(1/4)/sqrt(2)));
%      psi3 = log((1-sqrt(2)*r_3.^(1/4)+r_3.^(1/2))./(1+sqrt(2)*r_3.^(1/4)+r_3.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_3.^(1/4)/sqrt(2)),1+(r_3.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_3.^(1/4)/sqrt(2)),1-(r_3.^(1/4)/sqrt(2)));
%      psi4 = log((1-sqrt(2)*r_4.^(1/4)+r_4.^(1/2))./(1+sqrt(2)*r_4.^(1/4)+r_4.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_4.^(1/4)/sqrt(2)),1+(r_4.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_4.^(1/4)/sqrt(2)),1-(r_4.^(1/4)/sqrt(2)));
%      psi5 = log((1-sqrt(2)*r_5.^(1/4)+r_5.^(1/2))./(1+sqrt(2)*r_5.^(1/4)+r_5.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_5.^(1/4)/sqrt(2)),1+(r_5.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_5.^(1/4)/sqrt(2)),1-(r_5.^(1/4)/sqrt(2)));
% 
%     psi1 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1)))))))));
%     psi2 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2)))))))));
%     psi3 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_3))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_3)))))))));
%     psi4 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_4))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_4)))))))));
%     psi5 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_5))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_5)))))))));
 
     psi1 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_1)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_1)))));
     psi2 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_2)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_2)))));
     psi3 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_3)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_3)))));
     psi4 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_4)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_4)))));
     psi5 = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_5)))-atanh(exp(pi/4*1i)*sqrt(sqrt(r_5)))));

%       psi1 = 2*real(r_1.^(7/4)./(-r_1).^(7/4).*(atan((-r_1).^(1/4))-atanh((-r_1).^(1/4))));
%       psi2 = 2*real(r_2.^(7/4)./(-r_2).^(7/4).*(atan((-r_2).^(1/4))-atanh((-r_2).^(1/4))));
%       psi3 = 2*real(r_3.^(7/4)./(-r_3).^(7/4).*(atan((-r_3).^(1/4))-atanh((-r_3).^(1/4))));
%       psi4 = 2*real(r_4.^(7/4)./(-r_4).^(7/4).*(atan((-r_4).^(1/4))-atanh((-r_4).^(1/4))));
%       psi5 = 2*real(r_5.^(7/4)./(-r_5).^(7/4).*(atan((-r_5).^(1/4))-atanh((-r_5).^(1/4))));
elseif s==1/4
    psi1 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_1)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_1)))));
    psi2 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_2)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_2)))));
    psi3 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_3)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_3)))));
    psi4 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_4)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_4)))));
    psi5 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_5)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_5)))));
elseif s == 7/10
    psi1 = -2*real(atan(r_1.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_1).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_1).^(1/10))+atanh((-1)^(7/10)*(r_1).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_1.^(1/10))));
    psi2 = -2*real(atan(r_2.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_2).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_2).^(1/10))+atanh((-1)^(7/10)*(r_2).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_2.^(1/10))));
    psi3 = -2*real(atan(r_3.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_3).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_3).^(1/10))+atanh((-1)^(7/10)*(r_3).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_3.^(1/10))));
    psi4 = -2*real(atan(r_4.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_4).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_4).^(1/10))+atanh((-1)^(7/10)*(r_4).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_4.^(1/10))));
    psi5 = -2*real(atan(r_5.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_5).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_5).^(1/10))+atanh((-1)^(7/10)*(r_5).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_5.^(1/10))));
else
    psi1 = r_1.^(s).*hgeom2f1_eval(1,s,s+1,-r_1)/s;
    psi2 = r_2.^(s).*hgeom2f1_eval(1,s,s+1,-r_2)/s;
    psi3 = r_3.^(s).*hgeom2f1_eval(1,s,s+1,-r_3)/s;
    psi4 = r_4.^(s).*hgeom2f1_eval(1,s,s+1,-r_4)/s;
    psi5 = r_5.^(s).*hgeom2f1_eval(1,s,s+1,-r_5)/s;
end
J1 = (x.^3.*y.^2)';
J2 = (x.^3.*y.^2.*z)';
% 
% if ~isempty(find([J1.*d1.^(s-1).*psi1,J2.*d2.^(s-1).*psi2,J2.*d3.^(s-1).*psi3,J2.*d4.^(s-1).*psi4,J2.*d5.^(s-1).*psi5]>=10))
%     display('>100')
% end

ML = 4*areal*aream*sum(w_4D'.*(J1.*d1.^(s-1).*psi1+J2.*(d2.^(s-1).*psi2+d3.^(s-1).*psi3+d4.^(s-1).*psi4+d5.^(s-1).*psi5)));
end