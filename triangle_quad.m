function ML = triangle_quad(Bl,s,areal,p_4D,w_4D,xl,yl)

r_1 = max(0,(1- sum((Bl*[p_4D(:,1),p_4D(:,1).*(1-p_4D(:,2)+p_4D(:,2).*p_4D(:,3))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2))).*max(0,(1- sum((Bl*[p_4D(:,1).*(1-p_4D(:,2).*p_4D(:,3).*p_4D(:,4)),p_4D(:,1).*(1-p_4D(:,2))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2)))./((p_4D(:,1).*p_4D(:,2).*p_4D(:,3))'.^2.*sum((Bl*[p_4D(:,4)';ones(1,size(p_4D,1))]).^2));
r_2 = max(0,(1- sum((Bl*[p_4D(:,1),p_4D(:,1).*p_4D(:,2).*(1-p_4D(:,3)+p_4D(:,3).*p_4D(:,4))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2))).*max(0,(1- sum((Bl*[p_4D(:,1).*(1-p_4D(:,2).*p_4D(:,3)),p_4D(:,1).*p_4D(:,2).*(1-p_4D(:,3))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2)))./((p_4D(:,1).*p_4D(:,2).*p_4D(:,3))'.^2.*sum((Bl*[ones(1,size(p_4D,1));p_4D(:,4)']).^2));
r_3 = max(0,(1- sum((Bl*[p_4D(:,1),p_4D(:,1).*p_4D(:,2).*(1-p_4D(:,3))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2))).*max(0,(1- sum((Bl*[p_4D(:,1).*(1-p_4D(:,2).*p_4D(:,3).*p_4D(:,4)),p_4D(:,1).*(1-p_4D(:,3).*p_4D(:,4))]'+repmat([xl(1);yl(1)],1,size(p_4D,1))).^2)))./((p_4D(:,1).*p_4D(:,2).*p_4D(:,3))'.^2.*sum((Bl*[p_4D(:,4)';p_4D(:,4)'-1]).^2));

% f = @(t) t.^s./(t+1);
% 
% temp1 = arrayfun(@(x)[0:1/64:x],r_1,'UniformOutput',false);
% temp2 = arrayfun(@(x)[0:1/64:x],r_2,'UniformOutput',false);
% temp3 = arrayfun(@(x)[0:1/64:x],r_3,'UniformOutput',false);
% 
% f_1 = cellfun(f,temp1,'UniformOutput',false);
% f_2 = cellfun(f,temp2,'UniformOutput',false);
% f_3 = cellfun(f,temp3,'UniformOutput',false);
% 
% psi1 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_1,[1],ones(1,81)),f_1,'UniformOutput',false));
% psi2 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_2,[1],ones(1,81)),f_2,'UniformOutput',false));
% psi3 = cell2mat(cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),mat2cell(r_3,[1],ones(1,81)),f_3,'UniformOutput',false));

% psi1 = (r_1.^s)/s.*hypergeom([1,s],s+1,-r_1);
% psi2 = (r_2.^s)/s.*hypergeom([1,s],s+1,-r_2);
% psi3 = (r_3.^s)/s.*hypergeom([1,s],s+1,-r_3);
if s ==1/2
    psi1 = 2*atan(sqrt(r_1));
    psi2 = 2*atan(sqrt(r_2));
    psi3 = 2*atan(sqrt(r_3));
elseif s == 3/4
%     psi1 = (s^(-1)*r_1.^s).*hgeom2f1_eval(1,s,s+1,-r_1);
%     psi2 = (s^(-1)*r_2.^s).*hgeom2f1_eval(1,s,s+1,-r_2);
%     psi3 = (s^(-1)*r_3.^s).*hgeom2f1_eval(1,s,s+1,-r_3);
    
%      psi1 = log((1-sqrt(2)*r_1.^(1/4)+r_1.^(1/2))./(1+sqrt(2)*r_1.^(1/4)+r_1.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1+(r_1.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_1.^(1/4)/sqrt(2)),1-(r_1.^(1/4)/sqrt(2)));
%      psi2 = log((1-sqrt(2)*r_2.^(1/4)+r_2.^(1/2))./(1+sqrt(2)*r_2.^(1/4)+r_2.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1+(r_2.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_2.^(1/4)/sqrt(2)),1-(r_2.^(1/4)/sqrt(2)));
%      psi3 = log((1-sqrt(2)*r_3.^(1/4)+r_3.^(1/2))./(1+sqrt(2)*r_3.^(1/4)+r_3.^(1/2)))*sqrt(2)/2 + sqrt(2)*atan2((r_3.^(1/4)/sqrt(2)),1+(r_3.^(1/4)/sqrt(2)))+sqrt(2)*atan2((r_3.^(1/4)/sqrt(2)),1-(r_3.^(1/4)/sqrt(2)));

%      psi1 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_1)))))))));
%      psi2 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_2)))))))));
%      psi3 = 2*real(vpa(exp(vpa(pi/4*1i))*(atan(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_3))))))-atanh(vpa(exp(vpa(pi/4*1i))*sqrt(vpa(sqrt(vpa(r_3))))))))); 

      psi1 = 2*real(r_1.^(7/4)./(-r_1).^(7/4).*(atan((-r_1).^(1/4))-atanh((-r_1).^(1/4))));
      psi2 = 2*real(r_2.^(7/4)./(-r_2).^(7/4).*(atan((-r_2).^(1/4))-atanh((-r_2).^(1/4))));
      psi3 = 2*real(r_3.^(7/4)./(-r_3).^(7/4).*(atan((-r_3).^(1/4))-atanh((-r_3).^(1/4))));
elseif s == 1/4
    psi1 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_1)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_1)))));
    psi2 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_2)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_2)))));
    psi3 = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(r_3)))+atanh(exp(pi/4*1i)*sqrt(sqrt(r_3)))));    
elseif s == 7/10
    psi1 = -2*real(atan(r_1.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_1).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_1).^(1/10))+atanh((-1)^(7/10)*(r_1).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_1.^(1/10))));
    psi2 = -2*real(atan(r_2.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_2).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_2).^(1/10))+atanh((-1)^(7/10)*(r_2).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_2.^(1/10))));
    psi3 = -2*real(atan(r_3.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-r_3).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(r_3).^(1/10))+atanh((-1)^(7/10)*(r_3).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*r_3.^(1/10))));
else
    psi1 = r_1.^(s).*hgeom2f1_eval(1,s,s+1,-r_1)/s;
    psi2 = r_2.^(s).*hgeom2f1_eval(1,s,s+1,-r_2)/s;
    psi3 = r_3.^(s).*hgeom2f1_eval(1,s,s+1,-r_3)/s;
end
J = p_4D(:,1).^(2*s+1).*p_4D(:,2).^(2*s).*p_4D(:,3).^(2*s-1);

n1 = sum((Bl*[p_4D(:,4)';ones(1,size(p_4D,1))]).^2).^(s-1);
n2 = sum((Bl*[ones(1,size(p_4D,1));p_4D(:,4)']).^2).^(s-1);
n3 = sum((Bl*[p_4D(:,4)';p_4D(:,4)'-1]).^2).^(s-1);

% if ~isempty(find([J'.*n1.*psi1,J'.*n2.*psi2,J'.*n3.*psi3]>=10))
%     display('>100')
% end

ML = 8*areal*areal*sum(w_4D'.*J'.*(n1.*psi1+n2.*psi2+n3.*psi3));


end