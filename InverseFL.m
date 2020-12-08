% Inverse operator b
tic
s = 0.5;
mesh_file = 'UniCircle1.msh';

run helper.m;   % Precomputed quadrature points

kns = 2^(-2*s)/(gamma(s))^2/pi;

[R,p,t,bdrynodes,np,nt,nt_aux,nf,rawData] = import_mesh(mesh_file);
[BCentre,ECentre,E,VA,EA,VE,TE,BT,BP,E2,Bbdrynodes,E_dual,N_dual,T_dual,EA2,VA2,TR2] = Barycentric_mesh(mesh_file);

B_nn = size(BP,2); 
B_nt = size(BT,1);
B = zeros(B_nt,B_nt);
B_area = zeros(B_nt,1);

for i=1:B_nt
    a = BP( : ,BT(i,:));
    B_area(i) = 0.5.*abs(det([ a(:,1) - a(:,3) ...
        a(:,2) - a(:,3)]));
end

aux_ind = reshape( repmat( 1:3:3*B_nt , size(p_T_6,1) , 1 ) , [] , 1 );
% empty_vtx = zeros(2,3*B_nt);
% BBm = zeros(2,2*B_nt);

% Main loop
for l=1:B_nt
    %l
    edge = [ VA2{BT(l,1)} VA2{BT(l,2)} VA2{BT(l,3)} ];
    [nonempty, M, N] = unique( edge , 'first' );
    edge(M) = [];
    vertex = setdiff( nonempty , edge );
    ll = B_nt - l + 1 - sum( nonempty>=l );
    edge( edge<=l ) = [];
    vertex( vertex<=l ) = [];
    empty = setdiff( l:B_nt , nonempty );
    empty_vtx = BP( : , BT( empty(1:ll) , : )' );
    nodl = BT(l,:); 
    xl = BP(1 , nodl); yl = BP(2 , nodl);
    Bl = [xl(2)-xl(1) yl(2)-yl(1); xl(3)-xl(2) yl(3)-yl(2)]';
    %identical elements
    B(l, l) = B(l, l) + triangle_quad(Bl,s,B_area(l),p_4D,w_4D,xl,yl)/2; 
    
    if ~isempty(empty)

    BBm = reshape( [ empty_vtx( : , 2:3:3*ll )...
        - empty_vtx( : , 1:3:3*ll ) , empty_vtx( : , 3:3:3*ll )...
        - empty_vtx( : , 2:3:3*ll ) ] , [] , 2)' ;
    %vl = p_T_6*(Bl') + [ones(6,1).*xl(1) ones(6,1).*yl(1)];
    vl = Bl*p_T_6' + repmat([xl(1);yl(1)],1,size(p_T_6,1));
    vm = reshape(permute(reshape(p_T_6*BBm(:,1:2*ll),[size(p_T_6,1) 1 2 ll]),[1 4 3 2]),[size(p_T_6,1)*ll 2])...
        + empty_vtx(:,aux_ind(1:size(p_T_6,1)*ll))';
    
    norms = reshape(bsxfun(@plus,sum(vl'.^2,2),sum(vm.^2,2)')-2*vl'*vm', size(p_T_6,1)^2 , [] );
    numerator = reshape( max(1-sum(vl'.^2,2),0)*max((1-sum(vm.^2,2))',0), size(p_T_6,1)^2 , []);

    rs = numerator./(norms);
    if s==1/2
    Pint_final = 2*atan(sqrt(rs));
    elseif s == 3/4
        Pint_final = 2*real(exp(pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(rs)))-atanh(exp(pi/4*1i)*sqrt(sqrt(rs)))));
    elseif s == 1/4
        Pint_final = 2*real(-exp(3*pi/4*1i)*(atan(exp(pi/4*1i)*sqrt(sqrt(rs)))+atanh(exp(pi/4*1i)*sqrt(sqrt(rs)))));
    elseif s == 7/10
        Pint_final = -2*real(atan(rs.^(1/10))+(-1)^(1/10)*((-1)^(1/5)*atanh((-rs).^(1/10))+(-1)^(4/5)*atanh((-1)^(3/10)*(rs).^(1/10))+atanh((-1)^(7/10)*(rs).^(1/10))+ (-1)^(3/5)*atanh((-1)^(9/10)*rs.^(1/10))));
    else
%         Psi_int = @(t) t.^s./(t+1);
%         PINT = cellfun(@(x) Psi_int(0:1/64:x),rs,'UniformOutput',false);
%         Pint_final = cellfun(@(x,y) (x^s)/s - (x/(3*length(y)-1))*(y(1)+2*sum(y(3:2:end-2))+4*sum(y(2:2:end))+y(end)),rs,PINT,'UniformOutput',false);
          Pint_final = rs.^(s).*hgeom2f1_eval(1,s,s+1,-rs)/s;
    end
    B(l,empty) = 8*B_area(l)*B_area((empty(1:ll)))'.*((phiB*(Pint_final.*norms(:,1:ll).^(s-1))));
    end

    %vertex sharing
    for m = vertex
        nodm = BT(m,:);
        nod_com = intersect(nodl, nodm);
        B(l,m) = B(l,m) + 2.*vertex_quad(nodl,nodm,nod_com,BP,s,B_area(l),B_area(m),p_4D,w_4D);
    end

    %edge sharing
    for m = edge
        nodm=BT(m,:);
        nod_diff = [setdiff(nodl, nodm) setdiff(nodm, nodl)];
        B(l,m) = B(l,m)+ 2.*edge_quad(nodl,nodm,nod_diff,BP,s,B_area(l),B_area(m),p_4D,w_4D);
    end
end
B = kns*(B+B');
toc;