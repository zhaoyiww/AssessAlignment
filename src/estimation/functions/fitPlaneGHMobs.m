function [plane] = fitPlaneGHMobs(P, sig, p0)

%--------------------------------------------------------------------------
% Plane estimation with Gauß-Helmert model
% Included variance component estimation (Niemeier 2008, chap. 9)
%
% INPUT:    Observations of polar measurements P = [s hz v];
%           Angular and range precisions for stochastic model sig (sig_w,
%           sig_s_const, sig_s_prob)

% OUTPUT:   plane [struct]
%           - corr - corrections(range, vertical angle, horizontal angle, dist to plane)
%           - prec - estimated precisions for measurement groups
%           - parameters + standard deviation
%
% Author: Berit Schmitz
% Date of latest modification: 01.02.2018
% Date of latest modification: 2018-06-16 by Jannik 
%--------------------------------------------------------------------------

%% ========================== Adjustment ==================================

% Initialization
test = -1;
count = 1;
asq_r = 1;
asq_v = 1;
asq_h = 1;


% % Iteration loop for adjusting the stochastic model
% while test == -1
    
    % Umrechnung der Beobachtungen in polar Koordinaten
    hz = P(:,2);
    v = P(:,3);
    s = P(:,1);
    
    % Vektorisierung der Beobachtungen
    l = zeros(length(hz(:,1))*3,1);
    l(3:3:end) = hz;
    l(2:3:end) = v;
    l(1:3:end) = s;
    
    n = length(l);      % Anzahl Beobachtungen
    m = length(P(:,1)); % Anzahl Punkte
    u = 3;              % Anzahl Unbekannte
    
%     if count > 1
        sll = zeros(n,1);
        sll(1:3:end) = asq_r *(sig.s_const + sig.s_prop.*s).^2;
        sll(2:3:end) = asq_v * sig.w^2;
        sll(3:3:end) = asq_h * sig.w^2;
        Sll = spdiags(sll,0,n,n);
%     else
%         % Kovarianzmatrix der Beobachtungen
%         sllvec = zeros(n,1);
%         sllvec(1:3:end) = sig_s^2;
%         sllvec(2:3:end) = sig_w^2;
%         sllvec(3:3:end) = sig_w^2;
%         Sll = spdiags(sllvec,0,n,n);
%     end
    
    % Zusammenfassung zu drei Ebenenparametern p = [nx_, ny_, nz_]
    % Bebobachtungsgleichung: nx_*r*sin(v)*cos(hz) + ny_*r*sin(v)*sin(hz) + nz_*r*cos(v) -1 = 0
    
%     % Näherungswerte Ebenenparameter
%     P_ = rhv2xyz(P);
%     a = P_(1,:)-P_(round(end/2),:);
%     b = P_(1,:)-P_(end,:);
%     
%     % Normalenvektor
%     pn = cross(a,b);
%     % Distanz d
%     d = pn(1)*s(1)*sin(v(1)).*cos(hz(1))+pn(2)*s(1)*sin(v(1))*sin(hz(1))+pn(3)*s(1)*cos(v(1));
%     
    pn = p0(1:3); d = p0(4);
    % Reduktion auf 3 Parameter
    p0 = pn./d;
    
    % Initialisieren
    V = zeros(n,1);
    dps = 100;
    anz = 0;
    
    l0 = l;
    % Iterieren
    while abs(dps) > 1e-10
        
        hz = l(3:3:end);
        v = l(2:3:end);
        s = l(1:3:end);
        % Designmatrix A
        A = [s.*sin(v).*cos(hz) s.*sin(v).*sin(hz) s.*cos(v)];
        
        % Bedingungsmatrix B
        gs = p0(1).*sin(v).*cos(hz)+p0(2)*sin(v).*sin(hz)+p0(3)*cos(v);
        gv = p0(1)*s.*cos(v).*cos(hz)+p0(2)*s.*cos(v).*sin(hz)-p0(3)*s.*sin(v);
        ghz = -p0(1)*s.*sin(v).*sin(hz)+p0(2)*s.*sin(v).*cos(hz);
        B = sparse([1:m;1:m;1:m],[1:3:n;2:3:n;3:3:n],[gs';gv';ghz'],m,n)';
        
        % Widersprueche w
        fx = p0(1)*s.*sin(v).*cos(hz)+p0(2)*s.*sin(v).*sin(hz)+p0(3)*s.*cos(v)-1;
        w = -B'*V+fx;
        
        % Normalgleichung loesen
        Sllq = B'*Sll*B;
        dlq = -w;
        
        Nq = A'/Sllq*A;
        nq = A'/Sllq*dlq;
        dps = Nq\nq;
        
        % Update
        p0 = p0+dps; % Update Parameter
        vq = dlq-A*dps;
        V = Sll*B/Sllq*vq; % Verbesserungen
        l = l0 + V;     %  Beobachtungen
        
        anz = anz+1; % Anzahl Iterationen
        
        
        
        
        if anz == 50
            disp('no solution found')
            break
        end
        
        
        
        
    end
    
%     % Iterative Gewichtung der Gewichtsmatrix bis Globaltest angenommen wird
%     [ T,k ] = globalTest( V,Sll,n-3,0.95 );
    
    % Verbesserungen der Streckenmessung, Vertikalwinkel, Horizontalwinkel
    vs = V(1:3:end);
    vv = V(2:3:end);
    vh = V(3:3:end);
    
    
    % If global test rejected -> compute a posteriori standard deviations
    
    
%     % Estimate variance components
%     % Cofactor matrix of residuals
%     Q = Sll;
%     H = speye(m,m)/(B'*Q*B);
%     Qxx = speye(size(Nq))/Nq;
%     
%     
%     qvv1 = diag(Q*B*H*B'*Q);
%     At = A';
%     Va = Q*B*H;
%     Vb = H*B'*Q;
%     qvv2 = zeros(n,1);
%     
%     [r,c] = find(Va~=0);
%     
%     for i = 1 : floor(n/100) : n
%         D = A(c(i),:)*Qxx*At;
%         qvv2(i,1) = Va(r(i),c(i))*D*Vb(:,i);
%     end
%     
%     
%     % Diagonal of cofactor matrix of corrections
%     qvv = qvv1-qvv2;
%     
%     
%     % Weight matrix for each measurement group
%     Pdiag = diag(speye(n,n)/Sll);
%     P_r = spdiags(Pdiag(1:3:end),0,m,m);
%     P_v = spdiags(Pdiag(2:3:end),0,m,m);
%     P_h = spdiags(Pdiag(3:3:end),0,m,m);
%     
%     % Partial redundancy for range measurements
%     r_r = sum(qvv(1:3:end).*Pdiag(1:3:end));
%     % Partial redundancy for vertical angle
%     r_v = sum(qvv(2:3:end).*Pdiag(2:3:end));
%     % Partial redundancy for range measurements
%     r_h = sum(qvv(3:3:end).*Pdiag(2:3:end));
%     
%     % Estimate variance components
%     sq_r = vs'*P_r*vs/r_r;
%     sq_v = vv'*P_v*vv/r_v;
%     sq_h = vh'*P_h*vh/r_h;
%     
%     % Update of variance component estimation
%     asq_r = asq_r * sq_r;
%     asq_v = asq_v * sq_v;
%     asq_h = asq_h * sq_h;
%     
%     
%     
%     % Stop conditions for iteration loop
%     if abs(sq_r-1)<0.1 && abs(sq_v-1)<0.1 && abs(sq_h-1)<0.1
%         test = 1;
%         continue;
%     end
%     
%     % Update counter
%     count = count + 1;
%     
%     if count == 20
%         warning('Not possible to adjust stochastic model!')
%         warning('Stochastic model requires improvement')
%         disp(['VEC_range = ' num2str(sq_r) ', VEC_vert = ' num2str(sq_v)...
%             ', VEC_hor = ' num2str(sq_h)])
%         test = 1;
%     end
%     
% end  % End of adjusting loop





%% =================== Final Results ======================================

% Ruecksubstitution auf 4 Parameter
psq = p0;
d = 1/norm(psq);
p = [psq.*d; d];



%% Standardabweichung der Parameter
Spp = inv(Nq);

Sppq = Spp;
nx = psq(1); ny = psq(2); nz = psq(3);
wurzel = sqrt(nx^2+ny^2+nz^2);

F(1,1) = (wurzel-nx*nx/wurzel)/wurzel^2;
F(1,2) = (      -nx*ny/wurzel)/wurzel^2;
F(1,3) = (      -nx*nz/wurzel)/wurzel^2;
F(2,1) = (      -ny*nx/wurzel)/wurzel^2;
F(2,2) = (wurzel-ny*ny/wurzel)/wurzel^2;
F(2,3) = (      -ny*nz/wurzel)/wurzel^2;
F(3,1) = (      -nz*nx/wurzel)/wurzel^2;
F(3,2) = (      -nz*ny/wurzel)/wurzel^2;
F(3,3) = (wurzel-nz*nz/wurzel)/wurzel^2;
F(4,1) = (      -   nx/wurzel)/wurzel^2;
F(4,2) = (      -   ny/wurzel)/wurzel^2;
F(4,3) = (      -   nz/wurzel)/wurzel^2;

Spp    = F*Sppq*F';
spp    = sqrt(diag(Spp));



% ----- Endgültige Residuen -----------------------------------------------

% Residuen erst nach Rücktransformation der Parameter berechnen!
% erstmal Widersprüche

w = p(1).*P(:,1).*sin(P(:,3)).*cos(P(:,2)) + ...
    p(2).*P(:,1).*sin(P(:,3)).*sin(P(:,2)) + ...
    p(3).*P(:,1).*cos(P(:,3))-p(4);


% Residuen sind negative Widersprüche

vq = -w;

% ------- Final precisions for measurement elements -----------------------

% After variance component estimation
a_r = sqrt(full(Sll(1,1)));
a_v = sqrt(full(Sll(2,2)));
a_h = sqrt(full(Sll(3,3)));


% ---------------- Store important results --------------------------------

% corrections (range, vertical angle, horizontal angle, dist to plane)
plane.corr.r = vs;
plane.corr.v = vv;
plane.corr.h = vh;
plane.corr.res = vq;

% estimated precisions
% plane.prec.r = a_r;
% plane.prec.v = a_v;
% plane.prec.h = a_h;

% estimated parameters (parameters + standard deviation)
plane.p = p;
plane.spp = spp;
plane.n = p(1:3);
plane.d = p(4);

% Number of iteration loops
plane.count = count;
plane.anz = anz;



end