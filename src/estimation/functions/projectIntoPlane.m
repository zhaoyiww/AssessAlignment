function rhv = projectIntoPlane(rhv, sig, p)

% Beobachtungen
s = rhv(:,1);
hz = rhv(:,2);
v = rhv(:,3);
% Vektorisierung der Beobachtungen
l = zeros(length(hz(:,1))*3,1);
l(3:3:end) = hz;
l(2:3:end) = v;
l(1:3:end) = s;

n = length(l);      % Anzahl Beobachtungen
m = length(rhv(:,1)); % Anzahl Punkte
u = 3;              % Anzahl Unbekannte

% Kovarianzmatrix
sll = zeros(n,1);
sll(1:3:end) = (sig.s_const + sig.s_prop.*s).^2;
sll(2:3:end) = sig.w^2;
sll(3:3:end) = sig.w^2;
Sll = spdiags(sll,0,n,n);

% Reduktion auf 3 Parameter
pn = p(1:3); d = p(4);
p = pn./d;

% Bedingungsmatrix B
gs = p(1).*sin(v).*cos(hz)+p(2)*sin(v).*sin(hz)+p(3)*cos(v);
gv = p(1)*s.*cos(v).*cos(hz)+p(2)*s.*cos(v).*sin(hz)-p(3)*s.*sin(v);
ghz = -p(1)*s.*sin(v).*sin(hz)+p(2)*s.*sin(v).*cos(hz);
B = sparse([1:m;1:m;1:m],[1:3:n;2:3:n;3:3:n],[gs';gv';ghz'],m,n)';

% Widersprueche w
w = p(1)*s.*sin(v).*cos(hz)+p(2)*s.*sin(v).*sin(hz)+p(3)*s.*cos(v)-1;

% Normalgleichung loesen
Sllq = B'*Sll*B;

% Update
V = Sll*B/Sllq*-w; % Verbesserungen
ls = l + V;     %  Beobachtungen

rhv(:,1) = ls(1:3:end);
rhv(:,2) = ls(3:3:end);
rhv(:,3) = ls(2:3:end);




