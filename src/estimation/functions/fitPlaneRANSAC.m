function [RSC] = fitPlaneRANSAC(xyz,T,maxIter, minDist)
% input:
%       xyz...  [nx3] 3D-Punktewolke
%       T...    [double] Schwellwert: maximaler Abstand von der Ebene, der für Inlier erlaubt ist
%       maxIter... [int] maximale Anzahl an Iterationen
%       minDist... [double] mindest Abstand Punkte zur Ebenenschätzung

if nargin<4
    minDist = 0.05;
end

% Hilfsgrößen
m = size(xyz,1);

%Hilfsgrößen für Schleife
j = 1;
Ps = zeros(maxIter,4);
numinlier = zeros(maxIter,1);
Ind = zeros(maxIter,3);
while j <= maxIter
    % Auswahl der drei Punkt, Bedinung: sie müssen min. 1/4 der Targetgröße
    % oder 5 cm auseinander liegen
    
    % Hilfgrößen
    i = 1;
    d = 0;
    %while d < minDist
        
        % Auswahl
        ind = round(rand(3,1)*(m-1))+1;
        
        x12 = xyz(ind(2),:)'-xyz(ind(1),:)';
        x13 = xyz(ind(3),:)'-xyz(ind(1),:)';
        %x23 = xyz(ind(3),:)'-xyz(ind(2),:)';
        
        %d1 = norm(x12);
        %d2 = norm(x13);
        %d3 = norm(x23);
        %d = min([d1,d2,d3]);
        %i = i+1;
        %if i == 200;
         %   break
        %end
    %end
    % Normalenvektor berechenn
    n = cross(x12,x13);
    n = n/norm(n);
    d = n'*xyz(ind(1),:)';
    % Verbesserung zu allen anderen Punkten
    v = n(1)*xyz(:,1)+n(2)*xyz(:,2)+n(3)*xyz(:,3)-d;
    
    % Speichern der Werte für außerhalb der Schleife
    numinlier(j,1) = sum(abs(v)<=T);
    Ps(j,:) = [n',d];
    Ind(j,:) = ind';
    j = j+1;
end

% Finden der besten Lösung
best = find(numinlier == max(numinlier));
ps = Ps(best(1),:)'; % beste Parameter
numinlier_best = numinlier(best(1));
ind = Ind(best(1),:)';

v  = ps(1)*xyz(:,1)+ps(2)*xyz(:,2)+ps(3)*xyz(:,3)-ps(4);
idxInlier = (abs(v)<=T);
inlier = xyz(idxInlier);

RSC.ps          = ps;           % beste Parameter
RSC.v           = v;            % beste Verbesserungsquadratsumme
RSC.idxInlier   = idxInlier;    % Indizies der Inlier
RSC.inlier      = inlier;       % Koordianten, der Inlier
RSC.numinlier_best = numinlier_best;  % Anzahl der Inlier für die beste Lösng
RSC.numinlier_all  = numinlier;       % Anzahl der Inlier für jede Iteration
RSC.ind         = ind;          % Indizies der besten Punkte
RSC.best        = best(1);      % Indize der besten drei Punkt
RSC.Ps          = Ps;           % alle Parameter
RSC.Ind         = Ind;          % alle gewählten Indizies

end

