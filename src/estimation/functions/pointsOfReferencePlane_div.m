function [xyz,i,n,d,beta_ap,s_ap] = pointsOfReferencePlane_div(xyz,i, refp_ap, sot, sig)
% intput:   xyz...  Koordinaten
%           i...    Intensitäten
%           ap...   Näherungswert [1x3] karthesischer Punkt (xyz)
%           sot...  (size of Target) Ungefähre Größe des ZZ in m
%           sig...  Stochastisches Modell der Einzelpunktmessung
%
% output:   xyz...  Punkte, die zur Ebene gehören.
%
% modified: 2017-08-04
% modified: 2017-11-28 efficiency and threshold for outlier detection by
% plane
% modified: 2018-02-03
% modified: 2018-06-18 changing inputs, new threshold
%
% author: Jannik Janßen

%% Punkte mit gewissen Distanz zum Näherungswert verwerfen
d = sqrt((xyz(:,1)-refp_ap(:,1)).^2+(xyz(:,2)-refp_ap(:,2)).^2+(xyz(:,3)-refp_ap(:,3)).^2);
r_fp = (2*sig.omega + sig.gamma*norm(refp_ap))/2;

d_app = sqrt(refp_ap(1)^2+refp_ap(2)^2+refp_ap(3)^2);
if d_app > 10
    inlier  = (d<= 0.45*sot-r_fp);
else
    inlier  = (d<= 0.35*sot-r_fp);
end

xyz = xyz(inlier,:); 
i = i(inlier);

%% Ausreißereliminierung mittels Ransac
[RSC] = fitPlaneRANSAC(xyz,0.002,2000,0.3*sot); % Ebene durch mittlere Punkte
n = RSC.ps(1:3);
d = RSC.ps(4);

% Verbesserung zu allen anderen Punkten
v = n(1)*xyz(:,1)+n(2)*xyz(:,2)+n(3)*xyz(:,3)-d;

% Grenzwert für Ausreißer (da kleine Winkel: sig.w*strecke zulässig)
T = 3*sqrt((sig.w*norm(refp_ap))^2 + (sig.w*norm(refp_ap))^2 + (sig.s_const+sig.s_prop*norm(refp_ap))^2);

% Indizes der Punkte, die weiter als der T von der Ebene entfernt liegen.
inlier = (abs(v) <= T);

xyz = xyz(inlier,:); % Punkte eliminieren
i = i(inlier,:);

%% ungefähre Einfallswinkel und Distanz
beta_ap = acos((refp_ap*n)./(norm(n)*norm(refp_ap)))/pi*180;
if beta_ap > 90
    beta_ap = beta_ap-180;
end
beta_ap = abs(beta_ap);
s_ap = norm(refp_ap);

return
