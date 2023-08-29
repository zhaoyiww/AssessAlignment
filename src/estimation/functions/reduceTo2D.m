function [xy,para] = reduceTo2D(xyz,n,d)
%% Reduziert die Koordinaten mit Hilfe der Ebenenparamter
% input:    xyz [nx3] - 3D-Koordinaten
%           n [3x1] - Normalvektor
%           d [scalar] - Abstand zum Ursprung
% output:   xy [nx2] - 2D-Koordinaten
%           para [struct] - Rotationsmatrix und Translationsvektor
% 
% author: Jannik Janﬂen
% date: 2018-02-19

n = n/norm(n);

e1 = cross(n,[0;0;1]); e1 = e1/norm(e1);
e2 = cross(n,e1); e2 = e2/norm(e2);
e3 = cross(e1,e2);

R = [e1';e2';e3'];
xyz = R*xyz';
t = [0,0,d];
xyz = xyz';


xy = xyz(:,1:2);
para.R = R;
para.t = t;

return