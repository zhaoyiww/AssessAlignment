function refp = backTo3D(refp,para)
%% Tranformiert einen 2D-Punkt zurück in 3D
% input:    refp [1x2] - 2D-Koordinate
%           para [struct] - Rotationsmatrix und Translationsvektor
% output:   xyz [1x3] - 3D-Koordinate
%
% author: Jannik Janßen
% date: 2018-02-19

refp = [refp,0]+para.t;
refp = inv(para.R)*refp';
refp = refp';

return
