function xyz = rhv2xyz(rhv)
% Input:
%       - rhv: [nx3] Matrix with polar measurments
%
% Output:
%       - xyz: [nx3] Matrix with cartesian measurments
%
% Function calculates cartesian coordinates.
%
% author:   Jannik Janﬂen
% date:     2018-06-12

xyz = zeros(size(rhv));
xyz(:,1) = rhv(:,1).*sin(rhv(:,3)).*cos(rhv(:,2));
xyz(:,2) = rhv(:,1).*sin(rhv(:,3)).*sin(rhv(:,2));
xyz(:,3) = rhv(:,1).*cos(rhv(:,3));

end