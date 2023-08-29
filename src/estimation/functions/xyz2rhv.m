function rhv = xyz2rhv(xyz)
% Input:
%       - xyz: [nx3] Matrix with cartesian measurments
%
% Output:
%       - rhv: [nx3] Matrix with polar measurments
%
% Function calculates polar coordinates.
%
% author:   Jannik Janﬂen
% date:     2018-06-12

rhv = zeros(size(xyz));
rhv(:,1) = (xyz(:,1).*xyz(:,1)+xyz(:,2).*xyz(:,2)+xyz(:,3).*xyz(:,3)).^0.5;
rhv(:,2) = atan2(xyz(:,2),xyz(:,1));
rhv(:,3) = acos(xyz(:,3)./rhv(:,1));

end