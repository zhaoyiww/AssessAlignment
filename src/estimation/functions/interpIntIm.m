function [I,X,Y] = interpIntIm(xy,i,res)
% This function interpolate an intensity image from the given coordinates
% xy and intensity 
%
% modified: 2018-02-03
%
% author: Jannik Janßen

%% Grid erzeugen
[X,Y] = meshgrid([min(xy(:,1)):res:max(xy(:,1))],...
    [min(xy(:,2)):res:max(xy(:,2))]);

%% Wenn nur eine Farbe gescannt wurde
if length(unique(i)) == 1
    xy_new = [X(:),Y(:)];
    inb = zeros(size(xy_new,1),1);
    for k = 1:1:length(inb)
        
        d = sqrt((xy_new(k,1)-xy(:,1)).^2+(xy_new(k,2)-xy(:,2)).^2);
        if min(d) > sqrt(res.^2+res.^2)
            inb(k) = 1;
        end
    end
    inb = logical(inb);
    xy_new = xy_new(inb,:);
    
    xy = [xy;xy_new];
    i = [i;zeros(size(xy_new,1),1)];
end

%% Dateninterpolation/approximation
I = griddata(xy(:,1),xy(:,2),i,X,Y,'natural');

I(isnan(I)) = 0;

end