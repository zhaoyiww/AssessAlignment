function [c,mu] = classifyIntensity(i);
% Klassifikation von Schwarz und Weiß mittels kmeans
%
% Input:
%       - i: Intensitätswerte [n,1]
% Output:
%       - c: Klasse der Intensitätswerte [n,1], 0 oder 1
%       - mu: Erwartungswert der Klassen (sortiert)
% 
% author: Jannik Janßen
try
    try
        [clas,mu] = kmeans(i,2);
        c = zeros(size(clas));
        if mu(1)>mu(2)
            c(clas ==1) = 1;
        else
            c(clas ==2) = 1;
        end
    catch
        disp('2. Clusteringversuch!')
        [clas,mu] = kmeans(i,2);
        c = zeros(size(clas));
        if mu(1)>mu(2)
            c(clas ==1) = 1;
        else
            c(clas ==2) = 1;
        end
        
    end
catch
    disp('3. Clusteringversuch!')
    [clas,mu] = kmeans(i,2);
    c = zeros(size(clas));
    if mu(1)>mu(2)
        c(clas ==1) = 1;
    else
        c(clas ==2) = 1;
    end
    
end

mu = sort(mu);

return