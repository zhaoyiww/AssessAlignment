function T= generateTemplate(sot,res,typ)
% Funktion generiert ein Template eines Zielzeichen mit normalen
% Schachbrettmuster.
%
% input:
%       - sot [scalar] Größe des realen Zielzeichens in [m]
%       - res [scalar] optional, Auflösung des Templates [m]
% output:
%       - T [matrix] Template des Zielzeichens
%
% date:     2017-11-28
% modified: 2018-01-15
% modified: 2018-02-01
% author:   Jannik Janßen

if nargin == 1
    nquad = floor((0.45*sot)*1000); % Anzahl der Pixel für das halbe Template
else
    nquad = floor((0.45*sot)/res); % Anzahl der Pixel für das halbe Template
end

if nargin <= 2
    % Template initialisieren
    T = [ones(nquad,nquad),-ones(nquad,nquad);
        -ones(nquad,nquad),ones(nquad,nquad)];
    
    for row = 1:size(T,1)
        for col = 1:size(T,2)
            dist = sqrt((row-0.5-size(T,1)/2)^2 + (col-0.5-size(T,2)/2)^2);
            if dist > nquad
                T(row,col) = 0;
            end
        end
    end
else
    switch typ
        case 'cake'
            T = triu(ones(nquad,nquad));
            T = [T,imrotate(T,-90); imrotate(T,90),imrotate(T,180)];
            
            for row = 1:size(T,1)
                for col = 1:size(T,2)
                    dist = sqrt((row-0.5-size(T,1)/2)^2 + (col-0.5-size(T,2)/2)^2);
                    if dist > nquad
                        T(row,col) = 0;
                    end
                end
            end
            %             T = imrotate(T,22.5, 'nearest','crop');
        case 'classic'
            % Template initialisieren
            T = [ones(nquad,nquad),-ones(nquad,nquad);
                -ones(nquad,nquad),ones(nquad,nquad)];
            
            for row = 1:size(T,1)
                for col = 1:size(T,2)
                    dist = sqrt((row-0.5-size(T,1)/2)^2 + (col-0.5-size(T,2)/2)^2);
                    if dist > nquad
                        T(row,col) = 0;
                    end
                end
            end
    end
end


end