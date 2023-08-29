function TM = templatematchingalgorithm(I,sot,res,typ,d)
% author:   Jannik Janßen
% modified: 2018-01-09
% modified: 2018-01-15
% modified: 2018-02-03
% modified: 2018-02-15
% modified: 2018-03-01

%% Template erstellen
T = generateTemplate(sot,0.001,typ);

%% Templatematching - Schleife über alle Drehungen
modeToRot = 'bilinear';
bestoutput = zeros(1,5);
switch typ
    case 'cake'
        maxrot = 90;
    case 'classic'
        maxrot = 180;
end

if d > 10
    for k = 1:maxrot
        Trot = imrotate(T,k,modeToRot,'crop');
        output = dfttempmatch(I,Trot,100);

        corrVec(k,1) = output(1);
        alphaVec(k,1) = k;

        % remember output for smallest error
        if output(1) > bestoutput(1)
            bestoutput = [output,k];
        end
    end
else
    k = 0;
    for alpha = 1:0.25:maxrot
        k = k+1;
        Trot = imrotate(T,alpha,modeToRot,'crop');
        output = dfttempmatch(I,Trot,100);

        corrVec(k,1) = output(1);
        alphaVec(k,1) = alpha;

        % remember output for smallest error
        if output(1) > bestoutput(1)
            bestoutput = [output,k];
        end
    end 
end

%%
% % Plot für Paper
% figure
% plot(alphaVec,corrVec)
% xlabel('\alpha [°]')
% ylabel('correlation')
% xlim([0,180])
% ylim([0,1])
% set(gcf,'Paperposition',[0 0.1 6 4])
% set(gcf,'Papersize',[6 4]) % In Klammern Bounding Box in cm
% saveas(gcf, 'plots/s7_alpha','pdf') %'png'...
% 
% figure
% Trot = imrotate(T, bestoutput(5),modeToRot,'crop');
% [~,CC] = dfttempmatch(I,Trot,1);
% imshow(real(CC),[])
% set(gcf,'Paperposition',[0 0 4 4])
% set(gcf,'Papersize',[4 4]) % In Klammern Bounding Box in cm
% saveas(gcf, 'plots/s7_CC','pdf') %'png'...

%% Translationen (und beste Rotation)
corr = bestoutput(1);
dy = bestoutput(3);
dx = bestoutput(4);
alpha = bestoutput(5); % in Grad
T = imrotate(T,alpha,modeToRot,'crop');

%% Center in coordinatesystem of the image
tc_im = [dx,dy]+[0.5,0.5]+size(T)/2;

%% metdata
TM.T = T;
TM.alpha = alpha;
TM.dx = dx;
TM.dy = dy;
TM.tc_im = tc_im;
TM.corr = corr;
TM.corrVec = corrVec;
TM.alphaVec = alphaVec;
return