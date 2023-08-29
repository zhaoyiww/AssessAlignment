%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------- %
%               Iteration loop for Janniks TCE algorithm                  %
% ----------------------------------------------------------------------- %
% Author: Tomislav Medic, 19.01.2018, Bonn                                %
% Modified: Zhaoyi Wang, 30.01.2023                                       %
% ----------------------------------------------------------------------- %
% Description:                                                            %
%                                                                         %
% INPUT                                                                   %
% ------------                                                           %
% Code Requires two main inputs,                                          % 
% a) .txt File containing approximate values of target coordinates        %
% b) .txt or .pts files containing corresponding target point clouds.     %
% Naming of targets in both files should match, and names should not have %
% more than 5 to 6 characters (letters and numbers).                      %
%                                                                         %
% OUTPUT                                                                  %
% ------------                                                            %
% Code outputs one big matrix/structure called "Results" which contains   %
% following data for each target:                                         %
% Target ID, approximate coordinates, estimated coordinates,              %
% distance to target, incidence angle, mean ratio of black & white        %
% intensity, image correlation value, mean white intensity, mean black    %
% intensity, standard deviation of BestFitPlanes both for black & white   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all;

tStart = tic;

addpath('functions')

%% Loading file of approximated values
curr_proj_path = '../../data';
radiometry_types = ["XYZI"];
% radiometry_types = ["XYZI", "XYZRGB"];

%% loading data from different folder in loop
for radiometry = 1:length(radiometry_types)
    tic

    abs_path = fullfile(curr_proj_path, 'results');
    apr = load(fullfile(abs_path, 'extraction', 'approximate_3Dcoordinates.txt'));      
    
    %% Creating Unique numeric identifier based on the point name
    % Note: Point ID can have maximum of 5-6 letters/symbols! 
    ID = apr(:, 1);
    n_targets = length(ID);
    apr_n = apr(:, 2:end);
    
    
    %% Defining folder containing point clouds:
    
    % Reading names of all txt files in the folder
    % Folder should contain only point clouds of interest!

    PathName = fullfile(abs_path, 'extraction', char(radiometry_types(radiometry)));
    A = dir(PathName);
    A = rmfield(A,{'date','bytes','isdir','datenum'});
    A = A(3:end);
    A = struct2cell(A)';
    m = length(A);     % number of files in the target folder
    
    
    %% Creating empty "Results" structure
    Results = struct();
    % ID of each target
    Results.ID = A;
    % Target centers:
    % app - approximated, est_XYZ - BestFitPlane Cartesian coordinates,
    % est_rhv - BestFitPlane Polar Coordinates
    Results.TCE.app = repmat(apr_n,m,1);
    Results.TCE.est = zeros(m*n_targets,3);
    % Incidence angles - estimated using different best fit planes
    Results.IA = zeros(m*n_targets,1);
    % TCE to TLS distance - estimated using different best fit planes
    Results.d = zeros(m*n_targets,1);
    % Image correlations
    Results.img_corr = zeros(m*n_targets,1);
    % Intensity values (median black, white, ratio)
    Results.int.avg_white = zeros(m*n_targets,1);
    Results.int.avg_black = zeros(m*n_targets,1);
    Results.int.bw_ratio = zeros(m*n_targets,1);

    % Target Normal
    Results.n = zeros(m*n_targets,3);
    
    %% Main Algorithm
    
    % Settings
    %sD.c   = 1;         % constant outlier threshold [mm]
    %sD.ppm = 80;        % proportional outlier threshold [ppm]
    sot    = 0.15;      % diameter of the target [m] / 15 for T'n'T, 20 for Paper
    typ = 'classic';

    % Scanner noise from manufacturers specifications
    % Defined for Leica ScanStation P20
    %
    %     % Range noise is defined as an exponential function of the distance
    %     % Exponential: stdev_d = a * e^(bx);
    %     % Not linear: a[m] + b[ppm]!
    %     % Standard deviation of ranges white
    %       sdev.wh_ra = 0.0003425;
    %       sdev.wh_rb = 0.01474;
    %     % Standard deviation of ranges black
    %       sdev.bl_ra = 0.0006963;
    %       sdev.bl_rb = 0.02562;
    %     % Standard deviation of angluar measurements
    %       sdev.h = deg2rad(3/3600); % [ " ] - arc seconds
    %       sdev.v = deg2rad(3/3600); % [ " ] - arc seconds
    
    % Janniks noise definition
    sig.w = deg2rad(1/3600);
    sig.s_const = 0.4/1000;
    sig.s_prop = 0.07/1000;
    sig.omega = 1.4/1000;
    sig.gamma = 0.2/1000;
     
     % clear and create the directory for saving estimated results
    if exist(fullfile(abs_path, 'estimation', ['images_' char(radiometry_types(radiometry))]), 'dir') ~= 0
        rmdir(fullfile(abs_path, 'estimation', ['images_' char(radiometry_types(radiometry))]), 's');
        delete(fullfile(abs_path, 'estimation', ['estimated_3Dcoordinates_' char(radiometry_types(radiometry)) '.txt']));
        delete(fullfile(abs_path, 'estimation', ['estimated_corr_values_' char(radiometry_types(radiometry)) '.txt']));
    end
    mkdir(fullfile(abs_path, 'estimation', ['images_' char(radiometry_types(radiometry))]));
    
    % loads each point cloud, preforms Janniks alg., stores results, goes again
    hv1 = '';
    hv2 = 1;
    for times = 1:m    
        % Loading correct point cloud
        hv1 = A{times};
        PC = dlmread(fullfile(PathName,hv1),' ', 1, 0);
%             PC = dlmread(fullfile(PathName,hv1),'\t', 1, 0);
        PC = single(PC);
        
        for targets = 1:n_targets
                     
            % Copying whole PC as a single precision values
            xyz = PC(:,1:3);
            if radiometry == 1
                i = PC(:,4);
            else
                % convert rgb values to color grayscale value
                i = 0.2989 * PC(:,4) + 0.5870 * PC(:,5) + 0.1140* PC(:,6);
            end
    
            % Loading approximate value
            xyzap = apr_n(targets,:);
    
            %% Algorithm starts
 
            % 1. Delete all points which are more than 0.3m away from
            % approximate center of the target
            % NOTE:
            % for very big point clouds reduce window to 0.15m or less
            d_point_vs_center = sqrt((xyz(:,1) - xyzap(1)).^2 + (xyz(:,2) - xyzap(2)).^2 + (xyz(:,3) - xyzap(3)).^2);
            xyz(d_point_vs_center > 0.15, :) = [];
            i(d_point_vs_center > 0.15) = [];
            clear d_point_vs_center;

            % Recalculating remaining points as a double precision values
            xyz = double(xyz);
            i = double(i);
 
            % SAFETY SWITCH AGAINST WRONG APPROXIMATE VALUES
            if length(i) < 20
                test = strcat({'Wrong approximate value for Nr. '},{' '},hv1);
                test = test{:};
%                     disp(test)
                hv2 = hv2 +1;
                continue
            end
            try
                % 4. Rough estimate of the target plane + outlier removal
                [xyz_pl,i_pl,n,d,beta_ap,s_ap] = pointsOfReferencePlane_div(xyz,i, xyzap, sot, sig);
                %idx_pl = ismember(xyz,xyz_pl,'rows');
            catch
                test = strcat({'Wrong approximate value for Nr. '},{' '},hv1);
                test = test{:};
%                     disp(test)
                hv2 = hv2 +1;
                continue
            end

            %% Step 2: Analyse Intensities
            [c,mu] = classifyIntensity(i_pl);
            mu(1) = median(i_pl(c == 0));
            mu(2) = median(i_pl(c == 1));
            ratio = (mu(2)-mu(1))/4096;

            % Storing relevant results
            Results.int.avg_white(hv2) = mu(2);
            Results.int.avg_black(hv2) = mu(1);
            Results.int.bw_ratio(hv2) = ratio;

            % Remove systematically influenced points
            if d < 10
                middle = median(i_pl(c==1));
                deviation = mad(i_pl(c==1),1);
                limit = 2.5*deviation;
                outliers = i_pl < middle-limit | i_pl > middle + limit;
                c(outliers == 1) = 0;
                clear i_pl;
            end

            try
                % 7. a) BEST FIT PLANES by Jannik (polar) + TCE:
                % --------------------------------------------
                % Schtzung der Ebene ber GHM nur mit den "weien" Punkten

                %% Step 3: Fit best Plane
                % Calculate polar coordinates
                rhv_pl = xyz2rhv(xyz_pl);
                % GHM with polar coordinates(Berit)
                GHM = fitPlaneGHMobs(rhv_pl(c==1,:), sig, [n;d]);
                n = GHM.n;
                d = GHM.d;
                clear GHM rhv_pl xyz_pl;


                %% Estimation of the Target Centre

                %% Step 4: Select Points for Intensity Image
                % Points for the intensity image (150% of the target size)
                dap = sqrt((xyz(:,1)-xyzap(:,1)).^2+(xyz(:,2)-xyzap(:,2)).^2+(xyz(:,3)-xyzap(:,3)).^2);
                inlier  = (dap<= 0.75*sot);
                xyz_im = xyz(inlier,:); i_im = i(inlier,:);
                %idx_im = ismember(xyz,xyz_im,'rows');
                clear dap inlier;

                % --------------------------------------------------------
                % NOTE: I can remove this step for targets in close range
                % if necessary!!!
                % (point cloud already cutted anyway)
                % --------------------------------------------------------

                % Scaling the intensity between -1 and 1 by using mu
                idxl = (i_im<=mu(1));
                idxm = (i_im > mu(1) & i_im < mu(2));
                idxu = (i_im>=mu(2));

                i_im(idxl) = -1;
                i_im(idxm) = ((i_im(idxm)-mu(1))./(mu(2)-mu(1)).*2)-1;
                i_im(idxu) = 1;
                clear idxl idxm idxu;

                %% Step 6: Derive intensity image
                % Project points into plane
                rhv_im = xyz2rhv(xyz_im);
                rhv_im = projectIntoPlane(rhv_im, sig,[n;d]);
                xyz_im = rhv2xyz(rhv_im);
                clear rhv_im;

                % Reduce to 2D
                [xy_im,trafopara] = reduceTo2D(xyz_im,n,d);

                clear xyz_im;

                % intensity image
                [I,X,Y] = interpIntIm(xy_im,i_im, 0.001);

                %% Step 7: Match Template
                % Template Matching
                TM = templatematchingalgorithm(I,sot,0.001,typ,d);
                corr = TM.corr;

                %% Step 8: Calculate the estimated centre
                % Target centre in 2D
                tc2D(1,1) = X(1,1)+(TM.tc_im(1)-1)/1000;
                tc2D(1,2) = Y(1,1)+(TM.tc_im(2)-1)/1000;

%                     clear X Y I;

                % Target centre in 3D
                tc3D = backTo3D(tc2D,trafopara);

                %% Incidence angle and distance
                % incidence angle
                beta = acos((tc3D*n)./(norm(n)*norm(tc3D)))/pi*180;
                if beta > 90
                    beta = beta-180;
                end
                beta = abs(beta);
                % distance
                s = norm(tc3D);

                % Storing results

                Results.TCE.est(hv2,:) = tc3D;
                Results.IA(hv2) = beta;
                Results.d(hv2) = s;
                Results.img_corr(hv2) = corr;
                Results.n(hv2,1:3) = n;

%                     figure()
%                     imshow(I)
%                     hold on
%                     scatter(TM.tc_im(1),TM.tc_im(2), 500, 'red', '+')
%                     saveas(gcf, [fullfile(abs_path, 'estimation', ['images_' char(radiometry_types(radiometry))], '/'), num2str(targets), '.jpg']);
%                     close all  
            catch
                test = strcat({'Some unexpected problem with point Nr. '},{' '},hv1);
                test = test{:};
                disp(test)
                hv2 = hv2 +1;
                continue
            end
    
            % Saying to me that everything is all right
            test = strcat({'Estimated '},{' '},{num2str(hv2)},{' out of '},{' '},{num2str(m*n_targets)},{' target centers'});
            test = test{:};
%                 disp(test)
    
            % Clearing some variables
            clear Polar_results_bl Polar_results_wh s beta refp3D refp2D trafopara xy c xyz;
            clear xy_im i;
            hv2 = hv2 +1;
        end
        % save detected images
        target_name = strsplit(hv1, '.');
        target_name = target_name(1);
        figure()
        imshow(I)
        hold on
        scatter(TM.tc_im(1),TM.tc_im(2), 500, 'red', '+')
        saveas(gcf, [fullfile(abs_path, 'estimation', ['images_' char(radiometry_types(radiometry))], '/'), char(target_name), '.jpg']);
        close all  
    end
    
    %% Creating new files
    
    % For Klein Altendorf
    
    target_number = length(apr_n);
    
    target_ID = repmat((1:target_number)',m,1);
    %jc = [target_ID Results.TCE.est];
    Results.ID = target_ID;
    
    
    % Saving files
    %ename = FileName(1:end-4);
    % ename2 = strcat(ename,'_jc.txt');
    % save(ename2,'jc','-ascii','-double');
    
    
    %% Cleaning of uncessary variables
    
    %  clear apr apr_n beta beta_ap bi_mean c connection d d_b_and_w hv1 i ID IDn m metadata mu n;
    %  clear n_b_and_w PathName PCIDn ratio refp2D refp3D s s_ap sD sdev sot std_r_black std_r_white;
    %  clear times trafopara typ wi_mean xy xyz xyzap stdev3D_black stdev3D_white;
    %  clear wi_median bi_median A Cartesian_b Cartesian_w Polar_results_bl Polar_results_wh d_Jannik n_Jannik d_approximate;
    clear PC;
    
%         save Results Results;
    
    a = round((Results.TCE.app - Results.TCE.est)*1000,1);
    c = sqrt((a(:,1).^2)+(a(:,2).^2)+(a(:,3).^2));
    
    % ename6 = strcat(ename,'_XYZ.txt');
    % writetable(cell2table([IDn num2cell(Results.TCE.est)]),ename6,'delimiter',',','WriteVariableNames',0);
    
    % imshow(I)
    % hold on
    % scatter(TM.tc_im(1),TM.tc_im(2), 500, 'red', '+')
    coord = Results.TCE.est;
%         coord(all(coord==0, 2), :) = [];
    coord = coord(1:m+1:end, :);
    [n_coord, ~] = size(coord);
    corr_values = Results.img_corr;
%         corr_values(all(corr_values==0, 2), :) = [];
    corr_values = corr_values(1:m+1:end, :);
   
    % save data to txt file
    % fid = fopen('.\results\coordinates.txt', 'w');
    fid = fopen(fullfile(abs_path, 'estimation', ['estimated_3Dcoordinates_' char(radiometry_types(radiometry)) '.txt']), 'w');
    fprintf(fid, 'x \t y \t z \r\n');
    for i = 1:n_coord
        fprintf(fid, '%.4f\t%.4f\t%.4f\r\n', coord(i, :));
    end
    fclose(fid);

    fid = fopen(fullfile(abs_path, 'estimation', ['estimated_corr_values_' char(radiometry_types(radiometry)) '.txt']), 'w');
    fprintf(fid, 'corr_values [0, 1] \r\n');
    for i = 1:n_coord
        fprintf(fid, '%.4f\r\n', corr_values(i, :));
    end
    fclose(fid);

    toc
    disp(['running time for current loop: ', num2str(toc), ' seconds']);
end 
close all

tEnd = toc(tStart);
disp(['total running time: ', num2str(tEnd), ' seconds']);