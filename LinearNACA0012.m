function [ b1, b2, b3, b4, b5, b6, chord ] = LinearNACA0012()
    %% Import NACA 0012 profile
    %
    % The import of the airfoil data aids in the approximation of the structure
    % being analyzed. The data contains sufficient data points which allow for
    % accurate numerical integration and other discrete approximations.
    %
    Beta = atand(0.1/0.45);
    filename = 'user-000.csv';
    delimiter = ',';
    startRow = 10;
    endRow = 210;
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter',... 
        delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    xNACA0012 = dataArray{:, 1}/1.5939922481;
    yNACA0012 = dataArray{:, 2}/1.3429;
    clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
    filename = 'user-000.csv';
    delimiter = ',';
    startRow = 214;
    endRow = 314;
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter',...
        delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    xNACA0012_Camber = dataArray{:, 1}/1.5939922481;
    yNACA0012_Camber = dataArray{:, 2};
    clearvars filename delimiter startRow endRow formatSpec ...
        fileID dataArray ans;
    filename = 'user-000.csv';
    delimiter = ',';
    startRow = 318;
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
        'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);
    xNACA0012_Chord = dataArray{:, 1}/1.5939922481;
    yNACA0012_Chord = dataArray{:, 2};
    clearvars filename delimiter startRow formatSpec fileID dataArray ans;

    graphOut =0;

    if graphOut == 1;
        figure('Name','NACA0012')
        plot(xNACA0012,yNACA0012,xNACA0012_Camber,yNACA0012_Camber,...
            xNACA0012_Chord,yNACA0012_Chord)
        grid on
        grid minor
        axis([-1 101/1.5939922481 -25 25])
    end
    %% Top and Bottom NACA Linearization
    %
    % This process linearizes the NACA airfoil profile assuming the top and
    % bottom sections have no curvature.
    %
    chord = max(xNACA0012_Chord);
    initLinear = (1.249440000000000e+01)/1.5939922481; 
    endLinear  = (5.313950000000000e+01)/1.5939922481; 

    for ii = 1:numel(xNACA0012)
        if xNACA0012(ii) > initLinear & xNACA0012(ii) < endLinear
            if yNACA0012(ii) < 0
                replacement = yNACA0012(find(xNACA0012 == initLinear));
                yNACA0012(ii) = replacement(2);
            elseif yNACA0012(ii) > 0
                replacement = yNACA0012(find(xNACA0012 == initLinear));
                yNACA0012(ii) = replacement(1);
            end
        end
    end
    if graphOut == 1
        figure('Name','NACA0012 - Simplified')
        plot(xNACA0012,yNACA0012,xNACA0012_Camber,yNACA0012_Camber,...
            xNACA0012_Chord,yNACA0012_Chord)
        grid on
        grid minor
        axis([-1 101/1.5939922481 -25 25])
    end
    %% Calculate height
    %
    h = sum(abs(replacement));
    %% Calculate lengths

    b1 = 0.45*max(xNACA0012_Chord);
    b2 = b1;
    b4 = max(xNACA0012_Chord)*.25/cosd(Beta);
    b5 = b4;
    b6 = abs(h - 2*sind(Beta)*b4);
    %% Thickness Calculation

    t  = h/4;
    t1 = 0.4*t;
    t2 = t1;
    t3 = t2;
    t4 = t3;
    t5 = t4;
    t6 = 0.6*t;
    t7 = t6;
    %% Calculate length of the leading edge curvature

    cX = xNACA0012(find(xNACA0012<=(chord*.125)));
    cY = yNACA0012(find(xNACA0012<=(chord*.125)));
    b3 = 0;
    for ii = 2:numel(cX)
        b3 = b3 + sqrt(((cX(ii)-cX(ii-1))^2)+((cY(ii)-cY(ii-1))^2));
    end
end

