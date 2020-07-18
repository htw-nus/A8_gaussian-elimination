%% Gaussian Elimination Solver
% This script reads in a 3x3-Matrix and a 3x1-Vector
% Features:
    % Make the user repeat the input if the input
        % has the wrong size
        % is linear dependend
        % contains non-numeric values
    % transposes a 1x3-Vector into a 3x1-Vector
    % Valid Input sizes can be changed in the variables section as desired
% Output: clean 3x3-Matrix and a clean 3x1-Vector
% This Script solves the given system of linear equations if
% it has only one solution
% you can choose if you want to export the Input Data and the Solution

%% Clean Up

% clear command window
clc;
% clear variables
clearvars;
% close all windows
close all;

%% Welcome Text
% The text
welcome = sprintf(['This program solves a system of linear equations Ax = b\n'...
    'by using the gaussian elimination with row equilibration '...
    'and pivoting.\n\n']);
% Displaying the text
disp(welcome);


%% Frame conditions
% repeater to restart the solver if the solution is not valid
ScriptRepeater = 0;
while ~ScriptRepeater


%% Variables

% Define your Square Matrix Size here
m = 3;
% valid Matrix size
MatSize = [m,m];
% valid Vector sizes
VecSize = [m,1];
VecSizeAlt = flip(VecSize);
% repeater to re-prompt the user until a valid input is given
MatRepeater = 0;
% repeater to re-prompt the user until a valid input is given
VecRepeater = 0;


%% Error Messages

% not a Matrix of the desired size
err.WrongSizeMat = sprintf(['Error.\nInput must be a [' num2str(MatSize(1)) ']x'...
    '[' num2str(MatSize(2)) ']-Matrix.\n\n']);
% not a Vector of the desired size
err.WrongSizeVec = sprintf(['Error.\nInput must be a [' num2str(VecSize(1)) ']x'...
    '[' num2str(VecSize(2)) ']- or [' num2str(VecSizeAlt(1)) ']x'...
    '[' num2str(VecSizeAlt(2)) ']-Vector.\n\n']);
% non-numerical input
err.WrongType = sprintf(['Error.\nInput must be numeric only.\n\n']);
% linear dependent
err.LinDep = sprintf(['Error.\nThe given Matrix is linear dependent and\n'...
    'therefor has zero or an infinite amount of solutions.\n\n'...
    'Please enter a linear independent Matrix.\n\n']);
% linear dependent
err.InfSol = sprintf(['Error.\nThe given system of equations '...
    'has an infinite amount of solutions.\n\n'...
    'Please try another system of equations.\n\n']);


%% Input Matrix

% while 'repeater' is not TRUE
while ~MatRepeater
% ask the user to input a 3x3-Matrix with instructions
    MatUser = input(['Enter your [' num2str(MatSize(1)) ']x'...
        '[' num2str(MatSize(2)) ']-Matrix A (e.g. [1 2 3; 4 5 6; 7 8 9])'...
        '\n\n']);
    % check if input is not a 3x3-Matrix
    if ~isequal(size(MatUser), MatSize)
        % if TRUE make the user repeat the input
        MatRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongSizeMat);
    % Type checker
    elseif ~isnumeric(MatUser)
        % if Input is not numeric make the user repeat the imput
        MatRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongType);
    % check if the vectors are linear independent
    elseif det(MatUser) == 0
        % if input is linear dependent make the user repeat the input
        MatRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.LinDep);
    else
        % if FALSE exit while loop
        MatRepeater = 1;
    end
end

% clear command window
clc;

%% Input Vector
% while 'repeater' is not TRUE
while ~VecRepeater
% ask the user to input a 3x1-Vector with instructions
    VecUser = input(['Enter your [' num2str(VecSize(1)) ']x'...
    '[' num2str(VecSize(2)) ']- or [' num2str(VecSizeAlt(1)) ']x'...
    '[' num2str(VecSizeAlt(2)) ']-Vector b. (e.g. [9 8 7]): \n\n']);
    % check if input is not a 3x1- or 1x3-Vector
    if (~isequal(size(VecUser), VecSize) &...
            ~isequal(size(VecUser), VecSizeAlt))
        % if TRUE make the user repeat the input
        VecRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongSizeVec);
    elseif ~isnumeric(VecUser)
        % if Input is not numeric make the user repeat the imput
        VecRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongType);
    elseif size(VecUser) == [1,3]
        % if TRUE transpose the Vector
        VecUser = VecUser.';
        % and exit the while loop
        VecRepeater = 1;
    else
        % if FALSE exit the while loop
        VecRepeater = 1;
    end
end

% save input into a struct
UserInput.Mat = MatUser;
UserInput.Vec = VecUser;

% clear unnecessary vars
clearvars -except err ScriptRepeater UserInput m;
% clear command window
clc;


%% Row Scaling
% Calculation of the diagonal matrix
D = diag(1./sum(abs(UserInput.Mat),2));
% Scaling the Input with the diagonal matrix
ScalMat = D * UserInput.Mat;
ScalVec = D * UserInput.Vec;

%% Pivoting
% Determining the descending order of the absolute of the
% first element of each row
[~, order] = sort(abs(ScalMat(:,1)), 'descend');
% Ordering the scaled matrix and vector
OrdMat = ScalMat([order],:);
OrdVec = ScalVec([order],:);


%% Gaussian Elimination
% Source: https://www.youtube.com/watch?v=juDxSMpzv6c

% bringing the Matrix into upper triangular form
% initializing a new matrix to keep the previous results
GaussMat = OrdMat;
% and a new vector
GaussVec = OrdVec;
% for loop over m-1 elements
for i = 1:m-1
   % Calculating the Row Multiplier by dividing the i Element of the
   % row i+1:m through the i,i Element of the matrix
   Mult = GaussMat(i+1:m,i) / GaussMat(i,i);
   % Calculating the new rows of the matrix
   GaussMat(i+1:m,:) = GaussMat(i+1:m,:) - Mult*GaussMat(i,:);
   % and the vector
   GaussVec(i+1:m,:) = GaussVec(i+1:m,:) - Mult*GaussVec(i,:);
end

% removing rounding errors near zero
threshold=10^(-12);
% set near zeros to zero in the Matrix
GaussMat(abs(GaussMat)<threshold) = 0;
% set near zeros to zero in the Vector
GaussVec(abs(GaussVec)<threshold) = 0;
% initializing the solution vector x
X = zeros(m,1);
% determining the first element of the solution
X(m,:) = GaussVec(m,:)/GaussMat(m,m);
% for loop over the missing solutions
for i = m-1:-1:1
   % Calculating the i element of x by substracting product of the known
   % elements of x and the corresponding row of the Gauss Matrix
   % from the corresponding element of the Gauss Vector and dividing
   % the outcome by the matrix element i,i
   X(i,:) = (GaussVec(i,:) - GaussMat(i,i+1:m)*X(i+1:m,:))/GaussMat(i,i);
end

% clear unnecessary vars
clearvars -except err m ScriptRepeater UserInput X; 

% check if the gauss algorithm couldn't find a solution
% due to having an infinite amount of solutions
if sum(isnan(X))> 0
    % show error
    disp(err.InfSol)
    % and repeat the while-loop
    ScriptRepeater = 0;
else
    % Display the result
    disp('The result for x is:');
    disp(X);
    % else, leave the while-loop
    ScriptRepeater = 1;
end
end

%% Export Data
% Error Messages
% not a y/n Input
err.WrongInput = sprintf(['Error.\n I could not understand you.'...
    '\n\n']);
% repeater to re-prompt the user until a valid input is given
ExportRepeater = 0;
% while 'repeater' is not TRUE
while ~ExportRepeater
% ask the user if the data should be exported
    Export = input(['Do you want to export your Input and the Solution '...
        'as CSV? (y/n)\n\n'], 's');
    % save string in lowercase letters
    Export = lower(Export);
    % check if the user wants to export the data
    if Export == 'y'
        % if TRUE export the data
        ExportRepeater = 1;
        % clear command window
        clc;
        % Export the Input Matrix A
        writematrix(UserInput.Mat,'A.csv');
        % Export the Input Vector b
        writematrix(UserInput.Vec,'b.csv');
        % Export the Solution Vector x
        writematrix(X,'x.csv');
        % Inform the user about the export
        disp('The data has been exported.');
    % if he does not want to export the data
    elseif Export == 'n'
        % if TRUE export the data
        ExportRepeater = 1;
        % clear command window
        clc;
        % Inform the user about not exporting the data
        disp('The data has not been exported.');
    % else, the input is not valid
    else
        % if input is not valid make the user repeat the input
        ExportRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongInput);
    end
end

% clear unnecessary vars
clearvars -except err m UserInput X;


%% Plot
% Only show a plot if the third column contains no zeros
err.NoPlot = sprintf(['If you want to see a plot, please enter a Matrix\n'...
    'with only nonzero elements in the third column.\n\n']);
% Variable if User wants to see the plot
GraphIt = 0;
% Check if third column contains zeros
ZeroCheck = sum(UserInput.Mat(:,3) == 0);
% Only show a plot if m = 3
if m == 3
    % Only show a plot if the input can be solved for z
    if ZeroCheck == 0
        % repeater to re-prompt the user until a valid input is given
        plotRepeater = 0;
        % while 'repeater' is not TRUE
        while ~plotRepeater
        % ask the user if the data should be exported
            PlotIt = input(['\n\nDo you want to see a plot of your'...
                ' system of linear equations? (y/n)\n\n'], 's');
            % save string in lowercase letters
            PlotIt = lower(PlotIt);
            % check if the user wants to export the data
            if PlotIt == 'y'
                % if TRUE export the data
                plotRepeater = 1;
                % Variable to show the Plot
                GraphIt = 1;
                % clear command window
                clc;
                % Inform the user about the export
                disp('Your plot will open in an seperate window.');
            % if he does not want to export the data
            elseif PlotIt == 'n'
                % if TRUE export the data
                plotRepeater = 1;
                % clear command window
                clc;
                % Inform the user about not exporting the data
                disp('No plot will be shown.');
            % else, the input is not valid
            else
                % if input is not valid make the user repeat the input
                plotRepeater = 0;
                % clear command window
                clc;
                % Inform user about the error
                disp(err.WrongInput);
            end
        end
    else
        disp(err.NoPlot)
    end
else
    % else, don't show a plot
    GraphIt = 0;
end

% create the plot if the user wants it
if GraphIt == 1
    % Multiplication Factor as Diagonal Matrix
    DFactor = diag(1./UserInput.Mat(:,3));
    % Calculating the Equation matrix z = a + bx + cy
    Eqn = DFactor*[UserInput.Vec -UserInput.Mat(:,1:2)];
    % Calculating the plot area to make sure to show the intersaction X
    meshMax = max(max(2*ceil(abs(X(1:2)))),1);
    % set the mesh steps
    meshSteps = meshMax/15;
    % generate a meshgrid
    [x y] = meshgrid(-meshMax:meshSteps:meshMax);
    % reset the random seed to keep the colors stable
    rng default;
    % for loop 1:m
    for i = 1:m
        % Define equation z_i
        z = Eqn(i,1) + Eqn(i,2)*x + Eqn(i,3)*y;
        % generate the legend
        if Eqn(i,2) == 0
            % if the matrix entry is zero don't add text
            txtx2 = blanks(1);
        elseif Eqn(i,2) < 0
            % if it's below zero keep it as it is
            txtx2 = [blanks(1), num2str(round(Eqn(i,2),2)), 'x'];
        else
            % if it's greater than zero add a plus
            txtx2 = [' + ', num2str(round(Eqn(i,2),2)), 'x'];
        end
        if Eqn(i,3) == 0
            % if the matrix entry is zero don't add text
            txtx3 = blanks(1);
        elseif Eqn(i,3) < 0
            % if it's below zero keep it as it is
            txtx3 = [blanks(1), num2str(round(Eqn(i,3),2)), 'y'];
        else
            % if it's greater than zero add a plus
            txtx3 = [' +', num2str(round(Eqn(i,3),2)), 'y'];
        end
        % generate the full legend entry
        txt = ['z = ', num2str(round(Eqn(i,1),2)), txtx2, txtx3];
        % generate handles for the plot with legend text and a random color
        handles.HSurf(i) = [surf(x,y,z, 'DisplayName', txt,...
            'FaceColor', rand(1,3))];
        % set hold on to add multiple plots into one figure
        hold on;
    end
    % show the legend
    legend show;
    % set the x-Label to x
    xlabel('x');
    % set the y-Label to y
    ylabel('y');
    % set the z-Label to z and deactivate auto rotation
    zlabel('z', 'Rotation', 0);
    % set the size and position of the figure
    h=gcf;
    set(h, 'units','normalized','outerposition',[0.5 0.05 0.5 0.95]);
    %% Export the Plot
    % repeater to re-prompt the user until a valid input is given
    ExportRepeater = 0;
    % while 'repeater' is not TRUE
    while ~ExportRepeater
    % ask the user if the data should be exported
        Export = input(['\n\nDo you want to export the plot as SVG? '...
            '(y/n)\n\n'], 's');
        % save string in lowercase letters
        Export = lower(Export);
        % check if the user wants to export the plot
        if Export == 'y'
            % if TRUE export the plot
            ExportRepeater = 1;
            % clear command window
            clc;
            % rerun the plot in case the user has closed it and the handle
            % is lost
            % close all windows
            close all;
            % reset the random seed to keep the colors stable
            rng default;
            % for loop 1:m
            for i = 1:m
                % Define equation z_i
                z = Eqn(i,1) + Eqn(i,2)*x + Eqn(i,3)*y;
                % generate the legend
                if Eqn(i,2) == 0
                    % if the matrix entry is zero don't add text
                    txtx2 = blanks(1);
                elseif Eqn(i,2) < 0
                    % if it's below zero keep it as it is
                    txtx2 = [blanks(1), num2str(round(Eqn(i,2),2)), 'x'];
                else
                    % if it's greater than zero add a plus
                    txtx2 = [' + ', num2str(round(Eqn(i,2),2)), 'x'];
                end
                if Eqn(i,3) == 0
                    % if the matrix entry is zero don't add text
                    txtx3 = blanks(1);
                elseif Eqn(i,3) < 0
                    % if it's below zero keep it as it is
                    txtx3 = [blanks(1), num2str(round(Eqn(i,3),2)), 'y'];
                else
                    % if it's greater than zero add a plus
                    txtx3 = [' +', num2str(round(Eqn(i,3),2)), 'y'];
                end
                % generate the full legend entry
                txt = ['z = ', num2str(round(Eqn(i,1),2)), txtx2, txtx3];
                % generate handles for the plot with legend text and a random color
                handles.HSurf(i) = [surf(x,y,z, 'DisplayName', txt,...
                    'FaceColor', rand(1,3))];
                % set hold on to add multiple plots into one figure
                hold on;
            end
            % show the legend
            legend show;
            % set the x-Label to x
            xlabel('x');
            % set the y-Label to y
            ylabel('y');
            % set the z-Label to z and deactivate auto rotation
            zlabel('z', 'Rotation', 0);
            % set the size and position of the figure
            h=gcf;
            set(h, 'units','normalized','outerposition',[0.5 0.05 0.5 0.95])
            % Export the Plot as svg
            saveas(gcf,'SurfPlot.svg')
            % Inform the user about the export
            disp('The plot has been exported.');
        % if he does not want to export the data
        elseif Export == 'n'
            % if TRUE export the data
            ExportRepeater = 1;
            % clear command window
            clc;
            % Inform the user about not exporting the plot
            disp('The plot has not been exported.');
        % else, the input is not valid
        else
            % if input is not valid make the user repeat the input
            ExportRepeater = 0;
            % clear command window
            clc;
            % Inform user about the error
            disp(err.WrongInput);
        end
    end
end

clearvars -except err m UserInput X;