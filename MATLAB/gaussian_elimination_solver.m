%% Read in User Input
% This script reads in a 3x3-Matrix and a 3x1-Vector
% Features:
    % Make the user repeat the input if the input
        % has the wrong size
        % is linear dependend
        % contains non-numeric values
    % transposes a 1x3-Vector into a 3x1-Vector
    % Valid Input sizes can be changed in the variables section as desired
% Output: clean 3x3-Matrix and a clean 3x1-Vector


%% Clean Up

% clear command window
clc;
% clear variables
clearvars;
% close all windows
close all;


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
err.WrongSizeMat = ['Error.\nInput must be a [' num2str(MatSize(1)) ']x'...
    '[' num2str(MatSize(2)) ']-Matrix.\n\n'];
% not a Vector of the desired size
err.WrongSizeVec = ['Error.\nInput must be a [' num2str(VecSize(1)) ']x'...
    '[' num2str(VecSize(2)) ']- or [' num2str(VecSizeAlt(1)) ']x'...
    '[' num2str(VecSizeAlt(2)) ']-Vector.\n\n'];
% non-numerical input
err.WrongType = sprintf(['Error.\nInput must be numeric only.\n\n']);
% linear dependent
err.LinDep= sprintf(['Error.\nThe given Matrix is linear dependent and\n'...
    'therefor has an infinite amount of solutions.\n\n'...
    'Please enter a linear independent Matrix\n\n']);


%% Input Matrix

% while 'repeater' is not TRUE
while ~MatRepeater
% ask the user to input a 3x3-Matrix with instructions
    MatUser = input(['Enter your [' num2str(MatSize(1)) ']x'...
        '[' num2str(MatSize(2)) ']-Matrix (e.g. [1 2 3; 4 5 6; 7 8 9])'...
        '\n\n']);
    % check if input is not a 3x3-Matrix
    if size(MatUser) ~= MatSize
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
    '[' num2str(VecSizeAlt(2)) ']-Vector. (e.g. [1 2 3]): \n\n']);
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
clearvars -except UserInput;
% clear command window
clc;

%% Variables Merge

% Combining the user input into an system of linear equations
MatInput = [UserInput.Mat UserInput.Vec];

%% Spaltenpivotisierung

%% Zeilenskalierung

%% Solver

%% Plot

%% Export
