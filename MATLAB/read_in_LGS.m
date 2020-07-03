%% Read in User Input
% This script reads in a 3x3-Matrix and a 3x1-Vector
% Features:
    % Checks the Input Size and repeats the request if the size is not
    % valid
    % transposes a 1x3-Vector into a 3x1-Vector
    % Valid Input sizes can be changed in the variables section as desired
% Output: clean 3x3-Matrix and a clean 3x1-Vector


%% Clean Up

% clear command window
clc;
% clear variables
clear variables;
% close all windows
close all;


%% Variables

% repeater to re-prompt the user until a valid input is given
MatRepeater = 0;
% repeater to re-prompt the user until a valid input is given
VecRepeater = 0;
% valid Matrix size
MatSize = [3,3];
% valid Vector sizes
VecSize = [3,1];
VecSizeAlt = [1,3];


%% Error Messages

% not a Matrix of the desired size
err.WrongSizeMat = ['Error. Input must be a [' num2str(MatSize(1)) ']x'...
    '[' num2str(MatSize(2)) ']-Matrix.'];
% not a Vector of the desired size
err.WrongSizeVec = ['Error. Input must be a [' num2str(VecSize(1)) ']x'...
    '[' num2str(VecSize(2)) ']- or [' num2str(VecSizeAlt(1)) ']x'...
    '[' num2str(VecSizeAlt(2)) ']-Vector.'];


%% Input Matrix

% while 'repeater' is not TRUE
while ~MatRepeater
% ask the user to input a 3x3-Matrix with instructions
    Mat_User = input(['Enter your [' num2str(MatSize(1)) ']x'...
        '[' num2str(MatSize(2)) ']-Matrix (e.g. [1 2 3; 4 5 6; 7 8 9]) \n\n']);
    % check if input is not a 3x3-Matrix
    if size(Mat_User) ~= MatSize
        % if TRUE make the user repeat the input
        MatRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongSizeMat);
    else
        % if FALSE exit while loop
        MatRepeater = 1;
    end
end


%% Input Vector
% while 'repeater' is not TRUE
while ~VecRepeater
% ask the user to input a 3x1-Vector with instructions
    Vec_User = input(['\nEnter your [' num2str(VecSize(1)) ']x'...
    '[' num2str(VecSize(2)) ']- or [' num2str(VecSizeAlt(1)) ']x'...
    '[' num2str(VecSizeAlt(2)) ']-Vector. (e.g. [1 2 3]): \n\n']);
    % check if input is not a 3x1- or 1x3-Vector
    if (~isequal(size(Vec_User), VecSize) & ~isequal(size(Vec_User), VecSizeAlt))
        % if TRUE make the user repeat the input
        VecRepeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongSizeVec);
    elseif size(Vec_User) == [1,3]
        % if TRUE transpose the Vector
        Vec_User = Vec_User.';
        % and exit the while loop
        VecRepeater = 1;
    else
        % if FALSE exit the while loop
        VecRepeater = 1;
    end
end

clear err MatRepeater VecRepeater VecSizeAlt