%% Read in User Input
% requirements:
    % size = 3x3
    % no variables

%% Clean Up

% clear command window
clc;
% clear variables
clearvars;
% close all windows
close all;


%% Error Messages

% not a 3x3-Matrix
err.WrongSize = 'Error. Input must be a 3x3-Matrix';

%% Code

% repeater to re-prompt the user until a valid input is given
repeater = 0;

% while 'repeater' is not TRUE
while ~repeater
% ask the user to input a 3x3-Matrix with instructions
    Mat_User = input('Enter your 3x3-Matrix (e.g. [1 2 3; 4 5 6; 7 8 9]');
    % check if input is not a 3x3-Matrix
    if size(Mat_User) ~= [3,3]
        % if TRUE make the user repeat the input
        repeater = 0;
        % clear command window
        clc;
        % Inform user about the error
        disp(err.WrongSize);
    else
        % if FALSE exit while loop
        repeater = 1;
    end
end