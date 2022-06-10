function [] = run_selection(t)
%run_selection: run the modelKSelection function for a certain iteration of
%alpha and R values. Alpha and R have to be specified in the script
%manually.
%   t is a scalar with the range 1 to the total number of alpha and R
%   values.

addpath('../analysis tools'); 
addpath('../'); 
addpath('functions');
addpath('analysis');

% Parameters
data = "commships";
mermoved = "merged";
distance = "tracto2016";
model = "normal_GB_GB";
filename = "4_commships";
commonPath = sprintf("%s/%s/%s/%s", data, mermoved, distance, model);

alpha = 0;
na = numel(alpha);
R = 6;
nr = numel(R);

alpha = repelem(alpha, nr);
R = repmat(R', na, 1)';

% Create important folders
mkdir('Variables', data);
mkdir(sprintf("Variables/%s", data), mermoved);
mkdir(sprintf('Variables/%s/%s', data, mermoved), distance);
mkdir(sprintf('Variables/%s/%s/%s', data, mermoved, distance), model);
mkdir(sprintf('Variables/%s', commonPath), sprintf('al_%i', alpha(t)*100));
mkdir(sprintf('Variables/%s/%s', commonPath, sprintf('al_%i', alpha(t)*100)), "K");

% Run code
modelKSelection(commonPath, filename, alpha(t), R(t))

end