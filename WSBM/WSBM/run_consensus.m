function [] = run_consensus()
%run_consensus: Find the consensus partition for al alpha and R that the
%user has manually selected. Alpha and R have to be carefully put by the
%user.
%   No input variable

addpath('../analysis tools'); 
addpath('../'); 
addpath('functions');
addpath('analysis');

% Parameters
data = "commships";
mermoved = "merged";
distance = "tracto2016";
model = "normal_GB_GB";
filename = "3_commships";
commonPath = sprintf("%s/%s/%s/%s", data, mermoved, distance, model);
alpha=0;
R=2:15;

% Run code
for i=1:numel(R)
    for k=1:numel(alpha)
        mkdir(sprintf('Variables/%s/al_%i/', commonPath, alpha(k)*100), 'MAX');
        mkdir(sprintf('Variables/%s/al_%i/', commonPath, alpha(k)*100), 'CON');
        consensus(commonPath, filename, alpha(k), R(i))
    end
end

end
