function [] = modelKSelection(commonPath, filename, alpha, R)
%The function finds k communities in the selected 
%network. It saves the node-community assignation vector mu, logEvidence, 
%labels, and other parameters in a variable.

rng('shuffle')

summary = struct();

sprintf('al_%i_k_%i', alpha*100, R)

csvPath = sprintf("CSV/%s/%s.csv", commonPath, filename);

summary.alpha = alpha;
summary.W_Distr = 'LogNormal';
summary.E_Distr = 'Normal';
summary.iter = 100;
summary.R = R;
summary.parallel  = 0;
summary.numTrials = 100;
summary.mainMaxIter = 100;
summary.th = 0;
summary.name = sprintf('FLN_107_107_al_%i_th_%i_k_%i', summary.alpha*100, ...
    summary.th, summary.R);

A = readtable(csvPath, 'HeaderLines', 0, 'Delimiter', ',');
A = table2array(A);
A(A == 0) = nan;
A(A < summary.th/100) = nan;

summary.mu = cell(summary.iter,1);
summary.scores = zeros(summary.iter, 1);
summary.models = cell(summary.iter,1);

for i=1:summary.iter
    disp(i)
    [~, Model] = wsbm(A, summary.R, 'W_Distr', summary.W_Distr, ...
            'E_Distr',summary.E_Distr, 'numTrials', summary.numTrials, ...
            'alpha', summary.alpha, 'parallel', summary.parallel, ... 
            'mainMaxIter', summary.mainMaxIter, 'verbosity', 0);
    
    summary.scores(i) = Model.Para.LogEvidence;
    summary.mu{i,1} = Model.Para.mu;
    summary.models{i,1} = Model;
   
end

save(sprintf('Variables/%s/al_%i/K/k_%i.mat', commonPath, summary.alpha*100, summary.R), 'summary')

end
