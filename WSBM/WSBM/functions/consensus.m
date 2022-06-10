function [] = consensus(commonPath, filename, al, k)
    %%% Consensus partition from WSBM

    rng('shuffle')
    warning('off', 'MATLAB:MKDIR:DirectoryExists');

    %%% Create folders
    path = sprintf("CSV/labels/%s", commonPath);
    mkdir(path, sprintf('al_%i', al*100)) 
    path_al = sprintf('%s/al_%i', path, al*100);
    mkdir(path_al, 'K')
    mkdir(sprintf("%s/K", path_al), 'CON')
    mkdir(sprintf("%s/K", path_al), 'MAX')
    mkdir(path_al, 'LOGEV')
    mkdir(sprintf('%s/LOGEV', path_al), sprintf('k_%i', k))
    mkdir(path_al, 'NMI')
    mkdir(sprintf('%s/NMI', path_al), sprintf('k_%i', k))

    folderVar=sprintf('Variables/%s/al_%i/K', commonPath, al*100);

    listfiles = {dir(folderVar).name};
    
    flag = false;
    for i=1:numel(listfiles)
        K = strsplit(listfiles{i}, 'k_');
        if numel(K) == 2
            K = strsplit(K{2}, '.');
            K = str2double(K{1});
            if  K == k
                flag = true;
                
                load(sprintf('%s/%s', folderVar, listfiles{i}), 'summary')

                %%% Get partition with maximum logEv
                [maxlogev, imax] = max(summary.scores);
                [~, labels_max] = max(summary.mu{imax},[],1);
                labels_max = labels_max';

                %%% Find the highest logEv 50-percentil partitions
                idx_selected_scores = find(summary.scores >= prctile(summary.scores, 50));

                %%% Find the VI between those partitions
                nsum = numel(idx_selected_scores);

                vi_k_al = zeros(nsum);

                for ii=1:nsum
                 for jj=1:nsum
                         if ii < jj
                             vi_k_al(ii,jj) = varInfo(summary.mu{idx_selected_scores(ii)}, ...
                                 summary.mu{idx_selected_scores(jj)});
                         end
                  end

                end

                vi_k_al = vi_k_al + vi_k_al';

                total_distances = sum(vi_k_al,2); 

                %%% From this set, the partition with minimum VI

                [~, min_var] = min(total_distances);

                min_var = idx_selected_scores(min_var);

                %%% Continue
                idxs = idx_selected_scores(total_distances <= prctile(total_distances, 50));

                nidxs = numel(idxs);

                %%% Partition alignment
	            disp("* Warning: Be careful selecting the right fln or commship file")
                csvPath = sprintf("CSV/%s/%s.csv", commonPath, filename);
                
                A = readtable(csvPath, 'HeaderLines', 0, 'Delimiter', ',');
                A = table2array(A);
                A(A == 0) = nan;

                [nx, ~] = size(A);
                M = zeros(nx, nidxs);

                for ii=1:nidxs
                    [~, labels] = max(summary.mu{idxs(ii)},[],1);
                    M(:,ii) = labels;
                end

                m = M(:, idxs == min_var);
                m1 = M(:,1);
                M(:,1) = m;
                M(:, idxs == min_var) = m1;

                for ii=2:nidxs
                    M(:,ii) = munkres_consensus(M(:, 1), M(:,ii));
            %             M(:,i) = jaccard_consensus(M(:, 1), M(:,i));
                end

                %%% Frequency prior

                F = zeros(nx,k);
                binc = 1:k;
                for ii=1:nx
                    counts = hist(M(ii,:),binc);
                    F(ii,:) = counts;
                end
                F = F/nidxs;
                F =  F';

                %%% Find new partitions using the WSBM and this frequency prior

                [labels, model] = wsbm(A, k, 'W_Distr', 'LogNormal', ...
                        'E_Distr', 'Normal', 'alpha', summary.alpha, 'seed', F);

                mu_max = zeros([k nx]);
                mu_cons = zeros([k nx]);

                for jj=1:nx
                   mu_max(labels_max(jj), jj) = 1;
                   mu_cons(labels(jj), jj) = 1;
                end

                NMI = nmi(mu_max, mu_cons);
                LogEvidence = model.Para.LogEvidence;

                summary = summary.models{imax};
                save(sprintf('Variables/%s/al_%i/MAX/k_%i.mat', commonPath, al*100, k), 'summary')
                summary = model;
                save(sprintf('Variables/%s/al_%i/CON/k_%i.mat', commonPath, al*100, k), 'summary')
            
                writematrix(LogEvidence, sprintf('%s/LOGEV/k_%i/logev_cons.csv', path_al, k))
                writematrix(maxlogev, sprintf('%s/LOGEV/k_%i/logev_max.csv', path_al, k))
                writematrix(NMI, sprintf('%s/NMI/k_%i/nmi_max_cons.csv', path_al, k))
                writematrix(labels, sprintf('%s/K/CON/k_%i.csv', path_al, k))
                writematrix(labels_max, sprintf('%s/K/MAX/k_%i.csv', path_al, k))

                disp(folderVar)
                fprintf('K: %i \nMaxLogEv: %f \nConLogEv: %f\n', K, maxlogev, LogEvidence)

            end

        end
    end
    
    if ~flag
        disp('No K found')
    end
    
end

% h = heatmap(NMI_matrix);
%  
% xlabels = ["mk1 theta"; "mk1 beta"; "mk1 highbeta"; "mk1 gamma"; ...
%     "mk2 theta"; "mk2 beta"; "mk2 highbeta"; "mk2 gamma"];
% ylabels = ["precue"; "postcue"];
%  
% h.XDisplayLabels = xlabels;
% h.YDisplayLabels = ylabels;
