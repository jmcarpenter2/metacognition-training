function [ M_ratio ] = FigureSXB( results, plotFigs, exportFigs, figDir )
%FIGURESXB Compares M_ratio scores between real and simulated
%data and optionally plots and exports Figure SXB

allSesh = [1];

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
subjects = fieldnames(results.real);

% Initialize arrays
groups = {'real', 'sim'};
for g = 1:numel(groups)
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            M_ratio.(groups{g}).(session).perception.abstract.raw = [];
        end
    end
end

% Concatenate raw data
for g = 1:numel(groups)
    for sub = 1:numel(subjects)
        sessions = fieldnames(results.(groups{g}).(subjects{sub}));
        for sesh = 1:numel(sessions)
            if strncmp(sessions{sesh},'session',7)
                [tok, rem] = strtok(sessions{sesh}, '_');
                session = str2double(rem(2:end));
                if session == 1 || session == 10 || session == 100
                    for d = 1:numel(dom)
                        for s = 1:numel(stim)
                            if isfield(results.(groups{g}).(subjects{sub}).(sessions{sesh}).(dom{d}), stim{s})
                                M_ratio.(groups{g}).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(M_ratio.(groups{g}).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(groups{g}).(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.M_ratio);
                            end
                        end
                    end
                else % Sessions 2-9
    %                 if ~(strcmp(subjects{sub},'subject_098') && session == 9)
                        M_ratio.(groups{g}).(sessions{sesh}).perception.abstract.raw = vertcat(M_ratio.(groups{g}).(sessions{sesh}).perception.abstract.raw, results.(groups{g}).(subjects{sub}).(sessions{sesh}).perception.abstract.fit.M_ratio);
    %                 end
                end
            end
        end    
    end
end

% Take difference + calculate mean and standard error
for d = 1:numel(dom)
   for s = 1:numel(stim)
       M_ratio.diff.session_01.(dom{d}).(stim{s}).raw = M_ratio.sim.session_01.(dom{d}).(stim{s}).raw - M_ratio.real.session_01.(dom{d}).(stim{s}).raw;
       M_ratio.diff.session_01.(dom{d}).(stim{s}).mean = nanmean(M_ratio.diff.session_01.(dom{d}).(stim{s}).raw);
       M_ratio.diff.session_01.(dom{d}).(stim{s}).sem = nanstd(M_ratio.diff.session_01.(dom{d}).(stim{s}).raw)/sqrt(length(M_ratio.diff.session_01.(dom{d}).(stim{s}).raw));
   end
end
plots.M_ratio.diff.session_01.mean = [M_ratio.diff.session_01.perception.abstract.mean, M_ratio.diff.session_01.perception.words.mean;...
    M_ratio.diff.session_01.memory.abstract.mean, M_ratio.diff.session_01.memory.words.mean];
plots.M_ratio.diff.session_01.sem = [M_ratio.diff.session_01.perception.abstract.sem, M_ratio.diff.session_01.perception.words.sem;...
    M_ratio.diff.session_01.memory.abstract.sem, M_ratio.diff.session_01.memory.words.sem];



% Take mean and standard error
for g = 1:numel(groups)
    [plots.M_ratio.(groups{g}).learningCurve.mean, plots.M_ratio.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.M_ratio.(groups{g}).(session).mean = [M_ratio.(groups{g}).(session).perception.abstract.mean, M_ratio.(groups{g}).(session).perception.words.mean;...
                M_ratio.(groups{g}).(session).memory.abstract.mean, M_ratio.(groups{g}).(session).memory.words.mean];
            plots.M_ratio.(groups{g}).(session).sem = [M_ratio.(groups{g}).(session).perception.abstract.sem, M_ratio.(groups{g}).(session).perception.words.sem;...
                M_ratio.(groups{g}).(session).memory.abstract.sem, M_ratio.(groups{g}).(session).memory.words.sem];
        else % sessions 2-9
            M_ratio.(groups{g}).(session).perception.abstract.mean = nanmean(M_ratio.(groups{g}).(session).perception.abstract.raw);
            M_ratio.(groups{g}).(session).perception.abstract.sem = nanstd(M_ratio.(groups{g}).(session).perception.abstract.raw)/sqrt(length(M_ratio.(groups{g}).(session).perception.abstract.raw));
            plots.M_ratio.(groups{g}).learningCurve.mean = vertcat(plots.M_ratio.(groups{g}).learningCurve.mean, M_ratio.(groups{g}).(session).perception.abstract.mean);
            plots.M_ratio.(groups{g}).learningCurve.sem = vertcat(plots.M_ratio.(groups{g}).learningCurve.sem, M_ratio.(groups{g}).(session).perception.abstract.sem);
        end
    end
    
end

if plotFigs
    % Error Bar Comparison Plot
    M_ratioComparisonPlot = figure;
    set(gcf,'position', [200 200 750 300]);
    subplot(1,3,1); % Real
    [hBar hErrorbar] = barwitherr([plots.M_ratio.real.session_01.sem(1,:), plots.M_ratio.real.session_01.sem(2,:)],...
        [plots.M_ratio.real.session_01.mean(1,:), plots.M_ratio.real.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([0 2]);
    set(gca, 'fontsize', 14);
    ylabel('meta-d''/d''', 'fontsize', 14);
    title('Real Data', 'fontsize', 14);
    box off;
    
    subplot(1,3,2); % Conf Biased
    [hBar hErrorbar] = barwitherr([plots.M_ratio.sim.session_01.sem(1,:), plots.M_ratio.sim.session_01.sem(2,:)],...
        [plots.M_ratio.sim.session_01.mean(1,:), plots.M_ratio.sim.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([0 2]);
    set(gca, 'fontsize', 14);
    title('Confidence Biased Data', 'fontsize', 14);
    box off;
    
    subplot(1,3,3); % Difference
    [hBar hErrorbar] = barwitherr([plots.M_ratio.diff.session_01.sem(1,:), plots.M_ratio.diff.session_01.sem(2,:)],...
        [plots.M_ratio.diff.session_01.mean(1,:), plots.M_ratio.diff.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([-0.3 0.3]);
    set(gca, 'fontsize', 14);
    title('Difference', 'fontsize', 14);
    box off;
    
    if exportFigs
        export_fig(fullfile(figDir, 'FigureSXB.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', M_ratioComparisonPlot)
    end
    
end
end
    

    