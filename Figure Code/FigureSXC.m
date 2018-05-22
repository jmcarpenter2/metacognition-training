function [ QSR ] = FigureSXC( analysis, results, plotFigs, exportFigs, figDir )
%FIGURESXC Compares QSR scores between real and simulated data and
%optionally plots and exports Figure SXC

allSesh = [1];

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
subjects = fieldnames(results.real);
% Initialize arrays
groups = {'real', 'sim'};
for g = 1:numel(groups)
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        TxTQSR.(groups{g}).(session).raw = [];
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            QSR.(groups{g}).(session).perception.abstract.raw = [];
        end
    end
end

% Concatenate raw data
for g = 1:numel(groups)
    for sub = 1:numel(subjects)
        sessions = fieldnames(results.real.(subjects{sub}));
        for sesh = 1:numel(sessions)
            if strncmp(sessions{sesh},'session',7)
                [tok, rem] = strtok(sessions{sesh}, '_');
                session = str2double(rem(2:end));
                if session == 1 || session == 10 || session == 100
                    for d = 1:numel(dom)
                        for s = 1:numel(stim)
                            if isfield(results.(groups{g}).(subjects{sub}).(sessions{sesh}).(dom{d}), stim{s})
                                TxTQSR.(groups{g}).(sessions{sesh}).raw = vertcat(TxTQSR.(groups{g}).(sessions{sesh}).raw, analysis.(groups{g}).(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).QSR');
                                QSR.(groups{g}).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(QSR.(groups{g}).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(groups{g}).(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanQSR);
                            end
                        end
                    end
                else % Sessions 2-9
                    TxTQSR.(groups{g}).(sessions{sesh}).raw = vertcat(TxTQSR.(groups{g}).(sessions{sesh}).raw, analysis.(groups{g}).(subjects{sub}).(sessions{sesh}).perception.abstract.QSR');
                    QSR.(groups{g}).(sessions{sesh}).perception.abstract.raw = vertcat(QSR.(groups{g}).(sessions{sesh}).perception.abstract.raw, results.(groups{g}).(subjects{sub}).(sessions{sesh}).perception.abstract.meanQSR);
                end
            end
        end    
    end
end

% Take difference + calculate mean and standard error
for d = 1:numel(dom)
   for s = 1:numel(stim)
       QSR.diff.session_01.(dom{d}).(stim{s}).raw = QSR.sim.session_01.(dom{d}).(stim{s}).raw - QSR.real.session_01.(dom{d}).(stim{s}).raw;
       QSR.diff.session_01.(dom{d}).(stim{s}).mean = nanmean(QSR.diff.session_01.(dom{d}).(stim{s}).raw);
       QSR.diff.session_01.(dom{d}).(stim{s}).sem = nanstd(QSR.diff.session_01.(dom{d}).(stim{s}).raw)/sqrt(length(QSR.diff.session_01.(dom{d}).(stim{s}).raw));
   end
end
plots.QSR.diff.session_01.mean = [QSR.diff.session_01.perception.abstract.mean, QSR.diff.session_01.perception.words.mean;...
    QSR.diff.session_01.memory.abstract.mean, QSR.diff.session_01.memory.words.mean];
plots.QSR.diff.session_01.sem = [QSR.diff.session_01.perception.abstract.sem, QSR.diff.session_01.perception.words.sem;...
    QSR.diff.session_01.memory.abstract.sem, QSR.diff.session_01.memory.words.sem];


% Take mean and standard error
for g = 1:numel(groups)
    [plots.QSR.(groups{g}).learningCurve.mean, plots.QSR.(groups{g}).learningCurve.sem] = deal([]);
    [plots.TxTQSR.(groups{g}).allSessions.mean, plots.TxTQSR.(groups{g}).allSessions.sem] = deal([]);
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        TxTQSR.(groups{g}).(session).mean = nanmean(TxTQSR.(groups{g}).(session).raw, 1);
        TxTQSR.(groups{g}).(session).sem = nanstd(TxTQSR.(groups{g}).(session).raw, [], 1)/sqrt(size(TxTQSR.(groups{g}).(session).raw,1));
        plots.TxTQSR.(groups{g}).allSessions.mean = horzcat(plots.TxTQSR.(groups{g}).allSessions.mean, TxTQSR.(groups{g}).(session).mean);
        plots.TxTQSR.(groups{g}).allSessions.sem = horzcat(plots.TxTQSR.(groups{g}).allSessions.sem, TxTQSR.(groups{g}).(session).sem);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    QSR.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(QSR.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.QSR.(groups{g}).(session).mean = [QSR.(groups{g}).(session).perception.abstract.mean, QSR.(groups{g}).(session).perception.words.mean;...
                QSR.(groups{g}).(session).memory.abstract.mean, QSR.(groups{g}).(session).memory.words.mean];
            plots.QSR.(groups{g}).(session).sem = [QSR.(groups{g}).(session).perception.abstract.sem, QSR.(groups{g}).(session).perception.words.sem;...
                QSR.(groups{g}).(session).memory.abstract.sem, QSR.(groups{g}).(session).memory.words.sem];
        else % sessions 2-9
            QSR.(groups{g}).(session).perception.abstract.mean = nanmean(QSR.(groups{g}).(session).perception.abstract.raw);
            QSR.(groups{g}).(session).perception.abstract.sem = nanstd(QSR.(groups{g}).(session).perception.abstract.raw)/sqrt(length(QSR.(groups{g}).(session).perception.abstract.raw));
            plots.QSR.(groups{g}).learningCurve.mean = vertcat(plots.QSR.(groups{g}).learningCurve.mean, QSR.(groups{g}).(session).perception.abstract.mean);
            plots.QSR.(groups{g}).learningCurve.sem = vertcat(plots.QSR.(groups{g}).learningCurve.sem, QSR.(groups{g}).(session).perception.abstract.sem);
        end
    end
end

if plotFigs
    % Error Bar Comparison Plot
    QSRComparisonPlot = figure;
    set(gcf,'position', [200 200 750 300]);
    subplot(1,3,1); % Real
    [hBar hErrorbar] = barwitherr([plots.QSR.real.session_01.sem(1,:), plots.QSR.real.session_01.sem(2,:)],...
        [plots.QSR.real.session_01.mean(1,:), plots.QSR.real.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([.65 .85]);
    set(gca, 'fontsize', 14);
    ylabel('QSR', 'fontsize', 14);
    title('Real Data', 'fontsize', 14);
    box off;
    
    subplot(1,3,2); % Conf Biased
    [hBar hErrorbar] = barwitherr([plots.QSR.sim.session_01.sem(1,:), plots.QSR.sim.session_01.sem(2,:)],...
        [plots.QSR.sim.session_01.mean(1,:), plots.QSR.sim.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([.65 .85]);
    set(gca, 'fontsize', 14);
    title('Confidence Biased Data', 'fontsize', 14);
    box off;
    
    subplot(1,3,3); % Difference
    [hBar hErrorbar] = barwitherr([plots.QSR.diff.session_01.sem(1,:), plots.QSR.diff.session_01.sem(2,:)],...
        [plots.QSR.diff.session_01.mean(1,:), plots.QSR.diff.session_01.mean(2,:)]);
    set(gca, 'xtick', [1:4],'xticklabel', {'P, abstract';'P, words';'M, abstract';'M, words'}, 'XTickLabelRotation', 45);
    ylim([-0.1 0.1]);
    set(gca, 'fontsize', 14);
    title('Difference', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'FigureSXC.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', QSRComparisonPlot)
    end
         
end
end