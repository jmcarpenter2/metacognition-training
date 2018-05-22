function [ diff ] = FigureS1( results, plotFigs, exportFigs, figDir )
%FIGURES1 Runs group difficulty analysis and optionally plots and exports
%Figure S1

allSesh = 1:10;

dom = {'perception'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    diff.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            diff.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            [tok, rem] = strtok(sessions{sesh}, '_');
            session = str2double(rem(2:end));
            if session == 1 || session == 10 || session == 100
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        if isfield(results.(subjects{sub}).(sessions{sesh}).(dom{d}), stim{s})
                            diff.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(diff.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanDifficulty);
                        end
                    end
                end
            else % Sessions 2-9
                diff.(group).(sessions{sesh}).perception.trained.raw = vertcat(diff.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanDifficulty);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.diff.(groups{g}).learningCurve.mean, plots.diff.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    diff.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(diff.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    diff.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(diff.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(diff.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.diff.(groups{g}).(session).mean = [diff.(groups{g}).(session).perception.trained.mean, diff.(groups{g}).(session).perception.untrained.mean];
            plots.diff.(groups{g}).(session).sem = [diff.(groups{g}).(session).perception.trained.sem, diff.(groups{g}).(session).perception.untrained.sem];
        else % sessions 2-9
            diff.(groups{g}).(session).perception.trained.mean = nanmean(diff.(groups{g}).(session).perception.trained.raw);
            diff.(groups{g}).(session).perception.trained.sem = nanstd(diff.(groups{g}).(session).perception.trained.raw)/sqrt(length(diff.(groups{g}).(session).perception.trained.raw));
            plots.diff.(groups{g}).learningCurve.mean = vertcat(plots.diff.(groups{g}).learningCurve.mean, diff.(groups{g}).(session).perception.trained.mean);
            plots.diff.(groups{g}).learningCurve.sem = vertcat(plots.diff.(groups{g}).learningCurve.sem, diff.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    diffCurvePlot = figure;
    set(gcf,'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(-1, diff.group_1.session_01.perception.trained.mean, diff.group_1.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(0, diff.group_1.session_01.perception.untrained.mean, diff.group_1.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    hCurve = errorbar(2:9, plots.diff.group_1.learningCurve.mean, plots.diff.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, diff.group_1.session_10.perception.trained.mean, diff.group_1.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(12, diff.group_1.session_10.perception.untrained.mean, diff.group_1.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    plot([1,1],[100,200],':k');
    plot([10,10],[100,200],':k');
    plot([13.01,13.01],[100,200],':k');
    ylim([114, 124]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-0.5,2:9,11.5], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-2, 13]);
    xlabel('Session', 'fontsize', 14);
    ylabel('Difficulty Level', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained stimulus', 'P, untrained stimulus', 'location', 'n');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(-1, diff.group_2.session_01.perception.trained.mean, diff.group_2.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(0, diff.group_2.session_01.perception.untrained.mean, diff.group_2.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    hCurve = errorbar(2:9, plots.diff.group_2.learningCurve.mean, plots.diff.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, diff.group_2.session_10.perception.trained.mean, diff.group_2.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(12, diff.group_2.session_10.perception.untrained.mean, diff.group_2.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    plot([1,1],[100,200],':k');
    plot([10,10],[100,200],':k');
    plot([13.01,13.01],[100,200],':k');
    ylim([114, 124]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-0.5,2:9,11.5], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-2, 13]);
    xlabel('Session', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'FigureS1.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', diffCurvePlot)
    end
    
end
end