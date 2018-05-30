function [ meanConf ] = Figure4A( analysis, results, plotFigs, exportFigs, figDir )
%FIGURE4A Runs group meanConf analysis and optionally plots and exports
%Figure 4A

allSesh = 1:10;

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);

% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        confResp.(groups{g}).(session).raw = [];
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            meanConf.(groups{g}).(session).perception.trained.raw = [];
        end
    end
end

% Concatenate raw data
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    sessions = fieldnames(results.(subjects{sub}));
    meanConf.(group).(subjects{sub}).correct.raw = [];
    meanConf.(group).(subjects{sub}).incorrect.raw = [];
    for sesh = 1:numel(sessions)
        if strncmp(sessions{sesh},'session',7)
            [tok, rem] = strtok(sessions{sesh}, '_');
            session = str2double(rem(2:end));
            if session == 1 || session == 10 || session == 100
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        if isfield(results.(subjects{sub}).(sessions{sesh}).(dom{d}), stim{s})
                            confResp.(group).(sessions{sesh}).raw = vertcat(confResp.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp');
                            meanConf.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(meanConf.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).meanConf);
                            meanConf.(group).(subjects{sub}).correct.raw = vertcat(meanConf.(group).(subjects{sub}).correct.raw, nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc == 1)));
                            meanConf.(group).(subjects{sub}).incorrect.raw = vertcat(meanConf.(group).(subjects{sub}).incorrect.raw,  nanmean(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).confResp(analysis.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).acc == 0)));
                        end
                    end
                end
            else % Sessions 2-9
                confResp.(group).(sessions{sesh}).raw = vertcat(confResp.(group).(sessions{sesh}).raw, analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp');
                meanConf.(group).(sessions{sesh}).perception.trained.raw = vertcat(meanConf.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.meanConf);
                meanConf.(group).(subjects{sub}).correct.raw = vertcat(meanConf.(group).(subjects{sub}).correct.raw, nanmean(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.acc == 1)));
                meanConf.(group).(subjects{sub}).incorrect.raw = vertcat(meanConf.(group).(subjects{sub}).incorrect.raw,  nanmean(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.confResp(analysis.(subjects{sub}).(sessions{sesh}).perception.trained.acc == 0)));
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.meanConf.(groups{g}).learningCurve.mean, plots.meanConf.(groups{g}).learningCurve.sem] = deal([]);
    [plots.confResp.(groups{g}).allSessions.mean, plots.confResp.(groups{g}).allSessions.sem] = deal([]);
    [meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).incorrect.raw] = deal([]);
    subjects = fieldnames(meanConf.(groups{g}));
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            meanConf.(groups{g}).(subjects{sub}).correct.mean = nanmean(meanConf.(groups{g}).(subjects{sub}).correct.raw);
            meanConf.(groups{g}).(subjects{sub}).incorrect.mean = nanmean(meanConf.(groups{g}).(subjects{sub}).incorrect.raw);
            meanConf.(groups{g}).correct.raw = vertcat(meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).(subjects{sub}).correct.mean);
            meanConf.(groups{g}).incorrect.raw = vertcat(meanConf.(groups{g}).incorrect.raw, meanConf.(groups{g}).(subjects{sub}).incorrect.mean);
        end
    end
    meanConf.(groups{g}).correct.mean = nanmean(meanConf.(groups{g}).correct.raw);
    meanConf.(groups{g}).correct.sem = nanstd(meanConf.(groups{g}).correct.raw)/sqrt(length(meanConf.(groups{g}).correct.raw));
    meanConf.(groups{g}).incorrect.mean = nanmean(meanConf.(groups{g}).incorrect.raw);
    meanConf.(groups{g}).incorrect.sem = nanstd(meanConf.(groups{g}).incorrect.raw)/sqrt(length(meanConf.(groups{g}).incorrect.raw));
    [meanConf.(groups{g}).h, meanConf.(groups{g}).p, meanConf.(groups{g}).ci, meanConf.(groups{g}).stats] = ttest(meanConf.(groups{g}).correct.raw, meanConf.(groups{g}).incorrect.raw);
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        confResp.(groups{g}).(session).mean = nanmean(confResp.(groups{g}).(session).raw,1);
        confResp.(groups{g}).(session).sem = nanstd(confResp.(groups{g}).(session).raw,[],1)/sqrt(size(confResp.(groups{g}).(session).raw,1));
        plots.confResp.(groups{g}).allSessions.mean = horzcat(plots.confResp.(groups{g}).allSessions.mean, confResp.(groups{g}).(session).mean);
        plots.confResp.(groups{g}).allSessions.sem = horzcat(plots.confResp.(groups{g}).allSessions.sem, confResp.(groups{g}).(session).sem);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                    meanConf.(groups{g}).(session).(dom{d}).(stim{s}).std = nanstd(meanConf.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                end
            end
            plots.meanConf.(groups{g}).(session).mean = [meanConf.(groups{g}).(session).perception.trained.mean, meanConf.(groups{g}).(session).perception.untrained.mean;...
                meanConf.(groups{g}).(session).memory.trained.mean, meanConf.(groups{g}).(session).memory.untrained.mean];
            plots.meanConf.(groups{g}).(session).sem = [meanConf.(groups{g}).(session).perception.trained.sem, meanConf.(groups{g}).(session).perception.untrained.sem;...
                meanConf.(groups{g}).(session).memory.trained.sem, meanConf.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            meanConf.(groups{g}).(session).perception.trained.mean = nanmean(meanConf.(groups{g}).(session).perception.trained.raw);
            meanConf.(groups{g}).(session).perception.trained.sem = nanstd(meanConf.(groups{g}).(session).perception.trained.raw)/sqrt(length(meanConf.(groups{g}).(session).perception.trained.raw));
            plots.meanConf.(groups{g}).learningCurve.mean = vertcat(plots.meanConf.(groups{g}).learningCurve.mean, meanConf.(groups{g}).(session).perception.trained.mean);
            plots.meanConf.(groups{g}).learningCurve.sem = vertcat(plots.meanConf.(groups{g}).learningCurve.sem, meanConf.(groups{g}).(session).perception.trained.sem);
        end
    end
end

if plotFigs  
    meanConfCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(-2.25, meanConf.group_1.session_01.perception.trained.mean, meanConf.group_1.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(-1.5, meanConf.group_1.session_01.perception.untrained.mean, meanConf.group_1.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(-0.75, meanConf.group_1.session_01.memory.trained.mean, meanConf.group_1.session_01.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(0, meanConf.group_1.session_01.memory.untrained.mean, meanConf.group_1.session_01.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meanConf.group_1.learningCurve.mean, plots.meanConf.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, meanConf.group_1.session_10.perception.trained.mean, meanConf.group_1.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(11.75, meanConf.group_1.session_10.perception.untrained.mean, meanConf.group_1.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(12.5, meanConf.group_1.session_10.memory.trained.mean, meanConf.group_1.session_10.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(13.25, meanConf.group_1.session_10.memory.untrained.mean, meanConf.group_1.session_10.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    plot([1,1],[-10,10],':k');
    plot([10,10],[-10,10],':k');
    plot([14.26,14.26],[-10,10],':k');
    ylim([2 4]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-1.125,2:9,12.125], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-3.25, 14.25]);
    xlabel('Session', 'fontsize', 14);
    ylabel({'Metacognitive bias','(confidence level)'}, 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained', 'P, untrained stimuus', 'M, trained stimulus', 'M, untrained stimulus', 'location', 'n');
    set(leg, 'FontSize', 8); 
    legend boxoff; 
    box off;
    
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(-2.25, meanConf.group_2.session_01.perception.trained.mean, meanConf.group_2.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(-1.5, meanConf.group_2.session_01.perception.untrained.mean, meanConf.group_2.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(-0.75, meanConf.group_2.session_01.memory.trained.mean, meanConf.group_2.session_01.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(0, meanConf.group_2.session_01.memory.untrained.mean, meanConf.group_2.session_01.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    hCurve = errorbar(2:9, plots.meanConf.group_2.learningCurve.mean, plots.meanConf.group_2.learningCurve.sem,'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, meanConf.group_2.session_10.perception.trained.mean, meanConf.group_2.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(11.75, meanConf.group_2.session_10.perception.untrained.mean, meanConf.group_2.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(12.5, meanConf.group_2.session_10.memory.trained.mean, meanConf.group_2.session_10.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(13.25, meanConf.group_2.session_10.memory.untrained.mean, meanConf.group_2.session_10.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    plot([1,1],[-10,10],':k');
    plot([10,10],[-10,10],':k');
    plot([14.26,14.26],[-10,10],':k');
    ylim([2 4]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-1.125,2:9,12.125], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-3.25, 14.25]);
    xlabel('Session', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'Figure4A.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', meanConfCurvePlot)
    end
    
end
end