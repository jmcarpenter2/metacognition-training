function [ meanConf ] = FigureS3( analysis, results, plotFigs, exportFigs, figDir )
%FIGURES3 Runs group meanConf analysis and optionally plots and exports the
%Figure S3

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
    confRespPlot_twoSessions = figure;
    set(gcf, 'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    lineProps = struct('color', 'c', 'width', 2);
    mseb(1:numel([confResp.group_1.session_01.mean, confResp.group_1.session_02.mean]), [confResp.group_1.session_01.mean, confResp.group_1.session_02.mean], [confResp.group_1.session_01.sem, confResp.group_1.session_02.sem], lineProps); hold on;
    plot([108.5 108.5], [1 4], '-k', 'linewidth', 1);
    xlim([1 numel([confResp.group_1.session_01.mean, confResp.group_1.session_02.mean])]); ylim([1 4]);
    set(gca,'fontsize',14);
    t(1) = text(3,1.225, 'Session 1');
    t(2) = text(197,1.225, 'Session 2');
    set(t, 'fontsize', 14);
    title('Control Group');
    xlabel('Trials');
    ylabel('Confidence');
    subplot(1,2,2); % Experimental Group
    mseb(1:numel([confResp.group_2.session_01.mean, confResp.group_2.session_02.mean]), [confResp.group_2.session_01.mean, confResp.group_2.session_02.mean], [confResp.group_2.session_01.sem, confResp.group_2.session_02.sem], lineProps); hold on;
    plot([108.5 108.5], [1 4], '-k', 'linewidth', 1);
    xlim([1 numel([confResp.group_2.session_01.mean, confResp.group_2.session_02.mean])]); ylim([1 4]);
    set(gca,'fontsize',14);
    t(1) = text(3,1.225, 'Session 1');
    t(2) = text(197,1.225, 'Session 2');
    set(t, 'fontsize', 14);
    title('Experimental Group');
    xlabel('Trials');
    ylabel('Confidence');
    if exportFigs
        export_fig(fullfile(figDir, 'FigureS3.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', confRespPlot_twoSessions)
    end
end

