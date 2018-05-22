function [ da ] = Figure2( results, plotFigs, exportFigs, figDir)
%FIGURE2 Runs group d' analysis and optionally plots and exports Figure 2

allSesh = 1:10;

dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
subjects = fieldnames(results);
% subjects = setdiff(subjects,'subject_098');% exclude this subject cuz poor M-ratio estimation

% Initialize arrays
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    da.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            da.(groups{g}).(session).perception.trained.raw = [];
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
                            da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(da.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.da);
                        end
                    end
                end
            else % Sessions 2-9
                da.(group).(sessions{sesh}).perception.trained.raw = vertcat(da.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.da);
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    [plots.da.(groups{g}).learningCurve.mean, plots.da.(groups{g}).learningCurve.sem] = deal([]);
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        if sesh == 1 || sesh == 10 || sesh == 100
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    da.(groups{g}).(session).(dom{d}).(stim{s}).mean = nanmean(da.(groups{g}).(session).(dom{d}).(stim{s}).raw);
                    da.(groups{g}).(session).(dom{d}).(stim{s}).sem = nanstd(da.(groups{g}).(session).(dom{d}).(stim{s}).raw)/sqrt(length(da.(groups{g}).(session).(dom{d}).(stim{s}).raw));
                end
            end
            plots.da.(groups{g}).(session).mean = [da.(groups{g}).(session).perception.trained.mean, da.(groups{g}).(session).perception.untrained.mean;...
                da.(groups{g}).(session).memory.trained.mean, da.(groups{g}).(session).memory.untrained.mean];
            plots.da.(groups{g}).(session).sem = [da.(groups{g}).(session).perception.trained.sem, da.(groups{g}).(session).perception.untrained.sem;...
                da.(groups{g}).(session).memory.trained.sem, da.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            da.(groups{g}).(session).perception.trained.mean = nanmean(da.(groups{g}).(session).perception.trained.raw);
            da.(groups{g}).(session).perception.trained.sem = nanstd(da.(groups{g}).(session).perception.trained.raw)/sqrt(length(da.(groups{g}).(session).perception.trained.raw));
            plots.da.(groups{g}).learningCurve.mean = vertcat(plots.da.(groups{g}).learningCurve.mean, da.(groups{g}).(session).perception.trained.mean);
            plots.da.(groups{g}).learningCurve.sem = vertcat(plots.da.(groups{g}).learningCurve.sem, da.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs   
    daCurvePlot = figure;
    set(gcf, 'position', [200 200 900 300]);    
    subplot(1,2,1); % Control Group
    hBar(1) = errorbar(-2.25, da.group_1.session_01.perception.trained.mean, da.group_1.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(-1.5, da.group_1.session_01.perception.untrained.mean, da.group_1.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(-0.75, da.group_1.session_01.memory.trained.mean, da.group_1.session_01.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(0, da.group_1.session_01.memory.untrained.mean, da.group_1.session_01.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    hCurve = errorbar(2:9, plots.da.group_1.learningCurve.mean, plots.da.group_1.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, da.group_1.session_10.perception.trained.mean, da.group_1.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(11.75, da.group_1.session_10.perception.untrained.mean, da.group_1.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(12.5, da.group_1.session_10.memory.trained.mean, da.group_1.session_10.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(13.25, da.group_1.session_10.memory.untrained.mean, da.group_1.session_10.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    plot([1,1],[-10,10],':k');
    plot([10,10],[-10,10],':k');
    plot([14.26,14.26],[-10,10],':k');
    ylim([1.1, 1.8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-1.125,2:9,12.125], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-3.25, 14.25]);
    xlabel('Session', 'fontsize', 14);
    ylabel('d''', 'fontsize', 14);
    title('Control Group', 'fontsize', 14);
    leg = legend('P, trained stimulus', 'P, untrained stimulus', 'M, trained stimulus', 'M, untrained stimulus', 'location', 'n');
    set(leg, 'FontSize', 8); 
    legend boxoff; box off;
    
    subplot(1,2,2); % Experimental Group
    hBar(1) = errorbar(-2.25, da.group_2.session_01.perception.trained.mean, da.group_2.session_01.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(-1.5, da.group_2.session_01.perception.untrained.mean, da.group_2.session_01.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(-0.75, da.group_2.session_01.memory.trained.mean, da.group_2.session_01.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(0, da.group_2.session_01.memory.untrained.mean, da.group_2.session_01.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    hCurve = errorbar(2:9, plots.da.group_2.learningCurve.mean, plots.da.group_2.learningCurve.sem, 'linewidth', 2.5);
    set(hCurve, 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    hBar(1) = errorbar(11, da.group_2.session_10.perception.trained.mean, da.group_2.session_10.perception.trained.sem,'o', 'linewidth', 2); hold on; % P, TS
    hBar(2) = errorbar(11.75, da.group_2.session_10.perception.untrained.mean, da.group_2.session_10.perception.untrained.sem,'o', 'linewidth', 2); hold on; % P, UTS
    hBar(3) = errorbar(12.5, da.group_2.session_10.memory.trained.mean, da.group_2.session_10.memory.trained.sem,'o', 'linewidth', 2); hold on; % M, TS
    hBar(4) = errorbar(13.25, da.group_2.session_10.memory.untrained.mean, da.group_2.session_10.memory.untrained.sem,'o', 'linewidth', 2); hold on; % M, UTS
    set(hBar(1), 'color', [.5 0 0], 'markerfacecolor', [.5 0 0], 'markersize', 10);
    set(hBar(2), 'color', [.5 0 0], 'markeredgecolor', [.5 0 0], 'markersize', 10);
    set(hBar(3), 'color', [0 0 .5], 'markerfacecolor', [0 0 .5], 'markersize', 10);
    set(hBar(4), 'color', [0 0 .5], 'markeredgecolor', [0 0 .5], 'markersize', 10);
    plot([1,1],[-10,10],':k');
    plot([10,10],[-10,10],':k');
    plot([14.26,14.26],[-10,10],':k');
    ylim([1.1, 1.8]);
    set(gca, 'fontsize', 14);
    set(gca, 'XTick', [-1.125,2:9,12.125], 'XTickLabel', {'Pre',2:9,'Post'}, 'fontsize', 11);
    xlim([-3.25, 14.25]);
    xlabel('Session', 'fontsize', 14);
    title('Experimental Group', 'fontsize', 14);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'Figure2.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', daCurvePlot)
    end
    
end
end