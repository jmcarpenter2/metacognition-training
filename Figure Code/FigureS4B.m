function [ M_ratio ] = FigureS4B( results, plotFigs, exportFigs, figDir )
%FIGURES4B Runs group M_ratio analysis and optionally plots and exports the
%Figure S4B

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
                    M_ratio.(groups{g}).(session).(dom{d}).(stim{s}).raw = [];
                end
            end
        else  % sessions 2-9
            M_ratio.(groups{g}).(session).perception.trained.raw = [];
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
                            M_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw = vertcat(M_ratio.(group).(sessions{sesh}).(dom{d}).(stim{s}).raw, results.(subjects{sub}).(sessions{sesh}).(dom{d}).(stim{s}).fit.M_ratio);
                        end
                    end
                end
            else % Sessions 2-9
%                 if ~(strcmp(subjects{sub},'subject_098') && session == 9)
                    M_ratio.(group).(sessions{sesh}).perception.trained.raw = vertcat(M_ratio.(group).(sessions{sesh}).perception.trained.raw, results.(subjects{sub}).(sessions{sesh}).perception.trained.fit.M_ratio);
%                 end
            end
        end
    end    
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
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
            plots.M_ratio.(groups{g}).(session).mean = [M_ratio.(groups{g}).(session).perception.trained.mean, M_ratio.(groups{g}).(session).perception.untrained.mean;...
                M_ratio.(groups{g}).(session).memory.trained.mean, M_ratio.(groups{g}).(session).memory.untrained.mean];
            plots.M_ratio.(groups{g}).(session).sem = [M_ratio.(groups{g}).(session).perception.trained.sem, M_ratio.(groups{g}).(session).perception.untrained.sem;...
                M_ratio.(groups{g}).(session).memory.trained.sem, M_ratio.(groups{g}).(session).memory.untrained.sem];
        else % sessions 2-9
            M_ratio.(groups{g}).(session).perception.trained.mean = nanmean(M_ratio.(groups{g}).(session).perception.trained.raw);
            M_ratio.(groups{g}).(session).perception.trained.sem = nanstd(M_ratio.(groups{g}).(session).perception.trained.raw)/sqrt(length(M_ratio.(groups{g}).(session).perception.trained.raw));
            plots.M_ratio.(groups{g}).learningCurve.mean = vertcat(plots.M_ratio.(groups{g}).learningCurve.mean, M_ratio.(groups{g}).(session).perception.trained.mean);
            plots.M_ratio.(groups{g}).learningCurve.sem = vertcat(plots.M_ratio.(groups{g}).learningCurve.sem, M_ratio.(groups{g}).(session).perception.trained.sem);
        end
    end
    
end

if plotFigs
    % Scatter Plot
    MRatioScatterPlot = figure;
    set(gcf, 'position', [200 200 650 450]);
    subplot(2,2,1); % Perception, Control Group
    hMarker(1) = scatter(M_ratio.group_1.session_01.perception.trained.raw, M_ratio.group_1.session_10.perception.trained.raw, 100, 'o', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_1.session_01.perception.untrained.raw, M_ratio.group_1.session_10.perception.untrained.raw, 100, 'o', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Perception meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,2); % Perception, Experimental Group
    hMarker(1) = scatter(M_ratio.group_2.session_01.perception.trained.raw, M_ratio.group_2.session_10.perception.trained.raw, 100, 'd', 'markerfacecolor', [.5 0 0], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_2.session_01.perception.untrained.raw, M_ratio.group_2.session_10.perception.untrained.raw, 100, 'd', 'markerfacecolor', [1 0 0], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Perception meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,3); % Memory, Control Group
    hMarker(1) = scatter(M_ratio.group_1.session_01.memory.trained.raw, M_ratio.group_1.session_10.memory.trained.raw, 100, 'o', 'markerfacecolor', [0 0 0.5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_1.session_01.memory.untrained.raw, M_ratio.group_1.session_10.memory.untrained.raw, 100, 'o', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('CG: Memory meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'nw');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    subplot(2,2,4); % Memory, Experimental Group
    hMarker(1) = scatter(M_ratio.group_2.session_01.memory.trained.raw, M_ratio.group_2.session_10.memory.trained.raw, 100, 'd', 'markerfacecolor', [0 0 .5], 'markeredgecolor', 'k'); hold on;
    hMarker(2) = scatter(M_ratio.group_2.session_01.memory.untrained.raw, M_ratio.group_2.session_10.memory.untrained.raw, 100, 'd', 'markerfacecolor', [0 0 1], 'markeredgecolor', 'k');
    plot([-5 5], [-5 5], '-k', 'linewidth', 2);
    xlim([0 4]); ylim([0 4]);
    set(gca, 'fontsize', 14);
    xlabel('Pre'); ylabel('Post');
    title('EG: Memory meta-d''/d''');
    leg = legend('Trained', 'Untrained', 'location', 'se');
    set(leg, 'fontsize', 6);
    legend boxoff; box off;
    if exportFigs
        export_fig(fullfile(figDir, 'FigureS4B.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', MRatioScatterPlot)
    end
    
end
end
    

    