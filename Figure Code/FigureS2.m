function [ confDistr ] = FigureS2( analysis, plotFigs, exportFigs, figDir )
%FIGURES2 Runs group confidence distribution analysis and optionally
%plots and exports the Figure S2

evaluateSesh = [1,10];

% Initialize arrays
groups = {'group_1', 'group_2'};
dom = {'perception', 'memory'};
for g = 1:numel(groups)
    numSubjects.(groups{g}) = 0;
    for sesh = evaluateSesh
        session = sprintf('session_%.2d', sesh);
        confDistr.(groups{g}).(session).correct.raw = [];
        confDistr.(groups{g}).(session).incorrect.raw = [];
        for d = 1:numel(dom)
            confDistr.(groups{g}).(session).(dom{d}).correct.raw = [];
            confDistr.(groups{g}).(session).(dom{d}).incorrect.raw = [];
        end
        confDistr.(groups{g}).(session).perception.trained.correct.raw = [];
        confDistr.(groups{g}).(session).perception.trained.incorrect.raw = [];
    end
end

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
subjects = fieldnames(analysis);
% Concatenate raw data
for sub = 1:numel(subjects)
    if strncmp(subjects{sub},'subject',7)
        group = sprintf('group_%d', analysis.(subjects{sub}).group);
        numSubjects.(group) = numSubjects.(group) + 1;
        for sesh = evaluateSesh
            session = sprintf('session_%.2d', sesh);
            for d = 1:numel(dom)
                for s = 1:numel(stim)
                    if isfield(analysis.(subjects{sub}).(session).(dom{d}), stim{s})
                        confDistr.(group).(session).correct.raw = vertcat(confDistr.(group).(session).correct.raw, analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                        confDistr.(group).(session).incorrect.raw = vertcat(confDistr.(group).(session).incorrect.raw, (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                        confDistr.(group).(session).(dom{d}).correct.raw = vertcat(confDistr.(group).(session).(dom{d}).correct.raw, analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                        confDistr.(group).(session).(dom{d}).incorrect.raw = vertcat(confDistr.(group).(session).(dom{d}).incorrect.raw, (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                    end
                end
            end
            confDistr.(group).(session).perception.trained.correct.raw = vertcat(confDistr.(group).(session).perception.trained.correct.raw, analysis.(subjects{sub}).(session).perception.trained.accBinned);
            confDistr.(group).(session).perception.trained.incorrect.raw = vertcat(confDistr.(group).(session).perception.trained.incorrect.raw, (analysis.(subjects{sub}).(session).perception.trained.numTrialsBinned - analysis.(subjects{sub}).(session).perception.trained.accBinned));
        end  
    end
end

% Take mean and standard error
groups = {'group_1', 'group_2'};
for g = 1:numel(groups)
    for sesh = evaluateSesh
        session = sprintf('session_%.2d',sesh);
        confDistr.(groups{g}).(session).correct.mean = nanmean(confDistr.(groups{g}).(session).correct.raw,1);
        confDistr.(groups{g}).(session).correct.sem = nanstd(confDistr.(groups{g}).(session).correct.raw,0,1)/sqrt(numSubjects.(groups{g}));
        confDistr.(groups{g}).(session).incorrect.mean = nanmean(confDistr.(groups{g}).(session).incorrect.raw,1);    
        confDistr.(groups{g}).(session).incorrect.sem = nanstd(confDistr.(groups{g}).(session).incorrect.raw,0,1)/sqrt(numSubjects.(groups{g}));
        for d = 1:numel(dom)
            confDistr.(groups{g}).(session).(dom{d}).correct.mean = nanmean(confDistr.(groups{g}).(session).(dom{d}).correct.raw,1);
            confDistr.(groups{g}).(session).(dom{d}).correct.sem = nanstd(confDistr.(groups{g}).(session).(dom{d}).correct.raw,0,1)/sqrt(numSubjects.(groups{g}));
            confDistr.(groups{g}).(session).(dom{d}).incorrect.mean = nanmean(confDistr.(groups{g}).(session).(dom{d}).incorrect.raw,1);    
            confDistr.(groups{g}).(session).(dom{d}).incorrect.sem = nanstd(confDistr.(groups{g}).(session).(dom{d}).incorrect.raw,0,1)/sqrt(numSubjects.(groups{g}));
            
            plotConfDistr.(groups{g}).(session).(dom{d}) = [confDistr.(groups{g}).(session).(dom{d}).incorrect.mean/sum(confDistr.(groups{g}).(session).(dom{d}).incorrect.mean + confDistr.(groups{g}).(session).(dom{d}).correct.mean); confDistr.(groups{g}).(session).(dom{d}).correct.mean/sum(confDistr.(groups{g}).(session).(dom{d}).incorrect.mean + confDistr.(groups{g}).(session).(dom{d}).correct.mean)]';
        end
        confDistr.(groups{g}).(session).perception.trained.correct.mean = nanmean(confDistr.(groups{g}).(session).perception.trained.correct.raw,1);
        confDistr.(groups{g}).(session).perception.trained.correct.sem = nanstd(confDistr.(groups{g}).(session).perception.trained.correct.raw,0,1)/sqrt(numSubjects.(groups{g}));
        confDistr.(groups{g}).(session).perception.trained.incorrect.mean = nanmean(confDistr.(groups{g}).(session).perception.trained.incorrect.raw,1);
        confDistr.(groups{g}).(session).perception.trained.incorrect.sem = nanstd(confDistr.(groups{g}).(session).perception.trained.incorrect.raw,0,1)/sqrt(numSubjects.(groups{g}));
    end
end

if plotFigs
    confDistrPlot = figure;
    set(gcf,'position', [200 200 900 300]);
    subplot(1,2,1); % Control Group
    hBar = bar(1:4,[confDistr.group_1.session_01.incorrect.mean/sum(confDistr.group_1.session_01.incorrect.mean + confDistr.group_1.session_01.correct.mean); confDistr.group_1.session_01.correct.mean/sum(confDistr.group_1.session_01.incorrect.mean + confDistr.group_1.session_01.correct.mean)]'); hold on;
    hBar2 = bar(7:10,[confDistr.group_1.session_10.incorrect.mean/sum(confDistr.group_1.session_10.incorrect.mean + confDistr.group_1.session_10.correct.mean); confDistr.group_1.session_10.correct.mean/sum(confDistr.group_1.session_10.incorrect.mean + confDistr.group_1.session_10.correct.mean)]');
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,7:10],'xticklabel', [1:4,1:4]);
    xlim([0 11]);
    ylim([0 .5]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Control Group');
    t(1) = text(1.75,.3,'Pre');
    t(2) = text(7.75,.3,'Post');
    set(t, 'fontsize', 14);
    legend(hBar(2:-1:1), 'Correct','Incorrect','location','nw');
    legend boxoff;
    box off;
    
    subplot(1,2,2); % Experimental Group
    hBar = bar(1:4,[confDistr.group_2.session_01.incorrect.mean/sum(confDistr.group_2.session_01.incorrect.mean + confDistr.group_2.session_01.correct.mean); confDistr.group_2.session_01.correct.mean/sum(confDistr.group_2.session_01.incorrect.mean + confDistr.group_2.session_01.correct.mean)]'); hold on;
    hBar2 = bar(7:10,[confDistr.group_2.session_10.incorrect.mean/sum(confDistr.group_2.session_10.incorrect.mean + confDistr.group_2.session_10.correct.mean); confDistr.group_2.session_10.correct.mean/sum(confDistr.group_2.session_10.incorrect.mean + confDistr.group_2.session_10.correct.mean)]');
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(hBar2(1),'facecolor', 'r');
    set(hBar2(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4,7:10],'xticklabel', [1:4,1:4]);
    xlim([0 11]);
    ylim([0 .5]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Experimental Group');
    t(1) = text(1.75,.3,'Pre');
    t(2) = text(7.75,.3,'Post');
    set(t, 'fontsize', 14);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'FigureS2.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', confDistrPlot)
    end    

end