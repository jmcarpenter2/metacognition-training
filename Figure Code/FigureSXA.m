function [ confDistr ] = FigureSXA( analysis, plotFigs, exportFigs, figDir )
%FIGURESXA Compares confidence distributions between real and simulated
%data and optionally plots and exports Figure SXA

evaluateSesh = [1];

% Initialize arrays
groups = {'real', 'sim'};
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
    end
end

dom = {'perception', 'memory'};
stim = {'abstract', 'words'};
subjects = fieldnames(analysis.real);
% Concatenate raw data
for g = 1:numel(groups)
    for sub = 1:numel(subjects)
        if strncmp(subjects{sub},'subject',7)
            numSubjects.(groups{g}) = numSubjects.(groups{g}) + 1;
            for sesh = evaluateSesh
                session = sprintf('session_%.2d', sesh);
                for d = 1:numel(dom)
                    for s = 1:numel(stim)
                        if isfield(analysis.(groups{g}).(subjects{sub}).(session).(dom{d}), stim{s})
                            confDistr.(groups{g}).(session).correct.raw = vertcat(confDistr.(groups{g}).(session).correct.raw, analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                            confDistr.(groups{g}).(session).incorrect.raw = vertcat(confDistr.(groups{g}).(session).incorrect.raw, (analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                            confDistr.(groups{g}).(session).(dom{d}).correct.raw = vertcat(confDistr.(groups{g}).(session).(dom{d}).correct.raw, analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                            confDistr.(groups{g}).(session).(dom{d}).incorrect.raw = vertcat(confDistr.(groups{g}).(session).(dom{d}).incorrect.raw, (analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned - analysis.(groups{g}).(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned)); % ./ (analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).accBinned + analysis.(subjects{sub}).(session).(dom{d}).(stim{s}).numTrialsBinned));
                        end
                    end
                end
            end  
        end
    end
end

% Take mean and standard error
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
    end
end

if plotFigs
    confDistrPlot = figure;
    set(gcf,'position', [200 200 500 300]);
    subplot(1,2,1); % Control Group
    hBar = bar(1:4,[confDistr.real.session_01.incorrect.mean/sum(confDistr.real.session_01.incorrect.mean + confDistr.real.session_01.correct.mean); confDistr.real.session_01.correct.mean/sum(confDistr.real.session_01.incorrect.mean + confDistr.real.session_01.correct.mean)]'); hold on;
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4],'xticklabel', [1:4]);
    xlim([0 5]);
    ylim([0 .5]);
    set(gca,'fontsize',14);
    ylabel('Proportion Responded');
    xlabel('Confidence Rating');
    title('Real Data');
    legend(hBar(2:-1:1), 'Correct','Incorrect','location','nw');
    legend boxoff;
    box off;
    
    subplot(1,2,2); % Experimental Group
    hBar = bar(1:4,[confDistr.sim.session_01.incorrect.mean/sum(confDistr.sim.session_01.incorrect.mean + confDistr.sim.session_01.correct.mean); confDistr.sim.session_01.correct.mean/sum(confDistr.sim.session_01.incorrect.mean + confDistr.sim.session_01.correct.mean)]'); hold on;
    set(hBar(1),'facecolor', 'r');
    set(hBar(2),'facecolor', 'g');
    set(gca, 'xtick', [1:4],'xticklabel', [1:4]);
    xlim([0 5]);
    ylim([0 .5]);
    set(gca,'fontsize',14);
    xlabel('Confidence Rating');
    title('Confidence Biased Data');
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'FigureSXA.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', confDistrPlot)
    end    

end