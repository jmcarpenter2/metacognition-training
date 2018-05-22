function [ ] = Figure6AB( results, plotFigs, exportFigs, figDir )
%FIGURE6AB Performs analysis on correlatation between mean confidence
%(bias) and log(meta-d'/d')(efficiency) for only perception/trained
% Optionally plots and exports Figures 6A, 6B

allSesh = 1:10;

groups = {'group_1', 'group_2'};
dom = {'perception', 'memory'};
stim = {'trained', 'untrained'};
for g = 1:numel(groups)
    numSubjects.(groups{g}) = 0;
    M_ratioDerivativePeak.(groups{g}).idx.raw = [];
    metaDaDerivativePeak.(groups{g}).idx.raw = [];
    meanConfDerivativePeak.(groups{g}).idx.raw = [];
    M_ratioDerivative.(groups{g}).all.raw = [];
    metaDaDerivative.(groups{g}).all.raw = [];
    meanConfDerivative.(groups{g}).all.raw = [];
    M_ratioDerivative.(groups{g}).trainedSubjects.raw = [];
    metaDaDerivative.(groups{g}).trainedSubjects.raw = [];
    meanConfDerivative.(groups{g}).trainedSubjects.raw = [];
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        M_ratio.(groups{g}).(session).raw = [];
        metaDa.(groups{g}).(session).raw = [];
        meanConf.(groups{g}).(session).raw = [];
    end
end

subjects = fieldnames(results);
for sub = 1:numel(subjects)
    group = sprintf('group_%d', results.(subjects{sub}).group);
    numSubjects.(group) = numSubjects.(group) + 1;
    M_ratioDerivative.(group).(subjects{sub}) = [0];
    metaDaDerivative.(group).(subjects{sub}) = [0];
    meanConfDerivative.(group).(subjects{sub}) = [0];
    for sesh = allSesh
        session = sprintf('session_%.2d',sesh);
        if sesh > 1
            lastSession = sprintf('session_%.2d', lastSesh);
            M_ratioDerivative.(group).(subjects{sub}) = horzcat(M_ratioDerivative.(group).(subjects{sub}), results.(subjects{sub}).(session).perception.trained.fit.M_ratio - results.(subjects{sub}).(lastSession).perception.trained.fit.M_ratio);
            metaDaDerivative.(group).(subjects{sub}) = horzcat(metaDaDerivative.(group).(subjects{sub}), results.(subjects{sub}).(session).perception.trained.fit.meta_da - results.(subjects{sub}).(lastSession).perception.trained.fit.meta_da);
            meanConfDerivative.(group).(subjects{sub}) = horzcat(meanConfDerivative.(group).(subjects{sub}), results.(subjects{sub}).(session).perception.trained.meanConf - results.(subjects{sub}).(lastSession).perception.trained.meanConf);
        end
        M_ratio.(group).(session).raw = vertcat(M_ratio.(group).(session).raw, results.(subjects{sub}).(session).perception.trained.fit.M_ratio);
        metaDa.(group).(session).raw = vertcat(metaDa.(group).(session).raw, results.(subjects{sub}).(session).perception.trained.fit.meta_da);
        meanConf.(group).(session).raw = vertcat(meanConf.(group).(session).raw, results.(subjects{sub}).(session).perception.trained.meanConf); 
        lastSesh = sesh;
    end
    M_ratioDerivative.(group).all.raw = vertcat(M_ratioDerivative.(group).all.raw, M_ratioDerivative.(group).(subjects{sub}));
    metaDaDerivative.(group).all.raw = vertcat(metaDaDerivative.(group).all.raw, metaDaDerivative.(group).(subjects{sub}));
    meanConfDerivative.(group).all.raw = vertcat(meanConfDerivative.(group).all.raw, meanConfDerivative.(group).(subjects{sub}));
    
    [M_ratioDerivativePeak.(group).(subjects{sub}).val, M_ratioDerivativePeak.(group).(subjects{sub}).idx] = max(M_ratioDerivative.(group).(subjects{sub}));
    M_ratioDerivativePeak.(group).idx.raw = vertcat(M_ratioDerivativePeak.(group).idx.raw, M_ratioDerivativePeak.(group).(subjects{sub}).idx);
    [metaDaDerivativePeak.(group).(subjects{sub}).val, metaDaDerivativePeak.(group).(subjects{sub}).idx] = max(metaDaDerivative.(group).(subjects{sub}));
    metaDaDerivativePeak.(group).idx.raw = vertcat(metaDaDerivativePeak.(group).idx.raw, metaDaDerivativePeak.(group).(subjects{sub}).idx);
    [meanConfDerivativePeak.(group).(subjects{sub}).val, meanConfDerivativePeak.(group).(subjects{sub}).idx] = max(meanConfDerivative.(group).(subjects{sub}));
    meanConfDerivativePeak.(group).idx.raw = vertcat(meanConfDerivativePeak.(group).idx.raw, meanConfDerivativePeak.(group).(subjects{sub}).idx);
end

for g = 1:numel(groups)
    M_ratioDerivative.(groups{g}).mean = [0];
    metaDaDerivative.(groups{g}).mean = [0];
    meanConfDerivative.(groups{g}).mean = [0];
    M_ratioDerivative.(groups{g}).all.mean = nanmean(M_ratioDerivative.(groups{g}).all.raw,1);
    M_ratioDerivative.(groups{g}).all.sem = nanstd(M_ratioDerivative.(groups{g}).all.raw,1)/sqrt(numSubjects.(groups{g}));
    metaDaDerivative.(groups{g}).all.mean = nanmean(metaDaDerivative.(groups{g}).all.raw,1);
    metaDaDerivative.(groups{g}).all.sem = nanstd(metaDaDerivative.(groups{g}).all.raw,1)/sqrt(numSubjects.(groups{g}));
    meanConfDerivative.(groups{g}).all.mean = nanmean(meanConfDerivative.(groups{g}).all.raw,1);
    meanConfDerivative.(groups{g}).all.sem = nanstd(meanConfDerivative.(groups{g}).all.raw,1)/sqrt(numSubjects.(groups{g}));
    M_ratioDerivativePeak.(groups{g}).idx.mean = nanmean(M_ratioDerivativePeak.(groups{g}).idx.raw,1);
    M_ratioDerivativePeak.(groups{g}).idx.sem = nanstd(M_ratioDerivativePeak.(groups{g}).idx.raw,[],1)/sqrt(numSubjects.(groups{g}));
    metaDaDerivativePeak.(groups{g}).idx.mean = nanmean(metaDaDerivativePeak.(groups{g}).idx.raw,1);
    metaDaDerivativePeak.(groups{g}).idx.sem = nanstd(metaDaDerivativePeak.(groups{g}).idx.raw,[],1)/sqrt(numSubjects.(groups{g}));
    meanConfDerivativePeak.(groups{g}).idx.mean = nanmean(meanConfDerivativePeak.(groups{g}).idx.raw,1);
    meanConfDerivativePeak.(groups{g}).idx.sem = nanstd(meanConfDerivativePeak.(groups{g}).idx.raw,[],1)/sqrt(numSubjects.(groups{g}));
    for sesh = allSesh
        session = sprintf('session_%.2d', sesh);
        M_ratio.(groups{g}).(session).mean = nanmean(M_ratio.(groups{g}).(session).raw,1);
        M_ratio.(groups{g}).(session).sem = nanstd(M_ratio.(groups{g}).(session).raw,[],1)/sqrt(numSubjects.(groups{g}));
        metaDa.(groups{g}).(session).mean = nanmean(metaDa.(groups{g}).(session).raw,1);
        metaDa.(groups{g}).(session).sem = nanstd(metaDa.(groups{g}).(session).raw,[],1)/sqrt(numSubjects.(groups{g}));
        meanConf.(groups{g}).(session).mean = nanmean(meanConf.(groups{g}).(session).raw,1);
        meanConf.(groups{g}).(session).sem = nanstd(meanConf.(groups{g}).(session).raw,[],1)/sqrt(numSubjects.(groups{g}));
        if sesh > 1 % Sessions 2-10||100
            lastSession = sprintf('session_%.2d', lastSesh);
            M_ratioDerivative.(groups{g}).mean = horzcat(M_ratioDerivative.(groups{g}).mean, M_ratio.(groups{g}).(session).mean - M_ratio.(groups{g}).(lastSession).mean);
            metaDaDerivative.(groups{g}).mean = horzcat(metaDaDerivative.(groups{g}).mean, metaDa.(groups{g}).(session).mean - metaDa.(groups{g}).(lastSession).mean);
            meanConfDerivative.(groups{g}).mean = horzcat(meanConfDerivative.(groups{g}).mean, meanConf.(groups{g}).(session).mean - meanConf.(groups{g}).(lastSession).mean);
        end
        lastSesh = sesh;
    end
end

if plotFigs

    strategyMRatioTimeCoursePlot = figure;
    set(gcf, 'position', [200 200 450 300]);
    lineProps = struct('color',{'cm'},'width',3,'edgestyle',':');
    H = mseb(1:10, [meanConfDerivative.group_2.all.mean', M_ratioDerivative.group_2.all.mean']', [meanConfDerivative.group_2.all.sem', M_ratioDerivative.group_2.all.sem']',lineProps,0); hold on;
    plot([1 10], [0 0], '--k', 'linewidth', 1);
    set(gca,'xtick', 1:10, 'xticklabel', {'Pre', 2:9, 'Post'});
    xlim([0 11]); 
    set(gca,'fontsize',14);
    title('Experimental Group');
    xlabel('Session'); ylabel('Rate of Change');
    ylim([ -.5 1]);
    legend('Confidence', 'meta-d''/d''', 'location', 'nw');
    legend boxoff; box off;
    if exportFigs
        export_fig(fullfile(figDir, 'Figure6A.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', strategyMRatioTimeCoursePlot)
    end

    strategyMRatioPeakBar = figure;
    set(gcf, 'position', [200 200 450 300]);
    barwitherr([meanConfDerivativePeak.group_2.idx.sem, M_ratioDerivativePeak.group_2.idx.sem], [meanConfDerivativePeak.group_2.idx.mean, M_ratioDerivativePeak.group_2.idx.mean]);
    set(gca,'fontsize',14);
    set(gca,'xticklabel',{'Confidence', 'meta-d''/d'''}, 'ytick', 2:10, 'yticklabel', {2:9, 'Post'});
    title('Experimental Group');
    xlabel('Measures'); ylabel('Session Peak');
    xlim([0 3]); ylim([1 10]);
    [h, p, ci, stats] = ttest(meanConfDerivativePeak.group_2.idx.raw, M_ratioDerivativePeak.group_2.idx.raw);
    box off;
    if exportFigs
        export_fig(fullfile(figDir, 'Figure6B.pdf'), '-pdf', '-transparent', '-painters', '-nocrop', strategyMRatioPeakBar)
    end
    
end
end