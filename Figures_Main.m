% set paths for code we need
rootdir = fileparts(which('Figures_Main')); % Directory in which we put Figures_Main.m
addpath(rootdir); % Add rootdir path
addpath(fullfile(rootdir, 'Data')); % Add data path
addpath(fullfile(rootdir, 'Figure Code')); % Add Figure Code path
addpath(fullfile(rootdir, 'helper')); % Add helper functions path

% load anonymized, preprocessed data
load('results.mat');

% set plotting variables
plotFigs = true;
exportFigs = true;
figDir = 'Figures';

% Make figDir if it doesn't already exist
if exist(figDir, 'dir') ~= 7
    mkdir(figDir);
end

% Plot and export figures to figDir
Figure2(results, plotFigs, exportFigs, figDir);
Figure3(analysis, results, plotFigs, exportFigs, figDir);
Figure4A(analysis, results, plotFigs, exportFigs, figDir);
Figure4B(results, plotFigs, exportFigs, figDir);
Figure6AB(results, plotFigs, exportFigs, figDir);

% Plot and export supplementary figures to figDir
FigureS1(results, plotFigs, exportFigs, figDir);
FigureS2(analysis, plotFigs, exportFigs, figDir);
FigureS3(analysis, results, plotFigs, exportFigs, figDir);
FigureS4A(analysis, results, plotFigs, exportFigs, figDir);
FigureS4B(results, plotFigs, exportFigs, figDir);
FigureS5A(results, plotFigs, exportFigs, figDir);
FigureS5B(results, plotFigs, exportFigs, figDir);

% load simulated data for Figure SX
load('simComparison.mat')

% Plot and export supplementary Figure SX to figDir
FigureSXA(analysis, plotFigs, exportFigs, figDir);
FigureSXB(results, plotFigs, exportFigs, figDir);
FigureSXC(analysis, results, plotFigs, exportFigs, figDir);
