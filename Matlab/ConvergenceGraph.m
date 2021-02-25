function ConvergenceGraph(name, graph_name)

%BaseDirectory           = 'C:\Users\nomen\Documents\Thesis\tests\';
BaseDirectory           = 'C:\Users\nomen\Documents\Thesis\tests\';
Directory               = strcat(BaseDirectory, strcat(name, '\'));
Techniques              = {'PHASE_FUNCTION_ONLY', 'BRDF_ONLY', 'HYBRID', 'LIGHT_PATHS', 'LIGHT_PATHS_OCTO', 'LIGHT_PATHS_OCTO_GRADIENT', 'TEST_SHADER', 'REJECTION_SAMPLER', 'ONE_DIRECTIONAL', 'OCTO_GRADIENT_INVERSE'}; %folder name
ID                      = 1;
TechniqueIDs            = [];
Series                  = {'PhaseFunctionOnly', 'BRDFOnly', 'Hybrid', 'LightPaths', 'LightPathsOcto', 'LightPathsOctoGradient', 'TestShader', 'RejectionSampler', 'OneDirectional', 'OctoGradientInverse'}; %graph name
SPP                     = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];
MaxSPP                  = 14;%12;%9;%11;
Errors                  = zeros(1);
Timings                 = zeros(1);
ReferenceImageFile      = dir([Directory 'reference*.png']);
ReferenceImage          = double(imread(strcat(Directory, ReferenceImageFile.name))) / 255.0;
Figure                  = figure('Name', name);
AxisY                   = axes('Parent', Figure, 'YScale', 'log', 'YMinorTick', 'on');
Markers                 = { 'x', 'v', 'square', '^', '+', '*', 'o', 'pentagram', 'diamond', '^'};  
Results                 = cell(1, 1);
%ConvergenceGraph('TD_NoScatter_ShadingTest', 'Scattering Types')
for i = 1 : length(Series)
    Results(1, i + 1) = Series(i);
end

for i = 1 : length(SPP)
    Results(i + 1, 1) = {SPP(i)};
end

box(AxisY, 'on');
hold(AxisY, 'all');
xlim(AxisY, [0 SPP(MaxSPP)]);

for i = 1 : length(Techniques)
    TechniqueDir = sprintf('%s%s\\', Directory, Techniques{i});
    
    if exist(TechniqueDir, 'dir')
        TechniqueIDs(ID) = i;
        CSV = readmatrix(sprintf('%smeasurements.csv', TechniqueDir));
        
        for j = 1 : length(SPP)           
            ImageFileName = sprintf('%simages\\%d.png', TechniqueDir, SPP(j));

            if exist(ImageFileName, 'file')
                Image       = double(imread(ImageFileName)) / 255.0;
                %Statistics  = CompareImages(Image, ReferenceImage, 1.0);
                Statistics  = CompareImagesWithName(Image, ReferenceImage, 1.0, sprintf('%simages\\diff_%d.png', TechniqueDir, SPP(j)));
                Error       = Statistics.msesum;
                
                if j <= MaxSPP
                    Errors(j, ID) = Error;
                    Timings(j, ID) = CSV(j,2) / 1000;
                end
                
                Results(j + 1, ID + 1) = {Error};
            else
                warning('File does not exist: %s', ImageFileName);
            end
        end
        
        ID = ID + 1;
    else
        warning('Directory does not exist: %s', TechniqueDir);
    end
end

Plot = plot(SPP(1:MaxSPP), Errors);

for i = 1 : length(TechniqueIDs)
    TechniqueID = TechniqueIDs(i);
    set(Plot(i), 'DisplayName', Series{TechniqueID});
    set(Plot(i), 'Marker', Markers{TechniqueID});
end

title(graph_name, 'FontSize', 12);
xlabel('Samples per pixel', 'FontSize', 10);
ylabel('Mean-Squared Error', 'FontSize', 10);
legend(AxisY, 'show');

ConvergencePdfFilePath = sprintf('%sconvergence.pdf', Directory);

TI = get(gca, 'TightInset');

set(gca, 'Position', [TI(1) TI(2) 1 - TI(3) - TI(1) 1 - TI(4) - TI(2)]);
set(gca, 'units', 'centimeters')

P   = get(gca,'Position');
TI  = get(gca, 'TightInset');

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [P(3) + TI(1) + TI(3) P(4) + TI(2) + TI(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 P(3) + TI(1) + TI(3) P(4) + TI(2) + TI(4)]);

saveas(Figure, ConvergencePdfFilePath, 'pdf');

xlswrite(sprintf('%serrors.xlsx', Directory), Results);

FigureTimings = figure('Name', append(name, '-Timings'));
%box(AxisY, 'on');
%hold(AxisY, 'all');
%xlim(AxisY, [0 SPP(MaxSPP)]);
TimingsPlot = plot(SPP(1:MaxSPP), Timings);
for i = 1 : length(TechniqueIDs)
    TechniqueID = TechniqueIDs(i);
    set(TimingsPlot(i), 'DisplayName', Series{TechniqueID});
    set(TimingsPlot(i), 'Marker', Markers{TechniqueID});
end
title(graph_name, 'FontSize', 12);
xlabel('Samples per pixel', 'FontSize', 10);
ylabel('Time in seconds', 'FontSize', 10);
legend('show', 'Location', 'northwest');
TI = get(gca, 'TightInset');

set(gca, 'Position', [TI(1) TI(2) 1 - TI(3) - TI(1) 1 - TI(4) - TI(2)]);
set(gca, 'units', 'centimeters')

P   = get(gca,'Position');
TI  = get(gca, 'TightInset');

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [P(3) + TI(1) + TI(3) P(4) + TI(2) + TI(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 P(3) + TI(1) + TI(3) P(4) + TI(2) + TI(4)]);
saveas(FigureTimings, sprintf('%stimings.pdf', Directory) , 'pdf');
