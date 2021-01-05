function figHandle = s_createCoveragePlot_averaged(RFcov, name, fieldRange, norm, n)

% plotting subroutine for rmPlotCoverage. Broken off by ras 10/2009.
% And edited for scripting by JG 05/2016
% All you have to do is feed it a 128x128 coverage

% % % Load dummy filler data for kids study:
% fileDir = '/sni-storage/kalanit/biac2/kgs/projects/Longitudinal/FMRI/Retinotopy/results/pRF_figures/Coverage/plottingVariables';
% load(fullfile(fileDir,'dummyVariables.mat'));
vfc.cmap = 'jet';
vfc.fieldRange = fieldRange;
data.X = repmat(linspace(-fieldRange, fieldRange, 128), 128, 1);
data.Y = repmat(linspace(-fieldRange, fieldRange, 128), 128, 1)';
figHandle = gcf;


% normalize the color plots to 1
if norm == 1
    rfMax = max(RFcov(:));
else
	rfMax = 1;
end

img = RFcov ./ rfMax;

mask = makecircle(length(img));
img = img .* mask;
imagesc(data.X(1, :), data.Y(:, 1), img);
set(gca, 'YDir', 'normal');
grid on

colormap(vfc.cmap);
colorbar;

% start plotting
hold on;


% add polar grid on top
p.ringTicks = (1:3) / 3 * vfc.fieldRange;
p.color = [0.6, 0.6, 0.6];
polarPlot([], p);


% scale z-axis
caxis([0,1]);
% if norm == 1
%     caxis([0, 1]);
% else
%     if min(RFcov(:))>=0
%         caxis([0 ceil(max(RFcov(:)))]);
%     else
%         caxis([-1 1]); %* ceil(max(abs(RFcov(:)))));
%     end
% end

axis image; % axis square;
xlim([-vfc.fieldRange, vfc.fieldRange])
ylim([-vfc.fieldRange, vfc.fieldRange])


title([name, sprintf('\n'), 'N = ', n], 'FontSize', 16, 'Interpreter', 'none');


return;