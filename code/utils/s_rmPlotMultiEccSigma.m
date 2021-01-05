function data = s_rmPlotMultiEccSigma(vw, ROIlist, vethresh, eccthresh, sigthresh, voxthresh)
% rmPlotMultiEccSigma(vw, [ROIlist])
%
% Wrapper to plot pRF sigma vs eccentricity for multiple ROIs
%
%   vw: mrVista view struct
%   ROIlist: list of ROIs (if blank, call menu; if 0, plot all ROIs)
%
%   note: need to have an rm model and at least one ROI loaded into the
%           view.
%
% 2/2009: JW

%--------------------------
% VARIABLE CHECKS
%--------------------------
% check view struct
if notDefined('vw'), vw = getCurView; end

% check model
model = viewGet(vw, 'rmModel'); %#ok<NASGU>
if isempty('model'), vw = rmSelect(vw); end

% check ROIs
if (notDefined('ROIlist'))
    roiList = viewGet(vw, 'roinames');
    selectedROIs = find(buttondlg('ROIs to Plot', roiList));
elseif ROIlist == 0,
    selectedROIs = 1:length(viewGet(vw, 'ROIs'));
else
    selectedROIs = ROIlist;
end

nROIs = length(selectedROIs);
if (nROIs == 0), error('No ROIs selected'); end

%--------------------------
% PLOT
%--------------------------

% set up plot
graphwin = selectGraphWin;
figure(graphwin);
hold on;
set(graphwin, 'Color', 'w')
c = jet(nROIs);

% initialize a legend
legendtxt = cell(1, nROIs);

% initialize data struct
data = cell(1, nROIs);

% suppress individual plots from calls to rmPlotEccSigma
plotFlag = false;

% loop thru ROIs
for ii = 1:nROIs
    vw = viewSet(vw, 'curroi', selectedROIs(ii));
    data{ii} = s_rmPlotEccSigma(vw, [], [], [], plotFlag, vethresh, eccthresh, sigthresh, voxthresh);
    data{ii}.roi = viewGet(vw, 'roiname');
    legendtxt{ii} = data{ii}.roi;
    figure(graphwin);
    % plot the fit lines for each ROI (so we have one series per ROI to
    % make the legend nicer)
    if data{ii}.error == 0
        plot(data{ii}.xfit, data{ii}.yfit, '-', 'color', c(ii, :), 'LineWidth', 2)
    elseif data{ii}.error == 1 % If the ROI didn't have enough voxels after thresholding, plot a blank dot so we don't mess up the legend order
        plot(0, 0, 'w'); % Using white so that in the legend, it is clear the missing ROI has no fit associated with it
    end
end

legend(legendtxt);

% add the data points for each plot
for ii = 1:nROIs
    if data{ii}.error == 0
        errorbar(data{ii}.x, data{ii}.y, data{ii}.ysterr, 'x', 'color', c(ii, :));
    else
        plot(0, 0)
    end
end

ylabel('pRF size (sigma, deg)', 'fontsize', 16);
xlabel('Eccentricity (deg)', 'fontsize', 16);

data = cell2mat(data);

return
