%{
I want the user to be able to select a branch from all structs which are 
dde branches. Try getting everything with s = whos() then find all structs.
I want to be able to do this without having to resave the branch. It's
easier for telling the user info later.

%% Pick a branch
selection = inputParser;
selection.addParameter(

selection = input('Please select a branch by typing the name\n');
%}
%% Select a branch
selection = branch_stst;
nunst_selection = nunst_branch_stst;


%% Select a point along the branch

[ptNum, selectPlot] = locate_along_branch(selection, param, ...
    'nunst_color', nunst_selection);
disp('You have selected point: ')
disp(ptNum)

set(selectPlot,'PaperType','a4')
set(selectPlot,'PaperOrientation','landscape');
set(selectPlot,'PaperUnits','normalized');
set(selectPlot,'PaperPosition', [0 0 1 1]);
selectPlotFileName = [master_options.datadir_specific,'selectPlot.pdf'];
print(selectPlot,selectPlotFileName,'-dpdf')

%% Create a time series there

[~, solverPlot] = ...
    timeSeries_atBranchPt(selection.point(ptNum), ...
    [0, 1000], ...
    param, ...
    master_options, ...
    'plot', 1, ...
    'save', 0, ...
    'dde23_options', ddeset('RelTol',10^-10)); % ,'OutputFcn',@odeplot,'AbsTol',1e-9

set(solverPlot,'PaperType','a4')
set(solverPlot,'PaperOrientation','landscape');
set(solverPlot,'PaperUnits','normalized');
set(solverPlot,'PaperPosition', [0 0 1 1]);
solverPlotFileName = [master_options.datadir_specific,'solverPlot.pdf'];
print(solverPlot,solverPlotFileName,'-dpdf')



%% Combine and delete the old ones

unix(['pdftk ', selectPlotFileName,' ', solverPlotFileName, ' ', ...
'cat output ', ... 
master_options.datadir_specific, 'Selection_wSolver.pdf']);
unix(['rm ',selectPlotFileName]);
unix(['rm ',solverPlotFileName]);




