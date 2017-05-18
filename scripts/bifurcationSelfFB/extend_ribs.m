%% Extend certain phase 'ribs'
% Pick a few, useful phase ribs and continue them really, really far.
% First we plot to identify the useful loops



% % find highest nunst
% globalmax = 0;
% for i = 1: numel(branchPhaseStrDom)
%     localmax = max(branchPhaseStrDom(i,1).nunst);
%     if localmax > globalmax
%         globalmax = localmax;
%     end
% end
% % for j = 1:numel(ampPkgStrDom)
% %     branchAmpStrDom = ampPkgStrDom{j,1};
% %     for i = 1: numel(branchAmpStrDom)
% %         localmax = max(branchAmpStrDom(i,1).nunst);
% %         if localmax > globalmax
% %             globalmax = localmax;
% %         end
% %     end
% % end
% for i = 1: numel(branchStstStrDom)
%     localmax = max(branchStstStrDom(i,1).nunst);
%     if localmax > globalmax
%         globalmax = localmax;
%     end
% end
% 
% ribs
% for i = 1:numel(branchPhaseStrDom)
%     plot_branch3(branchPhaseStrDom(i,1), paramStrDom, ...
%         'nunst_color',{branchPhaseStrDom(i,1).nunst,  globalmax}, ...
%         'add_2_gcf',1, ...
%         'axes_indParam', ...
%         {paramStrDom.feed_phase.index, ...
%         paramStrDom.feed_ampli.index, ...
%         'x1'}, ...
%         'PlotStyle', { 'LineStyle', 'None', 'Marker', '.' } );
% end
% 
% %% Interesting ribs are numbered 19, 20, 21, 22
% figure();
% 
% ribs
% for i = [19,20,21,22]
%     plot_branch3(branchPhaseStrDom(i,1), paramStrDom, ...
%         'nunst_color',{branchPhaseStrDom(i,1).nunst,  globalmax}, ...
%         'add_2_gcf',1, ...
%         'axes_indParam', ...
%         {paramStrDom.feed_phase.index, ...
%         paramStrDom.feed_ampli.index, ...
%         'x1'}, ...
%         'PlotStyle', { 'LineStyle', 'None', 'Marker', '.' } );
% end

%% Continue them!!

extendedRibs = cell(4,1);
    save([datadir, ...
        'extendedRibs.mat'], ...
        'extendedRibs')
howmanypts = 3000;

for i = [19,20,21,22]
    curBranch = branchPhaseStrDom(i,1);
    
    % calculate
    curBranch = br_contn(funcs, curBranch, howmanypts);
    curBranch = br_rvers(curBranch);
    curBranch = br_contn(funcs, curBranch, howmanypts);
    
    % stability analysis
    [curBranch.nunst,~,~,curBranch.point] = GetRotStability( ...
        curBranch, ...
        funcs, 1);
    curBranch.indFold = find(abs(diff(curBranch.nunst))==1);
    curBranch.indHopf = find(abs(diff(curBranch.nunst))==2);
    curBranch.error = 0;
    
    % add to 
    extendedRibs{i,1} = curBranch;
    
    % Save Strong Dom
    save([datadir, ...
        'extendedRibs.mat'], ...
        'extendedRibs', ...
        '-append');
    
end

%% Let's do more in that region

% create container for phase stst branches
numPhaseBranches = 10;
pts4PhaseBranch = pts4PhaseBranch(18):pts4PhaseBranch(22);
branch_length = 3800;

transitionPhaseBranchesStrDom = repmat( ...
    struct( ...
    'method', struct, ...
    'parameter', struct, ...
    'point', struct, ...
    'nunst', 0, ...
    'indFold', 0, ...
    'indHopf', 0, ...
    'error', 0),...
    [numPhaseBranches,1]);

stepBoundStrDomPhase = { ...
    'step',pi/32, ...
    'max_step', ...
    [paramStrDom.feed_phase.index,1.75*pi/32, ...
    paramStrDom.omega1.index,0.1], ...
    'newton_max_iterations',10, ...
    'max_bound',[paramStrDom.feed_phase.index,20*pi], ...
    'min_bound', [paramStrDom.feed_phase.index,-20*pi], ...
    'halting_accuracy',1e-10, ...
    'minimal_accuracy',1e-8};

for i = 1:numPhaseBranches
    [branchStstStrDomPhase,~] = SetupStst( ...
        funcs, ...
        'contpar',[paramStrDom.feed_phase.index, paramStrDom.omega1.index], ...
        'corpar',[paramStrDom.omega1.index],...
        'x', branchStstStrDom.point(pts4PhaseBranch(i)).x, ...
        'parameter',branchStstStrDom.point(pts4PhaseBranch(i)).parameter,...
        opt_inputs{:},...
        stepBoundStrDomPhase{:});

    % calculate
    [branchStstStrDomPhase,~,~,~] = br_contn(funcs,branchStstStrDomPhase,branch_length);
    branchStstStrDomPhase = br_rvers(branchStstStrDomPhase);
    [branchStstStrDomPhase,~,~,~] = br_contn(funcs,branchStstStrDomPhase,branch_length);

    % stability analysis
    [branchStstStrDomPhase.nunst,~,~,branchStstStrDomPhase.point] = GetRotStability( ...
        branchStstStrDomPhase, ...
        funcs, 1);
    branchStstStrDomPhase.indFold = find(abs(diff(branchStstStrDomPhase.nunst))==1);
    branchStstStrDomPhase.indHopf = find(abs(diff(branchStstStrDomPhase.nunst))==2);
    branchStstStrDomPhase.error = 0;
    
    % add to 
    transitionPhaseBranchesStrDom(i,1) = branchStstStrDomPhase;
    
    % Save Strong Dom
    save([datadir, ...
        'extendedRibs.mat'], ...
        'transitionPhaseBranchesStrDom', ...
        '-append');
    
end
