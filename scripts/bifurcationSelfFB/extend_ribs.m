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

