%% strongDom, tau = 0.83 PLOT SCRIPT

branchplot = figure;
wrapplot = figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOPF
% SPLOTCH NEAR -10
for i = 1:numel(hopfTestNear025)
    % Add each hopf_branch
    if isa(hopfTestNear025(i).error,'double') && hopfTestNear025(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        figure(branchplot);
        plot_branch(hopfTestNear025(i), param, ...
            'add_2_gcf', 1, 'color','c');
        
        % wrap plot
        figure(wrapplot);
        wrapTemp = wrap_to_2pi(hopfTestNear025(i), param.feed_phase.index);
        plot_branch(wrapTemp, param, ...
            'add_2_gcf', 1, 'color','c');
        
    end
end

% STEEP LINES NEAR THE FOLDS
for i = 1:numel(hopfTestNear0311)
    % Add each hopf_branch
    if isa(hopfTestNear0311(i).error,'double') && hopfTestNear0311(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        figure(branchplot);
        plot_branch(hopfTestNear0311(i), param, ...
            'add_2_gcf', 1, 'color','c');
                
        % wrap plot
        figure(wrapplot);
        wrapTemp = wrap_to_2pi(hopfTestNear0311(i), param.feed_phase.index);
        plot_branch(wrapTemp, param, ...
            'add_2_gcf', 1, 'color','c');
        
    end
end

% ONE STEEP LINE NEAR -10, HIGH UP SPLOTCHES NEAR 0.4 FAMP
for i = 1:numel(hopfTestNear0389)
    % Add each hopf_branch
    if isa(hopfTestNear0389(i).error,'double') && hopfTestNear0389(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        figure(branchplot);
        plot_branch(hopfTestNear0389(i), param, ...
            'add_2_gcf', 1, 'color','c');
        
        % wrap plot
        figure(wrapplot);
        wrapTemp = wrap_to_2pi(hopfTestNear0389(i), param.feed_phase.index);
        plot_branch(wrapTemp, param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

% Long windy lines found from an already unstable soln at first on the
% right side
for i = 1:numel(longHopfsFoundUnstableRightSide)
    % Add each hopf_branch
    if isa(longHopfsFoundUnstableRightSide(i).error,'double') && ...
            longHopfsFoundUnstableRightSide(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        figure(branchplot);
        plot_branch(longHopfsFoundUnstableRightSide(i), param, ...
            'add_2_gcf', 1, 'color','c');
    
        % wrap plot
        figure(wrapplot);
        wrapTemp = wrap_to_2pi(longHopfsFoundUnstableRightSide(i), ...
            param.feed_phase.index);
        plot_branch(wrapTemp, param, ...
            'add_2_gcf', 1, 'color','c');
    end
end

%%%%%%%%%%%%%%%%%%%%% FOLD
% NICE LONG FOLD LINES!!
for i = 1:numel(foldNear0133)
    % Add each fold
    if isa(foldNear0133(i).error,'double') && foldNear0133(i).error == 0
        figure(branchplot);
        plot_branch(foldNear0133(i), param, ...
            'add_2_gcf', 1, 'color','r');
        
        % wrap plot
        figure(wrapplot);
        wrapTemp = wrap_to_2pi(foldNear0133(i), param.feed_phase.index);
        plot_branch(wrapTemp, param, ...
            'add_2_gcf', 1, 'color','r');
    end
end

% curved fold line on right side
% came from stdyBrnchFA_fromReg1(5).indFold(2)
figure(branchplot);
plot_branch(curveFold(1), param, ...
            'add_2_gcf', 1, 'color','r');
        
figure(wrapplot);
wrapTemp = wrap_to_2pi(curveFold(1), param.feed_phase.index);
plot_branch(wrapTemp, param, ...
    'add_2_gcf', 1, 'color','r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STST
% Plot stst branches
% for i = 1:numel(stdyBrnchFB_regSpace)
%     figure(branchplot);
%     plot_branch(stdyBrnchFB_regSpace(i), param, ...
%         'add_2_gcf', 1, 'color','g', ...
%         'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
%     
%     % wrap plot
%     figure(wrapplot);
%     wrapTemp = wrap_to_2pi(stdyBrnchFB_regSpace(i), param.feed_phase.index);
%     plot_branch(wrapTemp, param, ...
%         'add_2_gcf', 1, 'color','g', ...
%         'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
% end

%%

% vertical stst branches
for i = 2 %1:numel(stdyBrnchFA_fromReg1)
    % 1 - 7 on the right triangle, 7 pierces the fold
    % 8 is right in the middle, and it is long. => longMid
    plot_branch(stdyBrnchFA_fromReg1(i), param, ...
        'add_2_gcf', 1, 'color',[0.3 0.7 0], ...
        'axes_indParam',[param.feed_phase.index,param.feed_ampli.index]);
end

%%
plot_branch(stdyBrnchFA_fromReg1(5), param, ...
            'twoOmegaNunst', stdyBrnchFA_fromReg1(5).nunst)
        
figure(7);
        plot_branch(stdyBrnchFA_fromReg1(5), param, ...
            'twoOmegaNunst', stdyBrnchFA_fromReg1(5).nunst, ...
            'axes_indParam',[param.feed_phase.index,param.feed_ampli.index], ...
            'add_2_gcf', 1)



%% Polar
%%%%%%%%%%%%%%%%%%%%%%%% POLAR
branchplot = polar(0,0);

% Plot each testHopf
for i = 1:numel(hopfTest)
    % Add each hopf_branch
    if isa(hopfTest(i).error,'double') && hopfTest(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        plot_branch(hopfTest(i), param, ...
            'add_2_gcf', 1, 'color','c','polar',1);
    end
end

% Plot each testHopf
for i = 1:numel(hopfTest2)
    % Add each hopf_branch
    if isa(hopfTest2(i).error,'double') && hopfTest2(i).error == 0
        % Only plot hopf_branches that DO NOT have errors
        plot_branch(hopfTest2(i), param, ...
            'add_2_gcf', 1, 'color','c','polar',1);
    end
end


%%%%%%%%%%%%%%%%%%%%% FOLD POLAR
% NICE LONG FOLD LINES!!
for i = 1:numel(foldNear0133)
    % Add each fold
    if isa(foldNear0133(i).error,'double') && foldNear0133(i).error == 0
        plot_branch(foldNear0133(i), param, ...
            'add_2_gcf', 1, 'color','r','polar',1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STST POLAR
% Plot stst branches
for i = 1:numel(stdyBrnchFB_regSpace)
    
    plot_branch(stdyBrnchFB_regSpace(i), param, ...
                'add_2_gcf', 1, 'color','g', ...
                'axes_indParam',[param.feed_phase.index,param.feed_ampli.index], ...
                'polar',1);
end



