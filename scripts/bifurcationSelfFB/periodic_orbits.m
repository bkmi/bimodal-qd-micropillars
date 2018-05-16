%% Find the hopf bifurcations which destabilize a strDomStst branch 
destabilizing_hopfs_StrDom = struct;

% only the destabilzing bifurcations
for i = 1:numel(hopfPkgStrDom)
    % select the ones that destablize the branch
    [indDest,destabFold,destabHopf] = destabilize_finder(branchPhaseStrDom(i,1).nunst);
    
    if numel(destabHopf) ~= 0 && max(destabHopf) <= numel(hopfPkgStrDom{i})
        for j = destabHopf' % 1:numel(hopfPkgStrDom{i})
            if hopfPkgStrDom{i}(j,1).error == 0
                destabilizing_hopfs_StrDom(end).hopf = hopfPkgStrDom{i}(j,1);
                destabilizing_hopfs_StrDom(end).stst = branchPhaseStrDom(i,1);
                destabilizing_hopfs_StrDom(end).indDest = indDest;
                destabilizing_hopfs_StrDom(end+1).hopf = [];
            end
        end
    end
    
    if numel(destabFold) ~= 0 && max(destabFold) <= numel(foldPkgStrDom{i})
        for j = destabFold' % 1:numel(foldPkgStrDom{i})
            if foldPkgStrDom{i}(j,1).error == 0
                disp(i)
                disp(j)
            end
        end
    end    
end

%% Relative Periodic Orbit

[rw_phas_per, suc]=SetupPsol(funcs, ...
    destabilizing_hopfs_StrDom(1).stst, ...
    destabilizing_hopfs_StrDom(1).indDest(1), ...
    opt_inputs{:}, ...
    'print_residual_info', 0, ...
    'radius', 0.02);

disp(5)

if ~suc
    error('Hopf initialization failed');
end

% figure;clf
rw_phas_per=br_contn(funcs,rw_phas_per,120);
br_rvers(rw_phas_per);
rw_phas_per=br_contn(funcs,rw_phas_per,120);

%% stability

[rw_phas_per_nunst,dom_per,defect,rw_phas_per.point]=GetStability(rw_phas_per,...
    'exclude_trivial',true,'locate_trivial',@(p)[1,1],'funcs',funcs);
% plot with stability info
pp=arrayfun(@(p)p.parameter(paramStrDom.feed_phase.index),rw_phas_per.point);
Epow=@(x)sqrt(sum(x(1:2,:).^2,1));
Apmx=arrayfun(@(p)max(Epow(p.profile)),rw_phas_per.point);
Apmn=arrayfun(@(p)min(Epow(p.profile)),rw_phas_per.point);
figure;hold on
sel=@(x,i)x(rw_phas_per_nunst==i);
plot(sel(pp,0),[sel(Apmx,0);sel(Apmn,0)],'ko',...
    sel(pp,1),[sel(Apmx,1);sel(Apmn,1)],'ro',...
    sel(pp,2),[sel(Apmx,2);sel(Apmn,2)],'ro');
axis tight