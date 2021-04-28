function [w_av] = AvgOverRuns(ws,mode)
% averaging over inference runs of different seeds
% if mode is 'discard' then outlying runs are discarded
% otherwise all runs are included
    outliers=sum((ws-repmat(mean(ws),[size(ws,1) 1]))>repmat(std(ws),[size(ws,1) 1]),2);
    exclude=[];
    if strcmp(mode,'discard')
        exclude=(find((outliers-mean(outliers))>std(outliers)));
    end
    w_av=mean(ws(setdiff(1:end,exclude),:),1);
end

