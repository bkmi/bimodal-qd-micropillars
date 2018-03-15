function [ timeseries ] = turn_on( param, hist )
%Returns the turn on
ddesys = set_ddesys(param);

lags = param.tau_fb.value;
timeSpan = [0, 15];

timeseries = dde23( ...
    @(t,y,z)ddesys([y,z]),...
    lags,hist,timeSpan,ddeset('RelTol',10^-8));

end

