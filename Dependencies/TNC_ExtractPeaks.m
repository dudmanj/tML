function [events] = TNC_ExtractPeaks(signal,detect_thresh,min_duration,plotFlag)

% detect_thresh = 75;
% min_duration = 25;

signal_df = [diff(signal) 0];
starts_put = find(signal<detect_thresh & signal+signal_df>detect_thresh & signal_df>0);

if plotFlag
    figure(1); clf; 
    subplot(211);
    plot(signal); hold on;
    plot(signal_df);
    plot(starts_put,signal(starts_put),'r*');
    axis tight;
end

tmp = [0 diff(starts_put)];
tmp(1:2) = 0;
tmp(numel(tmp)-1:numel(tmp)) = 0;
starts = starts_put(find(tmp>min_duration));
events.starts = starts;

for j=1:numel(starts)
    if j==numel(starts)
        tmp2 = find(signal(starts(j)+1:numel(signal))>detect_thresh,1,'last');
        if numel(tmp2)==0
            tmp2=numel(signal);
        end
    else
        tmp2 = find(signal(starts(j):starts(j+1))>detect_thresh,1,'last');
        if numel(tmp2)==0
            tmp2=starts(j+1)-starts(j)-1;
        end
    end
    duration(j) = tmp2;
    stops(j) = starts(j)+tmp2;
    int(j) = trapz(signal(starts(j):stops(j)));
end

if plotFlag
    subplot(212);
    plot(signal,'.-'); hold on;
    plot(starts,signal(starts),'g*');
    plot(stops,signal(stops),'rs');
    axis tight;
end

events.starts = starts;
events.stops = stops;
events.dur = duration;
events.num = numel(starts);
events.int = int;
