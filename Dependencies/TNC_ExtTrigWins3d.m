function [sink] = TNC_ExtTrigWins3d(source,indices,window)

numEvents   = round(numel(indices));
range       = -window(1,1):window(1,2);
numRange    = round(numel(range));

if indices(numEvents)+window(1,2) > size(source,2)
    numEvents = numEvents-1;
end

for i=1:numEvents
%         sink.wins(i).win_dat    = source(:,indices(i)-window(1,1):indices(i)+window(1,2));
        sink.wins(:,:,i)        = source(:,indices(i)-window(1,1):indices(i)+window(1,2));
        sink.pv(:,i)            = trapz( source(:,indices(i)-window(1,1):indices(i)+window(1,2)) , 2 );
end

% sink.win_avg = sink.wins(1).win_dat;
% win_sq = sink.wins(1).win_dat.^2;

% for i=2:numEvents
%         sink.win_avg = sink.win_avg+sink.wins(i).win_dat;
%         win_sq = win_sq+(sink.wins(i).win_dat.^2);
%         sink.win_var = ( win_sq - (sink.win_avg.^2 / i) ) ./ (i-1);
% end
% sink.win_avg = sink.win_avg./numEvents;
sink.win_avg = mean(sink.wins,3);
sink.win_sem = std(sink.wins,[],3)./sqrt(numEvents);
sink.inds = indices(1:numEvents);
sink.xvals = range;