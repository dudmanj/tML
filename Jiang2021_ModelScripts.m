% Example simulation scripts that implement tML 

%% STF MODEL CODE FOR EXAMPLE SESSION FIGURE

clear *_pop_avg *_pop_all mu* om* targe*
cnt=0;
intercept = 'peak';

% Estimated optimal parameters from grid search
alphas = [ 0.8 ];
sig_speed = [ 20 ];

for alpha = alphas
    for sig_spd = sig_speed
        cnt = cnt+1;

                for qqq=1:15

                clear pos traj
                arena.dims = [-300 300 0 600];
                blocks = [80 60];

                arena.x = arena.dims([1 2 2 1]);
                arena.y = arena.dims([3 3 4 4]);
                smth_win = 20;
                target.cntr = [0 280];
                target.width = [50 50];
                target.x = [target.cntr(1)-target.width(1) target.cntr(1)+target.width(1) target.cntr(1)+target.width(1) target.cntr(1)-target.width(1)];
                target.y = [target.cntr(2)-target.width(2) target.cntr(2)-target.width(2) target.cntr(2)+target.width(2) target.cntr(2)+target.width(2)];

                target2.cntr = [0 380];
                target2.width = [50 50];
                target2.x = [target2.cntr(1)-target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)-target2.width(1)];
                target2.y = [target2.cntr(2)-target2.width(2) target2.cntr(2)-target2.width(2) target2.cntr(2)+target2.width(2) target2.cntr(2)+target2.width(2)];

                % Refining the noise model for the dilation parameter
                tmp = 2*randn(1,sum(blocks));
                clear tmp2;
                for jj=1:sum(blocks)
                tmp2(jj) = randperm(3,1)+1;
                end
                tmp_om = 0.1*randn(1,sum(blocks));

                tmpC = conv(tmp,[-ones(1,4) ones(1,4)]./8,'same'); % mice seem to run a longer trajectory after a short one and vice versa
                
                % draw arena
                figure(10); clf;
                        subplot(131); 
                        patch(arena.x,arena.y,[0.95 0.9 0.95]); hold on;
                        patch(target.x,target.y,[0.85 0.9 0.9]);
                        patch(target2.x,target2.y,[0.9 0.9 0.85]);

                mu_spd(1) = 330; 
                om_hd(1) = 0;
                beta=0.01;
                target_hit = zeros(1,sum(blocks)); 
                hit_time = ones(1,sum(blocks)); rew_rate=zeros(1,100); 


                % run a trial
                for k=1:sum(blocks) %:100

                    % Calculate the local reward rate
                    if k<=10
                         rew_rate(k) = 0.25 ;
                    else
                         rew_rate(k) = 0.1 + (1 - sum([target_hit(k-3).*0.1 target_hit(k-2).*0.3 target_hit(k-1).*0.6]));    
                    end

                    % Generate the speed profile as a function of time
                    speed = [TNC_CreateGaussian(250,50+randn(1),500,1) TNC_CreateGaussian(250,50+randn(1),500,1)];
                    if sign(tmp(k))>0
                        speed = (mu_spd(k) + (rew_rate(k) .* tmp2(k) .* sig_spd)) .* speed;
                    else
                        speed = (mu_spd(k) + (rew_rate(k) .* tmp2(k)./-1.75 .* sig_spd)) .* speed;
                    end

                    % Generate the heading profile as a function of time
                    heading = [0.001:2*pi/1000:2*pi] + sgolayfilt( conv( randn(size(speed)) , [0 ones(1,smth_win)./smth_win 0] , 'same' ) , 3 , 201 );
                    heading = heading + om_hd(k) + (rew_rate(k) .* tmp_om(k));

                    % Initialize variables
                    pos.x(1)=randperm(60,1)-30; pos.y(1)=0;

                    for t=2:1000

                        [dx,dy] = pol2cart(heading(t),speed(t));

                        pos.x(t) = pos.x(t-1)+dx;

                            if pos.x(t)<arena.dims(1)
                                pos.x(t)=arena.dims(1);                
                            elseif pos.x(t)>arena.dims(2)
                                pos.x(t)=arena.dims(2);
                            end

                        pos.y(t) = pos.y(t-1)+dy;

                            if pos.y(t)<arena.dims(3)
                                pos.y(t)=arena.dims(3);                
                            elseif pos.y(t)>arena.dims(4)
                                pos.y(t)=arena.dims(4);
                            end

                        % check if the target was intercepted
                        if k>blocks(1)
                            curr_target = target2;
                        else
                            curr_target = target;            
                        end

                        switch intercept % Our experimental data has been run with different detection rules, simulating either here
                            case 'peak'
                                if t==500
                                    if pos.y(t)>curr_target.cntr(2)-curr_target.width(2) & pos.y(t)<curr_target.cntr(2)+curr_target.width(2) & pos.x(t)>curr_target.cntr(1)-curr_target.width(1) & pos.x(t)<curr_target.cntr(1)+curr_target.width(1)
                                        target_hit(k) = 1;
                                        if hit_time(k)==1
                                            hit_time(k) = t;
                                        end
                                    end
                                end
                                
                            case 'anypt'
                                if pos.y(t)>curr_target.cntr(2)-curr_target.width(2) & pos.y(t)<curr_target.cntr(2)+curr_target.width(2) & pos.x(t)>curr_target.cntr(1)-curr_target.width(1) & pos.x(t)<curr_target.cntr(1)+curr_target.width(1)
                                    target_hit(k) = 1;
                                    if hit_time(k)==1
                                        hit_time(k) = t;
                                    end
                                end
                        end

                    end                   
                    
                    % Calculate the update deltas for speed (mu) and heading (om)
                    if k>1 & target_hit(k)==1

                        if sign(tmp(k))>0
                            d_mu = alpha .* (rew_rate(k) .* tmp2(k) .* sig_spd);
                        else
                            d_mu = alpha .* (rew_rate(k) .* tmp2(k)./-1.75 .* sig_spd);
                        end        

                        d_mu_h = ( beta .* (mu_spd(k)-mu_spd(1)) );
                        mu_spd(k+1) = mu_spd(k) + d_mu - d_mu_h;
                        if mu_spd(k+1) > 500
                            mu_spd(k+1) = 500;
                        end

                        d_om = alpha .* (rew_rate(k) .* tmp_om(k));
                        d_om_h = ( beta .* om_hd(k) );
                        om_hd(k+1) = om_hd(k) + d_om - d_om_h;

                    else % Can set an exploration term here, but most parsimonious is static

                        mu_spd(k+1) = mu_spd(k);
                        om_hd(k+1) = om_hd(k);

                    end

                    traj.x(k,:) = pos.x;
                    traj.y(k,:) = pos.y;

                    if k<=blocks(1)
                        subplot(131); 
                        plot(pos.x,pos.y,'LineWidth',2,'color',[0.5 0 0 0.025]); hold on;
                        plot(pos.x(hit_time(k)),pos.y(hit_time(k)),'*','LineWidth',2,'color',[0.5 0 0]); hold on;        
                    else
                        subplot(131); 
                        plot(pos.x,pos.y,'LineWidth',2,'color',[0 0 0.5 0.025]); hold on;
                        plot(pos.x(hit_time(k)),pos.y(hit_time(k)),'*','LineWidth',2,'color',[0 0 0.5]); hold on;        
                    end
                    subplot(1,3,2); hold off; plot(mu_spd(1:k),'o-','color',[0.5 0.5 0.5]); hold on; scatter(1:k,mu_spd(1:k),50,[1:k]>blocks(1),'filled'); colormap([ 1 0 0 ; 0 0 1 ]);
                    ylabel('Speed (mu)'); box off;

                    subplot(133); hold off; plot(rew_rate(1:k),'-');
                    ylabel('Reward rate (norm.)');  box off;
                    if qqq==1
                        drawnow;
                    end
                end

                transition(qqq) = blocks(1)+1;
                mu_pop_avg(qqq,:) = mu_spd;
                corr_pop_avg(qqq,:) = target_hit;

                end
                
                
scale = 1;                

% Plot outout to allow examination of simulation results
figure(302); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
subplot(121);
shadedErrorBar(1:blocks(1),mean(mu_pop_avg(:,1:blocks(1))),std(mu_pop_avg(:,1:blocks(1)))./sqrt(size(mu_pop_avg,1))./scale,{'color',cat_map(1,:)}); hold on;
shadedErrorBar(blocks(1)+1:size(mu_pop_avg,2),mean(mu_pop_avg(:,blocks(1)+1:size(mu_pop_avg,2))),std(mu_pop_avg(:,blocks(1)+1:size(mu_pop_avg,2)))./sqrt(size(mu_pop_avg,1))./scale,{'color',cat_map(2,:)}); hold on;
ylabel('Scaled amplitude');
xlabel('Trials');

subplot(122);
shadedErrorBar(1:blocks(1),1-mean(corr_pop_avg(:,1:blocks(1))),std(corr_pop_avg(:,1:blocks(1)))./sqrt(size(corr_pop_avg,1))./scale,{'color',cat_map(1,:)}); hold on;
shadedErrorBar(blocks(1)+1:size(corr_pop_avg,2),1-mean(corr_pop_avg(:,blocks(1)+1:size(corr_pop_avg,2))),std(corr_pop_avg(:,blocks(1)+1:size(corr_pop_avg,2)))./sqrt(size(corr_pop_avg,1))./scale,{'color',cat_map(2,:)}); hold on;
ylabel('Errors per trial');
xlabel('Trials');

figure(303); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
subplot(121);
for jj=1:blocks(1)
    plot(traj.x(jj,:),traj.y(jj,:),'color',[cat_map(1,:) 0.1]); hold on;
end
plot(mean(traj.x(1:blocks(1),:)),mean(traj.y(1:blocks(1),:)),'color',cat_map(1,:),'linewidth',2);
ylabel('Y');
xlabel('X');
axis([-300 300 0 600]);


subplot(122);
for jj=blocks(1)+1:size(traj.x,1)
    plot(traj.x(jj,:),traj.y(jj,:),'color',[cat_map(2,:) 0.1]); hold on;
end
plot(mean(traj.x(blocks(1)+1:end,:)),mean(traj.y(blocks(1)+1:end,:)),'color',cat_map(2,:),'linewidth',2);
ylabel('Y');
xlabel('X');
axis([-300 300 0 600]);

    end
end

%% NTF MODEL CODE FOR EXAMPLE SESSION FIGURE
clear *_pop_avg *_pop_all *_spd path *_hit
cnt=0;
corr_w=25;
detect_style = 'inter';

% Optimal parameter estimates from grid search comparison to experimental data
alphas = [ 0.6 ];
sig_speed = [ 40 ];

for alpha = alphas
    for sig_spd = sig_speed
        cnt = cnt+1;

                for qqq=1:15

                clear pos traj
                arena.dims = [-300 300 0 600];
                blocks = [40 40 40 40];

                arena.x = arena.dims([1 2 2 1]);
                arena.y = arena.dims([3 3 4 4]);
                smth_win = 20;
                
                target.cntr = [0 250];
                target.width = [300 50];
                target.x = [target.cntr(1)-target.width(1) target.cntr(1)+target.width(1) target.cntr(1)+target.width(1) target.cntr(1)-target.width(1)];
                target.y = [target.cntr(2)-target.width(2) target.cntr(2)-target.width(2) target.cntr(2)+target.width(2) target.cntr(2)+target.width(2)];

                target2.cntr = [0 450];
                target2.width = [300 50];
                target2.x = [target2.cntr(1)-target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)-target2.width(1)];
                target2.y = [target2.cntr(2)-target2.width(2) target2.cntr(2)-target2.width(2) target2.cntr(2)+target2.width(2) target2.cntr(2)+target2.width(2)];

                target3.cntr = [0 350];
                target3.width = [300 50];
                target3.x = [target3.cntr(1)-target3.width(1) target3.cntr(1)+target3.width(1) target3.cntr(1)+target3.width(1) target3.cntr(1)-target3.width(1)];
                target3.y = [target3.cntr(2)-target3.width(2) target3.cntr(2)-target3.width(2) target3.cntr(2)+target3.width(2) target3.cntr(2)+target3.width(2)];

                target4.cntr = [0 250];
                target4.width = [300 50];
                target4.x = [target4.cntr(1)-target4.width(1) target4.cntr(1)+target4.width(1) target4.cntr(1)+target4.width(1) target4.cntr(1)-target4.width(1)];
                target4.y = [target4.cntr(2)-target4.width(2) target4.cntr(2)-target4.width(2) target4.cntr(2)+target4.width(2) target4.cntr(2)+target4.width(2)];
                
                % Refining the noise model for the dilation parameter
                tmp = 2*randn(1,sum(blocks));
                clear tmp2;
                for jj=1:sum(blocks)
                    tmp2(jj) = randperm(5,1)+1;
                end
                tmp_om = 0.1*randn(1,sum(blocks));

                tmpC = conv(tmp,[-ones(1,4) ones(1,4)]./8,'same'); % mice seem to run a longer trajectory after a short one and vice versa

                mu_spd(1) = randperm(100,1)+275; 
                om_hd(1) = 0;
                beta=0.01;
                target_hit = zeros(1,sum(blocks)); 
                hit_time = ones(1,sum(blocks)); rew_rate=zeros(1,100); 


                % run a trial
                for k=1:sum(blocks) %:100

                    % estimate reward rate
                    if k<=10
                         rew_rate(k) = 0.25 ;
                    else
                         rew_rate(k) = 0.1 + (1 - sum([target_hit(k-3).*0.1 target_hit(k-2).*0.3 target_hit(k-1).*0.6]));    
                    end

                    % calc speed profile
                    speed = TNC_CreateGaussian(250,50+randn(1),500,1);
                    if sign(tmp(k))>0
                        speed = (mu_spd(k) + (rew_rate(k) .* tmp2(k) .* sig_spd)) .* speed;
                    else
                        speed = (mu_spd(k) + (rew_rate(k) .* tmp2(k)./-1.75 .* sig_spd)) .* speed;
                    end

                    % calc heading profile
                    heading = [0.001:2*pi/1000:pi] + sgolayfilt( conv( randn(size(speed)) , [0 ones(1,smth_win)./smth_win 0] , 'same' ) , 3 , 201 );
                    heading = heading + om_hd(k) + (rew_rate(k) .* tmp_om(k));

                    % initialize some variables for sim
                    pos.x(1)=randperm(60,1)-30; pos.y(1)=0;
                   
                    for t=2:500

                        [dx,dy] = pol2cart(heading(t),speed(t));

                        pos.x(t) = pos.x(t-1)+dx;

                            if pos.x(t)<arena.dims(1)
                                pos.x(t)=arena.dims(1);                
                            elseif pos.x(t)>arena.dims(2)
                                pos.x(t)=arena.dims(2);
                            end

                        pos.y(t) = pos.y(t-1)+dy;

                            if pos.y(t)<arena.dims(3)
                                pos.y(t)=arena.dims(3);                
                            elseif pos.y(t)>arena.dims(4)
                                pos.y(t)=arena.dims(4);
                            end

                    end
                    
                    path(k) = trapz(speed);
                    
                    switch find(k<=cumsum(blocks),1,'first')
                       case 1                           
                           curr_target = target;                           
                       case 2
                           curr_target = target2;                           
                       case 3
                           curr_target = target3;
                       case 4
                           curr_target = target4;
                    end
                    

                    switch detect_style
                        
                        case 'endpnt' %  (1) does trajectory end in target [idealized task design] 
                            if max(pos.y)>curr_target.cntr(2)-curr_target.width(2) & max(pos.y)<curr_target.cntr(2)+curr_target.width(2) & pos.x(find(pos.y==max(pos.y),1))>curr_target.cntr(1)-curr_target.width(1) & pos.x(find(pos.y==max(pos.y),1))<curr_target.cntr(1)+curr_target.width(1)
                                target_hit(k) = 1;
                                if hit_time(k)==1
                                    hit_time(k) = find(pos.y==max(pos.y),1);
                                end
                            end

                        case 'inter'    %  (2) does trajectory intercept for a little while [closer to actual design]
                            in_target = find(pos.y>curr_target.cntr(2)-curr_target.width(2) & pos.y<curr_target.cntr(2)+curr_target.width(2) & pos.x>curr_target.cntr(1)-curr_target.width(1) & pos.x<curr_target.cntr(1)+curr_target.width(1));
                            if numel(in_target)>100
                                target_hit(k) = 1;
                                if hit_time(k)==1
                                    hit_time(k) = in_target(1);
                                end
                            end
                     
                    end
                    
                    if k>1 & target_hit(k)==1

                        if sign(tmp(k))>0
                            d_mu = alpha .* (rew_rate(k) .* tmp2(k) .* sig_spd);
                        else
                            d_mu = alpha .* (rew_rate(k) .* tmp2(k)./-1.75 .* sig_spd);
                        end        

                        d_mu_h = ( beta .* (mu_spd(k)-mu_spd(1)) );
                        mu_spd(k+1) = mu_spd(k) + d_mu - d_mu_h;
                        if mu_spd(k+1) > 500
                            mu_spd(k+1) = 500;
                        end

                        d_om = alpha .* (rew_rate(k) .* tmp_om(k));
                        d_om_h = ( beta .* om_hd(k) );
                        om_hd(k+1) = om_hd(k) + d_om - d_om_h;

                    else

                        mu_spd(k+1) = mu_spd(k);
                        om_hd(k+1) = om_hd(k);

                    end

                    traj.x(k,:) = pos.x;
                    traj.y(k,:) = pos.y;
                    act_dist(k) = max(pos.y);

                end

                transition(qqq) = blocks(1)+1;
                mu_pop_avg(qqq,:) = mu_spd;
                pt_pop_avg(qqq,:) = path;
                corr_pop_avg(qqq,:) = target_hit;
                act_move(qqq,:) = act_dist;
                disp(find(target_hit==1 & [1:sum(blocks)]>60,1))

                end
                

figure(402); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
subplot(121);
for jj=1:4
    inds = [1:40]+(40*(jj-1));
    shadedErrorBar(inds,mean(pt_pop_avg(:,inds)),std(pt_pop_avg(:,inds))./sqrt(size(pt_pop_avg,1))./scale,{'color',cat_map(jj,:)}); hold on;
end
ylabel('Scaled amplitude');
xlabel('Trials');

subplot(122);
for jj=1:4
    inds = [1:40]+(40*(jj-1));
    shadedErrorBar(inds,1-mean(corr_pop_avg(:,inds)),std(corr_pop_avg(:,inds))./sqrt(size(corr_pop_avg,1))./scale,{'color',cat_map(jj,:)}); hold on;
end
ylabel('Errors per trial');
xlabel('Trials');

figure(403); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
for kk=1:4
    subplot(2,2,kk);
    inds = [1:40]+(40*(kk-1));
    for jj=inds
        plot(traj.x(jj,:),traj.y(jj,:),'color',[cat_map(kk,:) 0.1]); hold on;
    end
    plot(mean(traj.x(inds,:)),mean(traj.y(inds,:)),'color',cat_map(kk,:),'linewidth',2);
    ylabel('Y');
    xlabel('X');
    axis([-300 300 0 600]);
end

    end
end
