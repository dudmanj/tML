%% Consolidated scripts from Jiang et al 2022, Nature neuroscience

%% Script to load STF task data from files
% NOTE: These are provided for clarity about what is extracted. The final
% extracted data structures are provided rather than original files.

cnt = 1; clear allMice;
for animal = 1:3

    anim_num = animal;
    switch animal
        case 1
            dirBase = 'OpenField_Ai93/ANM373541/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 17;
        case 2
            dirBase = 'OpenField_Ai93/ANM377489/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 17;
        case 3
            dirBase = 'OpenField_Ai93/ANM461756/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 17;
    end
    
    thresh_joy = [235 330 330];
    
    % cycle through all data from this mouse
    for jj = 1:numel(dir_list)
    clear start_trials events;

        S=load([dirBase dir_list(jj).name '/OF_Analysis_CNMFE2_t2/SyncEvt_Analysis_v2.mat']);

        x_min = min(unique(S.stExp_All.vtPos(1,:)));
        x_avg = mean(S.stExp_All.vtPos(1,:));
        y_min = min(unique(S.stExp_All.vtPos(2,:)))
        if y_min<85
            y_min=85;
        end

        % causes a slightly weird error of being hard to compare goal location with trajectories
        tmp_pos = [ S.stExp_All.vtPos(1,:)-x_avg ;  S.stExp_All.vtPos(2,:)-y_min ];

        % need to shift position values by 5 frames to account for latency difference
        length = size(tmp_pos,2);
        norm_pos = [tmp_pos(:,2:length) zeros(2,1) ];

        norm_vel = sqrt( [0 diff(norm_pos(1,:)).^2] + [0 diff(norm_pos(2,:)).^2] );
        [norm_angle,norm_amp] = cart2pol(sgolayfilt(S.stExp_All.vtPos(1,:)-x_avg,3,11),sgolayfilt(S.stExp_All.vtPos(2,:)-y_min,3,11));

        [cat_map] = TNC_CreateRBColormap(8,'mapb');
        before = 34;
        after = 51;

        % compare error trials with successes to examine SPE differences
        detect_thresh = 100;
        min_duration = 17;
        [events] = TNC_ExtractPeaks(norm_amp,detect_thresh,min_duration,1);
        title(dir_list(jj).name);

        norm_spd = [0 abs(diff(norm_amp))];
        figure(1); clf; plot(norm_amp); hold on; plot(norm_pos(2,:));  plot(norm_pos(1,:));
        
        spk_kernel = TNC_CreateGaussian(50,1,100,1);
        spk_d_kernel = TNC_CreateGaussian(50,2,100,1);

        S.tSyncEvt_V = S.tSyncEvt_V;
        spes_cont = zeros(size(norm_amp));
        spes_cont(S.tSyncEvt_V) = 1;
        spes_cont = conv(spes_cont,spk_kernel,'same');

        events.corr = zeros(1,events.num);
        events.trig = zeros(1,events.num);
        events.blk = ones(1,events.num);
        blk_frames = max(S.trial_s1);

        for kk=1:events.num
            if anim_num==3

                if kk<events.num
                    tmp = find( S.stExp_All.stEvt.Solenoid.RiseFrm > events.starts(kk) & S.stExp_All.stEvt.Solenoid.RiseFrm < events.stops(kk) , 1);
                else
                    tmp = find( S.stExp_All.stEvt.Solenoid.RiseFrm > events.starts(kk) , 1);
                end
                if numel(tmp)>0
                    events.corr(kk) = 1;
                    events.trig(kk) = S.stExp_All.stEvt.Solenoid.RiseFrm(tmp) - events.starts(kk);
                end

            else
                
                if kk<events.num
                    tmp = find( S.stExp_All.stEvt.Collect.RiseFrm > events.starts(kk) & S.stExp_All.stEvt.Collect.RiseFrm < events.stops(kk) , 1);
                else
                    tmp = find( S.stExp_All.stEvt.Collect.RiseFrm > events.starts(kk) , 1);
                end
                if numel(tmp)>0
                    events.corr(kk) = 1;
                    events.trig(kk) = S.stExp_All.stEvt.Collect.RiseFrm(tmp) - events.starts(kk);
                end

            end
            
            if events.starts(kk)>blk_frames
                events.blk(kk) = 2;
            else
                events.blk(kk) = 1;
            end
            events.thresh(kk) = thresh_joy(events.blk(kk));
        end
        
        % errors per correct
        tmp = find(events.corr==1);
        events.err_per_corr = [tmp(1) diff(tmp)]-1;
        events.err_per_blk = events.blk(events.corr==1);
        
        % move amp per event
        [amp] = TNC_ExtTrigWins(norm_amp,events.starts,[before,after]);
        [ampS] = TNC_ExtTrigWins(norm_amp,events.stops,[before,after]);
        
        if anim_num==3
            [ampT] = TNC_ExtTrigWins(norm_amp,S.stExp_All.stEvt.Solenoid.RiseFrm,[before,after]);
            [posxT] = TNC_ExtTrigWins(norm_pos(1,:),S.stExp_All.stEvt.Solenoid.RiseFrm,[before,after]);
            [posyT] = TNC_ExtTrigWins(norm_pos(2,:),S.stExp_All.stEvt.Solenoid.RiseFrm,[before,after]);        
            [caT] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Sig_df,S.stExp_All.stEvt.Solenoid.RiseFrm,[before,after]);
            [spT] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Spikes,S.stExp_All.stEvt.Solenoid.RiseFrm,[before,after]);
        else
            [ampT] = TNC_ExtTrigWins(norm_amp,S.stExp_All.stEvt.Collect.RiseFrm,[before,after]);
            [posxT] = TNC_ExtTrigWins(norm_pos(1,:),S.stExp_All.stEvt.Collect.RiseFrm,[before,after]);
            [posyT] = TNC_ExtTrigWins(norm_pos(2,:),S.stExp_All.stEvt.Collect.RiseFrm,[before,after]);            
            [caT] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Sig_df,S.stExp_All.stEvt.Collect.RiseFrm,[before,after]);
            [spT] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Spikes,S.stExp_All.stEvt.Collect.RiseFrm,[before,after]);
        end
        
        [posx] = TNC_ExtTrigWins(norm_pos(1,:),events.starts,[before,after]);
        [posy] = TNC_ExtTrigWins(norm_pos(2,:),events.starts,[before,after]);            
        [spd] = TNC_ExtTrigWins(norm_spd,events.starts,[before,after]);
        [angle] = TNC_ExtTrigWins(norm_angle,events.starts,[before,after]);
        
        % spes per event
        [spe] = TNC_ExtTrigWins(spes_cont,events.starts,[before,after]); 
        [speS] = TNC_ExtTrigWins(spes_cont,events.stops,[before,after]); 
        
        % activity per event
        [ca] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Sig_df,events.starts,[before,after]);

        % compute performance on event basis
        clear *_ln *_id;
        for ii=1:size(amp.wins,1)
            [~,maxi] = max(amp.wins(ii,before:end));
            pthe_ln(1,ii) = trapz(amp.wins(ii,before:before+maxi-1));
        end

        for ii=1:size(ampT.wins,1)
            pth_ln(1,ii) = trapz(ampT.wins(ii,1:before));
        end
        
        if anim_num==3
            blk_id = S.stExp_All.stEvt.Solenoid.RiseFrm>blk_frames;
            blk_id = blk_id(1:size(ampT.wins,1))+1;            
        else
            blk_id = S.stExp_All.stEvt.Collect.RiseFrm>blk_frames;
            blk_id = blk_id(1:size(ampT.wins,1))+1;
        end
        
        for pp=unique(events.blk)
            events.perf.corr(pp) = mean(events.corr(events.blk==pp));
        end
        
        figure(17); clf; 
        subplot(131);
        boxplot(pth_ln,blk_id,'notch','on');
        subplot(132);
%         boxplot(pthe_ln,events.blk,'notch','on');
        subplot(133);
        plot(events.perf.corr); axis([0 3 0 1]);
        
        % train a decoder on the amp/angle data
        smth_spk = zeros(size(S.stExp_All.stSig.Spikes));
        for jj=1:size(S.stExp_All.stSig.Spikes,1)
            smth_spk(jj,:) = conv(S.stExp_All.stSig.Spikes(jj,:),spk_d_kernel,'same');
        end            
        [input_struct] = TNC_ExtTrigWins3d(smth_spk,events.starts,[before,after]);
        [target_struct] = TNC_ExtTrigWins3d([norm_amp;norm_angle],events.starts,[before,after]);
        input_labels.tag = ones(1,size(S.stExp_All.stSig.Spikes,1));
        input_labels.labels = {'dCA1'}; clear decoder;
        folds = 50;
        batch_size = 75; %leave 5 out cross validation for 50 folds
        
            for k=1:folds

                tmp_inds = randperm(size(input_struct.wins,3),batch_size+10);
                pv_data = []; pos_data=[]; pv_out_data=[]; pos_out_data=[];

                % concatenate training data for regression
                for p =1:batch_size
                    pv_data = [pv_data input_struct.wins(:,:,tmp_inds(p))];
                    pos_data = [pos_data target_struct.wins(:,:,tmp_inds(p))];
                end   

                % ignore cells with almost no activity
                pv_data(sum(pv_data,2)<1,:) = 0;
%                 pv_data(sum(pv_data,2)<median(sum(pv_data,2)),:) = 0;
                
                % calculate unregularized regression
                decode_w = pinv(pv_data')*pos_data';

                % parcel out each decoder dimension for later convenience
                for n=1:size(decode_w,2)
                    decoder.decode_w(:,k,n) = decode_w(:,n);
                end
                fprintf('%d ',k);
            end
            
            disp('  ');
            for n=1:size(decode_w,2)
                decoder.cm_w(:,n) = mean(decoder.decode_w(:,:,n),2);
            end

            test_eval_num = 30;
            tmp1 = find(events.blk==1);
            tmp2 = find(events.blk==2);
            tmp_inds = [tmp1(randperm(numel(tmp1),test_eval_num)) tmp2(randperm(numel(tmp2),test_eval_num))];
            pv_data = []; pos_data=[];

            % concatenate training data for regression
            for p = 1:numel(tmp_inds)
                pv_data = [pv_data input_struct.wins(:,:,tmp_inds(p))];
                pos_data = [pos_data target_struct.wins(:,:,tmp_inds(p))];
            end   
            
            pred_pos = pv_data'*decoder.cm_w;
            shuff_pos = pv_data'*decoder.cm_w(randperm(size(decoder.cm_w,1)),:);
            figure(9); clf; 
            subplot(211);
            plot(pred_pos(:,1));
            hold on;
            plot(pos_data(1,:),'Linewidth',2);
            subplot(212);
            plot(pred_pos(:,2));
            hold on;
            plot(pos_data(2,:),'Linewidth',2);
            
            mdl = fitlm(pred_pos(1:size(pred_pos,1)/2,1),pos_data(1,1:size(pred_pos,1)/2)');
            events.decode_corr(1,1) = mdl.Rsquared.Adjusted;
            mdl = fitlm(shuff_pos(1:size(pred_pos,1)/2,1),pos_data(1,1:size(pred_pos,1)/2)');
            events.decode_corr(1,2) =mdl.Rsquared.Adjusted;
            mdl = fitlm(pred_pos(size(pred_pos,1)/2:end,1),pos_data(1,size(pred_pos,1)/2:end)');
            events.decode_corr(2,1) = mdl.Rsquared.Adjusted;
            mdl = fitlm(shuff_pos(size(pred_pos,1)/2:end,1),pos_data(1,size(pred_pos,1)/2:end)');
            events.decode_corr(2,2) = mdl.Rsquared.Adjusted;
            events.decode_corr            
            
        % store all mice data
        allMice.ev(cnt).events = events;
        allMice.anim_num(cnt) = anim_num;
        allMice.blk_frames(cnt) = blk_frames;
        
        allMice.trials(cnt).t = S.stExp_All.stEvt.Trial_S.FallFrm;
        
        allMice.amp(cnt).amp = amp;
        allMice.amp(cnt).ampS = ampS;
        allMice.posx(cnt) = posx;
        allMice.posy(cnt) = posy;
        allMice.posxT(cnt) = posxT;
        allMice.posyT(cnt) = posyT;
        allMice.amp(cnt).before = before;
        allMice.amp(cnt).after = after;
        allMice.amp(cnt).spd = spd;
        allMice.amp(cnt).angle = angle;

        allMice.act(cnt).decoder.cm_w = decoder.cm_w;
        allMice.act(cnt).spe = spe;
        allMice.act(cnt).speS = speS;
        allMice.act(cnt).ca = ca;
        allMice.act(cnt).caT = caT;
        allMice.act(cnt).spT = spT;
        
        allMice.p_tune(cnt) = ranksum(pth_ln(blk_id==1),pth_ln(blk_id==2));
        allMice.p_corr(cnt,:) = events.perf.corr;
        allMice.p_amp(cnt,1) = mean(pth_ln(blk_id==1));
        allMice.p_amp(cnt,2) = mean(pth_ln(blk_id==2));
        
        allMice.pl(cnt).pth_ln = pth_ln;
        allMice.pl(cnt).pthe_ln = pthe_ln;
        allMice.pl(cnt).blk_id = blk_id;
        allMice.iti(cnt,:) = hist(diff(events.starts(events.blk<3))./frames_per_sec*1000,100:100:10000);
        allMice.iti_b(cnt,:) = 100:100:10000;
        
        allMice.data(cnt).spks = S.stExp_All.stSig.Spikes;
        allMice.data(cnt).spe = S.tSyncEvt_V;
        allMice.data(cnt).df = S.stExp_All.stSig.Sig_df;
        allMice.data(cnt).amp = norm_amp;
        allMice.data(cnt).ang = norm_angle;
        allMice.data(cnt).pos = norm_pos;

        % increment counter
        cnt = cnt+1;
        
    end
end

% normalize path lengths across mice (data collection bit depth)
allMice.p_ampN = allMice.p_amp;
for jj=unique(allMice.anim_num)
    
    allMice.p_ampN(allMice.anim_num==jj,:) = allMice.p_amp(allMice.anim_num==jj,:) ./ mean(allMice.p_amp(allMice.anim_num==jj,1));
    
end

crit_sess = find(allMice.p_tune'<0.1 & allMice.p_corr(:,2)>0.5);

disp(['Showing data from n=' num2str(numel(crit_sess)) ' sessions | N=' num2str(numel(unique(allMice.anim_num(crit_sess)))) ' mice']);

figure(32); clf;

% amplitude
subplot(131);
plot([0.5 1.5],[1 1],'k-'); hold on;
plot([1.5 2.5],[1.4 1.4],'k-'); hold on;
boxplot(allMice.p_ampN(crit_sess,:));
axis([0 3 0.75 2.75]); ylabel('Norm. path length');

% p(corr)
subplot(132);
boxplot(allMice.p_corr(crit_sess,1:2));
axis([0 3 0 1]);  ylabel('Prob. corr.');

% timing
subplot(133);
plot(mean(allMice.iti_b(crit_sess,:)),sgolayfilt(sum(allMice.iti(crit_sess,:),1),3,11));
axis([2000 10000 0 100]);

%% Script to load NTF task data from files

cnt = 1; use_amp=0;
for animal = 1:5

    anim_num = animal;
    
    switch animal
        case 1
            dirBase = 'VAO_Ai93/ANM373541/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 17;
        case 2
            dirBase = 'VAO_Ai93/ANM377489/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 17;
        case 3
            dirBase = 'VAO_GP43/ANM261311/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 10;
        case 4
            dirBase = 'VAO_GP43/ANM261313/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 10;
        case 5
            dirBase = 'VAO_GP43/ANM274690/'
            dir_list = dir([dirBase '2*']);
            frames_per_sec = 10;
    end    
    
        % cycle through all data from this mouse
    for jj = 1:numel(dir_list)

        S=load([dirBase dir_list(jj).name '/SyncEvt_Analysis.mat']);
        clear x y norm_pos norm_vel 
        x_min = min(unique(S.stExp_All.vtPos(1,:)));
        x_avg = mean(S.stExp_All.vtPos(1,:));
        y_min = min(unique(S.stExp_All.vtPos(2,:)));
        y_avg = mode(S.stExp_All.vtPos(2,:));

        if anim_num>2 
            norm_pos = [ sgolayfilt(S.stExp_All.vtPos(1,:)-x_avg,3,5) ;  sgolayfilt(S.stExp_All.vtPos(2,:)-y_avg,3,5) ].*2;
        else
            norm_pos = [ sgolayfilt(S.stExp_All.vtPos(1,:)-x_avg,3,5) ;  sgolayfilt(S.stExp_All.vtPos(2,:)-y_avg,3,5) ];
        end
        [norm_angle,norm_amp] = cart2pol(sgolayfilt(S.stExp_All.vtPos(1,:)-x_avg,3,5),sgolayfilt(S.stExp_All.vtPos(2,:)-y_avg,3,5));
        norm_spd = [0 abs(diff(norm_amp))];
        figure(1); clf; plot(norm_amp); hold on; plot(norm_pos(2,:));  plot(norm_pos(1,:));
        
        spk_kernel = TNC_CreateGaussian(50,1,100,1);

        S.tSyncEvt_V = S.tSyncEvt_V;
        spes_cont = zeros(size(norm_amp));
        spes_cont(S.tSyncEvt_V) = 1;
        spes_cont = conv(spes_cont,spk_kernel,'same');

        if numel(S.stExp_All.clBlocks)==4
            thresh_joy = [50 60 55 50];
            blk_lup = [1 2 3 4];
        else
            thresh_joy = [50 50 60 60 55 55 50 50];
            blk_lup = [1 1 2 2 3 3 4 4];
        end            

        % Threshold crosses
        thr_x = find(S.stExp_All.vtState==3 & [0 diff(S.stExp_All.vtState)]==1 & [1:numel(S.stExp_All.vtState)]>2*frames_per_sec);

        min_duration = frames_per_sec*2;
        if use_amp
                detect_thresh = 20;
                [events] = TNC_ExtractPeaks(norm_amp,detect_thresh,min_duration,1);
        else           
            if anim_num>2
                detect_thresh = 2;
                [events] = TNC_ExtractPeaks(sgolayfilt(medfilt1(norm_spd),3,15),detect_thresh,min_duration,1);
            else
                detect_thresh = 5;
                [events] = TNC_ExtractPeaks(sgolayfilt(medfilt1(norm_spd),3,21),detect_thresh,min_duration,1);
            end
        end
    
        events.corr = zeros(1,events.num);
        events.trig = zeros(1,events.num);
        events.blk = ones(1,events.num);
        blk_frames = cumsum(S.stExp_All.vtSessionFrmCount);

        for kk=1:events.num
            
            % did the program find a threshold crossing during the movement
            tmp = find( S.stExp_All.vtState(events.starts(kk):events.stops(kk))==3 , 1 , 'first');            
            if numel(tmp)>0
                events.corr(kk) = 1;
                events.trig(kk) = tmp;
            end

            % what block is this trial?
            tmp = find(events.starts(kk)>blk_frames,1,'last');
            if numel(tmp)==0
                events.blk(kk) = 1;
            else
                events.blk(kk) = blk_lup(tmp+1);
            end
            
        end

        % errors per correct
        tmp = find(events.corr==1);
        events.err_per_corr = [tmp(1) diff(tmp)]-1;
        events.err_per_blk = events.blk(events.corr==1);
        
        % move amp per event
        [amp_x] = TNC_ExtTrigWins(norm_pos(1,:),events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        [amp_y] = TNC_ExtTrigWins(norm_pos(2,:),events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        
        % also calc this aligned to trigger
        [amp_xT] = TNC_ExtTrigWins(norm_pos(1,:),thr_x,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        [amp_yT] = TNC_ExtTrigWins(norm_pos(2,:),thr_x,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        
        [amp] = TNC_ExtTrigWins(norm_amp,events.starts,[frames_per_sec frames_per_sec]);
        [ampT] = TNC_ExtTrigWins(norm_amp,thr_x,[frames_per_sec frames_per_sec]);
        [caT] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Spikes,events.starts,[frames_per_sec frames_per_sec]);
           
        % compute performance on event basis
        clear *_ln *_id;
        for ii=1:size(amp.wins,1)
            [~,maxi] = max(amp.wins(ii,frames_per_sec:end));
            pthe_ln(1,ii) = trapz(amp.wins(ii,frames_per_sec:frames_per_sec+maxi-1));
        end

        for pp=unique(events.blk)
            events.perf.corr(pp) = mean(events.corr(events.blk==pp));
            events.perf.pl(pp) = median(pthe_ln(1,events.blk==pp));
        end
        
        
        % compute tuning
       if anim_num<3 
            for ii=1:size(ampT.wins,1)
                pth_ln(1,ii) = trapz(ampT.wins(ii,frames_per_sec-4:frames_per_sec+2));
            end
            blk_id = [];
            for jj=1:numel(S.stExp_All.clLogLags_Frm)
                blk_id = [blk_id jj*ones(1,size(S.stExp_All.clLogLags_Frm{jj},1))];
            end
            blk_id = blk_id(1:size(ampT.wins,1));
       else
            for ii=1:size(ampT.wins,1)
                pth_ln(1,ii) = trapz(ampT.wins(ii,frames_per_sec-3:frames_per_sec+1));
            end
            blk_id = []; 
            for jj=1:numel(S.stExp_All.clLogLags_Frm)
                blk_id = [blk_id blk_lup(jj)*ones(1,size(S.stExp_All.clLogLags_Frm{jj},1))];
            end
            blk_id = blk_id(1:size(ampT.wins,1));
       end
        figure(17); clf; 
        subplot(131);
        boxplot(pth_ln,blk_id,'notch','on');
        subplot(132);
        boxplot(pthe_ln,events.blk,'notch','on');
        subplot(133);
        plot(events.perf.corr); axis([0 5 0 1]);
        
       
       
        figure(10); clf; plot(amp.avg); hold on; plot(ampT.avg);
            
        [spd] = TNC_ExtTrigWins(norm_spd,events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        [angle] = TNC_ExtTrigWins(norm_angle,events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        
        % spes per event
        [spe] = TNC_ExtTrigWins(spes_cont,events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]); 
        
        % activity per event
        [ca] = TNC_ExtTrigWins3d(S.stExp_All.stSig.Spikes,events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        
        % train a decoder on the amp/angle data
        smth_spk = zeros(size(S.stExp_All.stSig.Spikes));        
        spk_d_kernel = TNC_CreateGaussian(50,1*frames_per_sec/17,100,1);
        for jj=1:size(S.stExp_All.stSig.Spikes,1)
            smth_spk(jj,:) = conv(S.stExp_All.stSig.Spikes(jj,:),spk_d_kernel,'same');
        end            
        [input_struct] = TNC_ExtTrigWins3d(smth_spk,events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        [target_struct] = TNC_ExtTrigWins3d([norm_amp;norm_angle],events.starts,[round(1.5*frames_per_sec) 2*frames_per_sec]);
        clear decoder;
        folds = 50;
        batch_size = 60; %leave 5 out cross validation for 50 folds
        
            for k=1:folds

                tmp_inds = randperm(size(input_struct.wins,3),batch_size+10);
                pv_data = []; pos_data=[]; pv_out_data=[]; pos_out_data=[];

                % concatenate training data for regression
                for p =1:batch_size
                    pv_data = [pv_data input_struct.wins(:,:,tmp_inds(p))];
                    pos_data = [pos_data target_struct.wins(:,:,tmp_inds(p))];
                end   

                % calculate unregularized regression
                decode_w = pinv(pv_data')*pos_data';

                % parcel out each decoder dimension for later convenience
                for n=1:size(decode_w,2)
                    decoder.decode_w(:,k,n) = decode_w(:,n);
                end
                fprintf('%d ',k);
            end
            
            disp('  ');
            for n=1:size(decode_w,2)
                decoder.cm_w(:,n) = mean(decoder.decode_w(:,:,n),2);
            end

            test_eval_num = 20;
            tmp1 = find(events.blk==1 | events.blk==4);
            tmp2 = find(events.blk==2 | events.blk==3);
            tmp_inds = [tmp1(randperm(numel(tmp1)-1,test_eval_num)) tmp2(randperm(numel(tmp2(1:end-1))-2,test_eval_num))];
            pv_data = []; pos_data=[];

            % concatenate training data for regression
            for p = 1:numel(tmp_inds)
                pv_data = [pv_data input_struct.wins(:,:,tmp_inds(p))];
                pos_data = [pos_data target_struct.wins(:,:,tmp_inds(p))];
            end   
            
            pred_pos = pv_data'*decoder.cm_w;
            shuff_pos = pv_data'*decoder.cm_w(randperm(size(decoder.cm_w,1)),:);
            figure(9); clf; 
            subplot(211);
            plot(pred_pos(:,1));
            hold on;
            plot(pos_data(1,:),'Linewidth',2);
            plot(shuff_pos(:,1),'--');
            subplot(212);
            plot(pred_pos(:,2));
            hold on;
            plot(pos_data(2,:),'Linewidth',2);
            plot(shuff_pos(:,2),'--');
            
            mdl = fitlm(pred_pos(1:size(pred_pos,1)/2,1),pos_data(1,1:size(pred_pos,1)/2)');
            events.decode_corr(1,1) = mdl.Rsquared.Adjusted;
            mdl = fitlm(shuff_pos(1:size(pred_pos,1)/2,1),pos_data(1,1:size(pred_pos,1)/2)');
            events.decode_corr(1,2) =mdl.Rsquared.Adjusted;
            mdl = fitlm(pred_pos(size(pred_pos,1)/2:end,1),pos_data(1,size(pred_pos,1)/2:end)');
            events.decode_corr(2,1) = mdl.Rsquared.Adjusted;
            mdl = fitlm(shuff_pos(size(pred_pos,1)/2:end,1),pos_data(1,size(pred_pos,1)/2:end)');
            events.decode_corr(2,2) = mdl.Rsquared.Adjusted;
            events.decode_corr

        % store all mice data
        allMice.ev(cnt).events = events;
        allMice.anim_num(cnt) = anim_num;
        
        allMice.trials(cnt).t = S.stExp_All.stEvt.Trial_S.FallFrm;
        
        allMice.p_corr(cnt) = mean(events.corr(events.blk==1 | events.blk==2));
        
        allMice.amp(cnt).amp_x = amp_x;
        allMice.amp(cnt).amp_y = amp_y;
        allMice.amp(cnt).amp_xT = amp_xT;
        allMice.amp(cnt).amp_yT = amp_yT;
        allMice.amp(cnt).amp = amp;
        allMice.amp(cnt).spd = spd;
        allMice.amp(cnt).angle = angle;

        allMice.act(cnt).spe = spe;
        allMice.act(cnt).ca = ca;
        allMice.act(cnt).caT = caT;
        
        allMice.p_tune(cnt) = ranksum(pth_ln(blk_id==1),pth_ln(blk_id==2));
        allMice.pe_tune(cnt) = ranksum(pthe_ln(events.blk==1),pthe_ln(events.blk==2));
        allMice.p_corr(cnt,1:numel(events.perf.corr)) = events.perf.corr;
        allMice.p_amp(cnt,1) = mean(pth_ln(blk_id==1));
        allMice.p_amp(cnt,2) = mean(pth_ln(blk_id==2));
        
        allMice.pl(cnt).pth_ln = pth_ln;
        allMice.pl(cnt).pthe_ln = pthe_ln;
        allMice.pl(cnt).blk_id = blk_id;
        allMice.iti(cnt,:) = hist(diff(events.starts(events.blk<3))./frames_per_sec*1000,100:100:10000);
        allMice.iti_b(cnt,:) = 100:100:10000;
        
        allMice.data(cnt).spks = S.stExp_All.stSig.Spikes;
        allMice.data(cnt).df = S.stExp_All.stSig.Sig_df;        
        allMice.data(cnt).spe = S.tSyncEvt_V;
        allMice.data(cnt).amp = norm_amp;
        allMice.data(cnt).ang = norm_angle;
        allMice.data(cnt).pos = norm_pos;

        % increment counter
        cnt = cnt+1;
        
    end
end

% normalize path lengths across mice (data collection bit depth)
allMice.p_ampN = allMice.p_amp;
for jj=unique(allMice.anim_num)
    
    allMice.p_ampN(allMice.anim_num==jj,:) = allMice.p_amp(allMice.anim_num==jj,:) ./ mean(allMice.p_amp(allMice.anim_num==jj,1));
    
end

% crit_sess = allMice.p_tune<0.05;

crit_sess = find(allMice.p_tune'<0.1 & mean(allMice.p_corr,2)>0.5);

disp(['Showing data from n=' num2str(numel(crit_sess)) ' sessions | N=' num2str(numel(unique(allMice.anim_num(crit_sess)))) ' mice']);

figure(31); clf;

% amplitude
subplot(131);
plot([0.5 1.5],[1 1],'k-'); hold on;
plot([1.5 2.5],[1.2 1.2],'k-'); hold on;
boxplot(allMice.p_ampN(crit_sess,:));
axis([0 3 0.75 1.75]); ylabel('Norm. path length');

% p(corr)
subplot(132);
boxplot(allMice.p_corr(crit_sess,1:2));
axis([0 3 0 1]);  ylabel('Prob. corr.');

% timing
subplot(133);
plot(mean(allMice.iti_b(crit_sess,:)),sgolayfilt(sum(allMice.iti(crit_sess,:),1),3,11));
axis([2000 10000 0 100]);

%% tML Model script example

clear *_pop_avg *_pop_all mu* om* targe* fit_grid_search
cnt=0;
intercept = 'anypt'

% rk = TNC_CreateGaussian(50,0.5,100,1);
% rk = [0 0.5 0.5 0];
rk = [0 1 0];

% Optimum of refined at cnt=12 or a=0.8 s=20
alphas = [ 0.1 0.2 0.4 0.6 0.8 ];
sig_speed = [ 10 15 20 25 30 50 ];
sig_hd = 0.25;

for alpha = alphas
    for sig_spd = sig_speed
        cnt = cnt+1

                for qqq=1:10

                clear pos traj
                arena.dims = [-300 300 0 600];
                blocks = [80 60];

                arena.x = arena.dims([1 2 2 1]);
                arena.y = arena.dims([3 3 4 4]);
                smth_win = 20;
                target.cntr = [0 350];
%                 target.width = [50 50]; % STF task
                target.width = [300 50]; % NTF task
                target.x = [target.cntr(1)-target.width(1) target.cntr(1)+target.width(1) target.cntr(1)+target.width(1) target.cntr(1)-target.width(1)];
                target.y = [target.cntr(2)-target.width(2) target.cntr(2)-target.width(2) target.cntr(2)+target.width(2) target.cntr(2)+target.width(2)];

                target2.cntr = [0 430];
%                 target2.width = [50 50]; % STF task
                target2.width = [300 50]; % NTF task
                target2.x = [target2.cntr(1)-target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)+target2.width(1) target2.cntr(1)-target2.width(1)];
                target2.y = [target2.cntr(2)-target2.width(2) target2.cntr(2)-target2.width(2) target2.cntr(2)+target2.width(2) target2.cntr(2)+target2.width(2)];

                % Noise models for the dilation parameter
                tmp = -0.5 + (conv(rand(1,sum(blocks)),rk,'same')).^2;
                tmp(1) = 0;
                clear tmp2;
                for jj=1:sum(blocks)
%                     tmp2(jj) = randperm(20,1)+1;
                    tmp2(jj) = randperm(20,1)+1;
                end
                tmp_om = 0.01*randn(1,sum(blocks));
                tmp_om(1) = 0;
                
                % draw arena
                figure(10); clf;
                        subplot(131); 
                        patch(arena.x,arena.y,[0.95 0.9 0.95]); hold on;
                        patch(target.x,target.y,[0.85 0.9 0.9]);
                        patch(target2.x,target2.y,[0.9 0.9 0.85]);

                mu_spd(1) = 230+randperm(80,1); 
                om_hd(1) = 0;
                beta=0.01;
                target_hit = zeros(1,sum(blocks)); 
                hit_time = ones(1,sum(blocks)); rew_rate=zeros(1,100); 


                % run a trial
                for k=1:sum(blocks) %:100

                    % estimate reward rate
                    if k<=3
                         rew_rate(k) = 0 ;
                    else
                         rew_rate(k) = (1 - sum([target_hit(k-3).*0.1 target_hit(k-2).*0.3 target_hit(k-1).*0.6]));    
                    end

                    % calc speed profile
                    speed = [TNC_CreateGaussian(250,55+randn(1),500,1) TNC_CreateGaussian(250,55+randn(1),500,1)];                    
                    speed = ( mu_spd(k) + (tmp(k) .* sig_spd) + (rew_rate(k) .* tmp2(k)) ) .* speed;
                    obs_spd(k) = mu_spd(k) + (tmp(k) .* sig_spd) + (rew_rate(k) .* tmp2(k) .* sig_spd./10) ;

                    % calc heading profile
                    heading = [0.001:2*pi/1000:2*pi] + sgolayfilt( conv( randn(size(speed)) , [0 ones(1,smth_win)./smth_win 0] , 'same' ) , 3 , 201 );
                    if rand(1)<0.5
                        heading = 2*pi - heading;
                    end
                    this_hd = om_hd(k) + (rew_rate(k) .* tmp_om(k) .* sig_hd);
                    heading = heading + this_hd;

                    % initialize some variables for sim
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

                        switch intercept
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
                    
                    if k>1 & target_hit(k)==1
                        
                    % NOTE: need to switch heading update learning to this:
                    % calculate angle the way it is calculated in data
                    [obs_hd(k),~] = cart2pol(pos.y(150),pos.x(150));
                            
                        d_mu = alpha .* (obs_spd(k) - mu_spd(k));
                        d_mu_h = ( beta .* (mu_spd(k) - mu_spd(1)) );
                        
                        d_om = alpha .* (this_hd - om_hd(k));
                        d_om_h = ( beta .* om_hd(k)-om_hd(1) );
                        
                        mu_spd(k+1) = mu_spd(k) + d_mu - d_mu_h;
                        om_hd(k+1) = om_hd(k) + d_om - d_om_h;


                    else % how to explore given that no reward was obtained

                        d_mu = alpha .* (obs_spd(k) - mu_spd(k));
                        mu_spd(k+1) = mu_spd(k) + d_mu;
                        om_hd(k+1) = om_hd(k);
                        if mu_spd(k+1) > 500
                            mu_spd(k+1) = 500;
                        end
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
%                         subplot(1,3,2); hold off; plot(mu_spd(1:k),'o-','color',[0.5 0.5 0.5]); hold on; plot(100*target_hit(1:k)+275,'*'); scatter(1:k,mu_spd(1:k),50,[1:k]>blocks(1),'filled'); colormap([ 1 0 0 ; 0 0 1 ]);
                        subplot(1,3,2); hold off; plot(mu_spd(1:k),'o-','color',[0.5 0.5 0.5]); hold on; plot(10*target_hit(1:k)+275,'*'); scatter(1:k,obs_spd(1:k),50,[1:k]>blocks(1),'filled'); colormap([ 1 0 0 ; 0 0 1 ]);

                    if k==sum(blocks)
                        legend('Speed (mean)','Success');
                        box off;
                    end

                    subplot(133); hold off; plot(om_hd(1:k),'o-','color',[0.5 0.5 0.5]); hold on; plot(target_hit(1:k),'*'); plot(rew_rate(1:k),'-');
                    if qqq==1
                        drawnow;
                    end
                end

                transition(qqq) = blocks(1)+1;
                mu_pop_avg(qqq,:) = mu_spd;
                obs_pop_avg(qqq,:) = obs_spd;
                corr_pop_avg(qqq,:) = target_hit;

                end
                
                
scale = 1;                
                                
figure(302); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
subplot(121);
% shadedErrorBar(1:blocks(1),mean(mu_pop_avg(:,1:blocks(1))),std(mu_pop_avg(:,1:blocks(1)))./sqrt(size(mu_pop_avg,1))./scale,{'color',cat_map(1,:)}); hold on;
% shadedErrorBar(blocks(1)+1:size(mu_pop_avg,2),mean(mu_pop_avg(:,blocks(1)+1:size(mu_pop_avg,2))),std(mu_pop_avg(:,blocks(1)+1:size(mu_pop_avg,2)))./sqrt(size(mu_pop_avg,1))./scale,{'color',cat_map(2,:)}); hold on;
shadedErrorBar(1:blocks(1),mean(obs_pop_avg(:,1:blocks(1))),std(obs_pop_avg(:,1:blocks(1)))./sqrt(size(obs_pop_avg,1))./scale,{'color',cat_map(1,:)}); hold on;
shadedErrorBar(blocks(1)+1:size(obs_pop_avg,2),mean(obs_pop_avg(:,blocks(1)+1:size(obs_pop_avg,2))),std(obs_pop_avg(:,blocks(1)+1:size(obs_pop_avg,2)))./sqrt(size(obs_pop_avg,1))./scale,{'color',cat_map(2,:)}); hold on;
ylabel('Scaled amplitude');
xlabel('Trials');

subplot(122);
% shadedErrorBar(1:size(corr_pop_avg,2),1-mean(corr_pop_avg),std(corr_pop_avg)./sqrt(size(corr_pop_avg,1)),{'color',[0 0.33 0.67]}); hold on;
shadedErrorBar(1:blocks(1),mean(corr_pop_avg(:,1:blocks(1))),std(corr_pop_avg(:,1:blocks(1)))./sqrt(size(corr_pop_avg,1))./scale,{'color',cat_map(1,:)}); hold on;
shadedErrorBar(blocks(1)+1:size(corr_pop_avg,2),mean(corr_pop_avg(:,blocks(1)+1:size(corr_pop_avg,2))),std(corr_pop_avg(:,blocks(1)+1:size(corr_pop_avg,2)))./sqrt(size(corr_pop_avg,1))./scale,{'color',cat_map(2,:)}); hold on;
ylabel('Errors per trial');
xlabel('Trials');

figure(303); clf;
[cat_map] = TNC_CreateRBColormap(8,'mapb');
axis tight;
subplot(121);
for jj=1:blocks(1)
    plot(traj.x(jj,:),traj.y(jj,:),'color',[cat_map(1,:) 0.1]); hold on;
end
% plot(mean(traj.x(1:blocks(1),:)),mean(traj.y(1:blocks(1),:)),'color',cat_map(1,:),'linewidth',2);
ylabel('Y');
xlabel('X');
axis([-300 300 0 600]);


subplot(122);
for jj=blocks(1)+1:size(traj.x,1)
    plot(traj.x(jj,:),traj.y(jj,:),'color',[cat_map(2,:) 0.1]); hold on;
end
% plot(mean(traj.x(blocks(1)+1:end,:)),mean(traj.y(blocks(1)+1:end,:)),'color',cat_map(2,:),'linewidth',2);
ylabel('Y');
xlabel('X');
axis([-300 300 0 600]);

        clear out_mean;
        for zz=1:size(obs_pop_avg,1)
            out_mean(zz,:) = (scaler+obs_pop_avg(zz,blocks(1)-14:blocks(1)+56)) ./ (scaler+mean(obs_pop_avg(zz,blocks(1)-14:blocks(1))));
        end

        fit_grid_search(cnt).alpha      = alpha;
        fit_grid_search(cnt).sig_spd    = sig_spd;
        fit_grid_search(cnt).algn_avg   = mean(out_mean,1);
        fit_grid_search(cnt).algn_err   = std(out_mean,[],1);        
        fit_grid_search(cnt).rmse_stf   = sqrt( mean( (fit_grid_search(cnt).algn_avg-mean(to_fit.stf.aligned_pth_ln)).^2 ) );
        fit_grid_search(cnt).rmse_ntf   = sqrt( mean( (fit_grid_search(cnt).algn_avg-mean(to_fit.ntf.aligned_pth_ln)).^2 ) );
        
        figure(305); clf; % compare to data
        subplot(211);
        shadedErrorBar(-before:after,fit_grid_search(cnt).algn_avg,fit_grid_search(cnt).algn_err,{'color',[1 0.67 0]});
         hold on;
        plot(-before:after,mean(to_fit.stf.aligned_pth_ln),'color',[0.5 0.5 0.5]); hold on;

        subplot(212);
        shadedErrorBar(-before:after,fit_grid_search(cnt).algn_avg,fit_grid_search(cnt).algn_err,{'color',[1 0.67 0]});
         hold on;
         plot(-before:after,mean(to_fit.ntf.aligned_pth_ln),'color',[0.5 0.5 0.5]);

    end
end


%% After running 
clear all_rmse;
[mse] = TNC_CreateRBColormap(8,'yb');
for gg=1:numel(fit_grid_search)
    all_rmse_stf(find(fit_grid_search(gg).alpha==alphas),find(fit_grid_search(gg).sig_spd==sig_speed)) = fit_grid_search(gg).rmse_stf;
    all_rmse_ntf(find(fit_grid_search(gg).alpha==alphas),find(fit_grid_search(gg).sig_spd==sig_speed)) = fit_grid_search(gg).rmse_ntf;
    indices(find(fit_grid_search(gg).alpha==alphas),find(fit_grid_search(gg).sig_spd==sig_speed)) = gg;
end
figure(50); 
subplot(211); imagesc(alpha,sig_speed,all_rmse_stf); colormap(mse);
subplot(212); imagesc(alpha,sig_speed,all_rmse_ntf); colormap(mse);
[i,j] = find(all_rmse_stf==min(min(all_rmse_stf)))
[iN,jN] = find(all_rmse_ntf==min(min(all_rmse_ntf)))



%         figure(305); clf; % compare to data
%         subplot(211);
%         shadedErrorBar(-before:after,fit_grid_search(indices(i,j)).algn_avg,fit_grid_search(indices(i,j)).algn_err,{'color',[235 0 139]/255});
%          hold on;
%         plot(-before:after,mean(to_fit.stf.aligned_pth_ln),'color',[0.5 0.5 0.5]); hold on; box off;
% 
%         subplot(212);
%         shadedErrorBar(-before:after,fit_grid_search(indices(iN,jN)).algn_avg,fit_grid_search(indices(iN,jN)).algn_err,{'color',[235 0 139]/255});
%          hold on;
%          plot(-before:after,mean(to_fit.ntf.aligned_pth_ln),'color',[0.5 0.5 0.5]); box off;

        figure(305); clf; % compare to data
        subplot(211);
        plot([0 0],[0.8 1.35],'k-');
         hold on;  plot([-before after],[1 1],'k-');
        shadedErrorBar(-before:after,mean(to_fit.stf.aligned_pth_ln),std(to_fit.stf.aligned_pth_ln)./sqrt(size(to_fit.stf.aligned_pth_ln,1)),{'color',[0 0 0]}); hold on; box off;
        plot(-before:after,fit_grid_search(indices(i,j)).algn_avg,'color',[235 0 139]/255);
        

        subplot(212);
        plot([0 0],[0.8 1.35],'k-');
         hold on; plot([-before after],[1 1],'k-');
         shadedErrorBar(-before:after,mean(to_fit.ntf.aligned_pth_ln),std(to_fit.ntf.aligned_pth_ln)./sqrt(size(to_fit.ntf.aligned_pth_ln,1)),{'color',[0 0 0]}); box off;
        plot(-before:after,fit_grid_search(indices(iN,jN)).algn_avg,'color',[235 0 139]/255);
