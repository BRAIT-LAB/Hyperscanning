clc;clear all

% path ='E:\BBI\EEGData\Form\exp1_picture\ppl\*.set';
% addpath(path);
savename = {'resting','pol','poh','ppl','pph','no7','in7','now','inw','mol','moh','mpl','mph'};
for h = 1:13
    mkdir(['E:\analysis\results\epoch\ISC_all_mat\' savename{h} ]);
end

for k = 1:13
    path = ['E:\analysis\results\epoch\' savename{k} '\*.set'];
    namelist = dir(path);
    len = length(namelist);
    for i = 1:len
        tempfilepath = [ path(1:end-5) namelist(i).name]
        filename{i}= tempfilepath;
    end
    folder = dir(fullfile(path));
    folder = {folder.name};
    channel = {'Fp1','F3','F7','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz',...
        'P3','P7','O1','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FC6',...
        'FC2','F4','F8','Fp2','Fz','Oz','FT9','FT10'};
    locfile = 'E:\BBI\processing\ISC code\ele\BioSemi32.loc';
    comp_no = (1:1:32);
    Fbp1 = 8; Fbp2 = 13;   % 带通滤波
    eloc = readlocs('E:\BBI\EEGData\bbi.locs');
    eloc_label = {eloc.labels};
    % Fbp1 = 8; Fbp2 = 13;   % 带通滤波
    %% 独白基线
    for i = 1:2:40
        
        % 男生基线
        EEG_M = pop_loadset(filename{i});
        EEG_M = pop_resample(EEG_M,250);
        M_label = {EEG_M.chanlocs.labels};
        [C,IA,IB] = intersect(M_label,eloc_label,'stable');
        EEG_M = pop_select(EEG_M,'channel',M_label(IA));
        EEG_M.chanlocs = eloc(IB);
        missingchan = setdiff([1:32],IB);
        for ch = 1:length(missingchan)
            EEG_M.chanlocs = cat(2,EEG_M.chanlocs,eloc(missingchan(ch)));
        end
        M_label_new = {EEG_M.chanlocs.labels};
        [~,badchans] = setdiff(M_label_new,M_label);
        EEG_M = pop_interp(EEG_M,badchans);
        if strcmp('Fz',M_label_new(badchans))
            [~,Fz_idx,IC] = intersect(M_label_new,{'Fz'});
            EEG_M.data(Fz_idx,:)= 0;
        end
        % 去除50Hz工频干扰
        EEG_M = pop_eegfilt(EEG_M,49,51,[],1,1);
        %             % 0-75Hz滤波
        %             EEG_M = pop_eegfilt(EEG_M,0,75,[],[],1);
        % 重参考
        EEG_M = pop_reref(EEG_M,[]);
        % ICA 去除伪迹
        EEG_M = pop_runica(EEG_M,'icatype','runica');
        [kurt,rej_no] = rejkurt(EEG_M.icaactivations,1.64,[],1);
        rej_comp = comp_no(rej_no);
        EEG_M = pop_subcomp(EEG_M,rej_comp);
        % 带通滤波
        EEG_M_delta = pop_eegfilt(EEG_M,1,3,[],[],1);
        EEG_M_theta = pop_eegfilt(EEG_M,4,7,[],[],1);
        EEG_M_alpha = pop_eegfilt(EEG_M,8,12,[],[],1);
        EEG_M_beta = pop_eegfilt(EEG_M,13,30,[],[],1);
        EEG_M_gamma = pop_eegfilt(EEG_M,31,49,[],[],1);
        
        % 女生基线
        EEG_W = pop_loadset(filename{i+1});
        W_label = {EEG_W.chanlocs.labels};
        [C,IA,IB] = intersect(W_label,eloc_label,'stable');
        EEG_W = pop_select(EEG_W,'channel',W_label(IA));
        EEG_W.chanlocs = eloc(IB);
        missingchan = setdiff([1:32],IB);
        for ch = 1:length(missingchan)
            EEG_W.chanlocs = cat(2,EEG_W.chanlocs,eloc(missingchan(ch)));
        end
        W_label_new = {EEG_W.chanlocs.labels};
        [~,badchans] = setdiff(W_label_new,W_label);
        EEG_W = pop_interp(EEG_W,badchans);
        if strcmp('Fz',W_label_new(badchans))
            [~,Fz_idx,IC] = intersect(W_label_new,{'Fz'});
            EEG_W.data(Fz_idx,:)= 0;
        end
        % 去除50Hz工频干扰
        EEG_W = pop_eegfilt(EEG_W,49,51,[],1,1);
        % 重参考
        EEG_W = pop_reref(EEG_W,[]);
        %             % 0-75Hz滤波
        %             EEG_W = pop_eegfilt(EEG_W,0,75,[],[],1);
        % ICA 去除伪迹
        EEG_W = pop_runica(EEG_W,'icatype','runica');
        [kurt,rej_no] = rejkurt(EEG_W.icaactivations,1.64,[],1);
        rej_comp = comp_no(rej_no);
        EEG_W = pop_subcomp(EEG_W,rej_comp);
        % 带通滤波
        EEG_W_delta = pop_eegfilt(EEG_W,1,3,[],[],1);
        EEG_W_theta = pop_eegfilt(EEG_W,4,7,[],[],1);
        EEG_W_alpha = pop_eegfilt(EEG_W,8,12,[],[],1);
        EEG_W_beta = pop_eegfilt(EEG_W,13,30,[],[],1);
        EEG_W_gamma = pop_eegfilt(EEG_W,31,49,[],[],1);
        
        Locs_M = {EEG_M.chanlocs.labels};
        Locs_W = {EEG_W.chanlocs.labels};
        [~,idx_M,idx_W] = intersect(Locs_M,Locs_W,'stable');
        DATA_M = double(EEG_M.data');
        DATA_W = double(EEG_W.data');
        DATA_M = DATA_M(:,idx_M);
        DATA_W = DATA_W(:,idx_W);
        clear EEG_M EEG_W
        %% alpha频段
        DATA_M_alpha = double(EEG_M_alpha.data');
        DATA_W_alpha = double(EEG_W_alpha.data');
        DATA_M_alpha = DATA_M_alpha(:,idx_M);
        DATA_W_alpha = DATA_W_alpha(:,idx_W);
        s1 = size(DATA_M_alpha,1);
        s2 = size(DATA_W_alpha,1);
        min_s = min(s1,s2);
        DATA_M_alpha = DATA_M_alpha(1:min_s,:);
        DATA_W_alpha = DATA_W_alpha(1:min_s,:);
        DATA_alpha = cat(3,DATA_M_alpha,DATA_W_alpha);
        EEGData_alpha.X = DATA_alpha;
        EEGData_alpha.fs = 250;
        EEGData_alpha.badchannels = {};
        EEGData_alpha.eogchannels = [];
        save('E:\BBI\processing\ISC code\EEGData_alpha','-struct','EEGData_alpha')
        clear EEG_M_alpha DATA_M_alpha DATA_W_alpha DATA_alpha EEGData_alpha
        
        datafile = 'E:\BBI\processing\ISC code\EEGData_alpha.mat';
        gamma = 0.4;
        Nsec = 1;
        plotfig = 0;
        [ISC_alpha,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
        ISC_alpha_all(i,:) = ISC_alpha(1);
        %             ISC_persecond_alpha_all(i-2,j,m-2,:) = ISC_persecond_alpha(1,:);
        %             A_alpha_all(i-2,j,m-2,:) = A_alpha(1,:);
        
        %% beta频段
        DATA_M_beta = double(EEG_M_beta.data');
        DATA_W_beta = double(EEG_W_beta.data');
        DATA_M_beta = DATA_M_beta(:,idx_M);
        DATA_W_beta = DATA_W_beta(:,idx_W);
        s1 = size(DATA_M_beta,1);
        s2 = size(DATA_W_beta,1);
        min_s = min(s1,s2);
        DATA_M_beta = DATA_M_beta(1:min_s,:);
        DATA_W_beta = DATA_W_beta(1:min_s,:);
        
        DATA_beta = cat(3,DATA_M_beta,DATA_W_beta);
        EEGData_beta.X = DATA_beta;
        EEGData_beta.fs = 250;
        EEGData_beta.badchannels = {};
        EEGData_beta.eogchannels = [];
        save('E:\BBI\processing\ISC code\EEGData_beta','-struct','EEGData_beta')
        clear EEG_M_beta DATA_M_beta DATA_W_beta DATA_beta EEGData_beta
        
        datafile = 'E:\BBI\processing\ISC code\EEGData_beta.mat';
        gamma = 0.4;
        Nsec = 1;
        plotfig = 0;
        [ISC_beta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
        ISC_beta_all(i,:) = ISC_beta(1);
        %             ISC_persecond_beta_all(i-2,j,m-2,:) = ISC_persecond_beta(1,:);
        %             A_beta_all(i-2,j,m-2,:) = A_beta(1,:);
        
        %% delta频段
        DATA_M_delta = double(EEG_M_delta.data');
        DATA_W_delta = double(EEG_W_delta.data');
        DATA_M_delta = DATA_M_delta(:,idx_M);
        DATA_W_delta = DATA_W_delta(:,idx_W);
        s1 = size(DATA_M_delta,1);
        s2 = size(DATA_W_delta,1);
        min_s = min(s1,s2);
        DATA_M_delta = DATA_M_delta(1:min_s,:);
        DATA_W_delta = DATA_W_delta(1:min_s,:);
        
        DATA_delta = cat(3,DATA_M_delta,DATA_W_delta);
        EEGData_delta.X = DATA_delta;
        EEGData_delta.fs = 250;
        EEGData_delta.badchannels = {};
        EEGData_delta.eogchannels = [];
        save('E:\BBI\processing\ISC code\EEGData_delta','-struct','EEGData_delta')
        clear EEG_M_delta DATA_M_delta DATA_W_delta DATA_delta EEGData_delta
        
        datafile = 'E:\BBI\processing\ISC code\EEGData_delta.mat';
        gamma = 0.4;
        Nsec = 1;
        plotfig = 0;
        [ISC_delta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
        ISC_delta_all(i,:) = ISC_delta(1);
        %             ISC_persecond_delta_all(i-2,j,m-2,:) = ISC_persecond_delta(1,:);
        %             A_delta_all(i-2,j,m-2,:) = A_delta(1,:);
        
        %% gamma频段
        DATA_M_gamma = double(EEG_M_gamma.data');
        DATA_W_gamma = double(EEG_W_gamma.data');
        DATA_M_gamma = DATA_M_gamma(:,idx_M);
        DATA_W_gamma = DATA_W_gamma(:,idx_W);
        s1 = size(DATA_M_gamma,1);
        s2 = size(DATA_W_gamma,1);
        min_s = min(s1,s2);
        DATA_M_gamma = DATA_M_gamma(1:min_s,:);
        DATA_W_gamma = DATA_W_gamma(1:min_s,:);
        
        DATA_gamma = cat(3,DATA_M_gamma,DATA_W_gamma);
        EEGData_gamma.X = DATA_gamma;
        EEGData_gamma.fs = 250;
        EEGData_gamma.badchannels = {};
        EEGData_gamma.eogchannels = [];
        save('E:\BBI\processing\ISC code\EEGData_gamma','-struct','EEGData_gamma')
        clear EEG_M_gamma DATA_M_gamma DATA_W_gamma DATA_gamma EEGData_gamma
        
        datafile = 'E:\BBI\processing\ISC code\EEGData_gamma.mat';
        gamma = 0.4;
        Nsec = 1;
        plotfig = 0;
        [ISC_gamma,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
        ISC_gamma_all(i,:) = ISC_gamma(1);
        %             ISC_persecond_gamma_all(i-2,j,m-2,:) = ISC_persecond_gamma(1,:);
        %             A_gamma_all(i-2,j,m-2,:) = A_gamma(1,:);
        
        %% theta频段
        DATA_M_theta = double(EEG_M_theta.data');
        DATA_W_theta = double(EEG_W_theta.data');
        DATA_M_theta = DATA_M_theta(:,idx_M);
        DATA_W_theta = DATA_W_theta(:,idx_W);
        s1 = size(DATA_M_theta,1);
        s2 = size(DATA_W_theta,1);
        min_s = min(s1,s2);
        DATA_M_theta = DATA_M_theta(1:min_s,:);
        DATA_W_theta = DATA_W_theta(1:min_s,:);
        
        DATA_theta = cat(3,DATA_M_theta,DATA_W_theta);
        EEGData_theta.X = DATA_theta;
        EEGData_theta.fs = 250;
        EEGData_theta.badchannels = {};
        EEGData_theta.eogchannels = [];
        save('E:\BBI\processing\ISC code\EEGData_theta','-struct','EEGData_theta')
        clear EEG_M_theta DATA_M_theta DATA_W_theta DATA_theta EEGData_theta
        
        datafile = 'E:\BBI\processing\ISC code\EEGData_theta.mat';
        gamma = 0.4;
        Nsec = 1;
        plotfig = 0;
        [ISC_theta,ISC_persubject,ISC_persecond,W,A] = isceeg(datafile,locfile,gamma,Nsec,plotfig);
        ISC_theta_all(i,:) = ISC_theta(1);
        %             ISC_persecond_theta_all(i-2,j,m-2,:) = ISC_persecond_theta(1,:);
        %             A_theta_all(i-2,j,m-2,:) = A_theta(1,:);
        
        
    end
    

    
    save(['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\ISC_alpha_all.mat'],'ISC_alpha_all');
    save(['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\ISC_beta_all.mat'],'ISC_beta_all');
    save(['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\ISC_delta_all.mat'],'ISC_delta_all');
    save(['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\ISC_gamma_all.mat'],'ISC_gamma_all');
    save(['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\ISC_theta_all.mat'],'ISC_theta_all');
    
    file_name = {'ISC_alpha_all','ISC_beta_all','ISC_delta_all','ISC_gamma_all','ISC_theta_all'};
    for i= 1:length(file_name)
        file_path= ['E:\analysis\results\epoch\ISC_all_mat\' savename{k} '\' file_name{i}  '.mat'];
        load(file_path)
        data=[];
        for j=1:40
            eval( ['data(', num2str(j), [',:)=' file_name{i} '(' num2str(j) ',1);']] );
        end
        %     writetable(T,['E:\analysis\results\epoch\resting' file_name{i} '.xlsx'],'WriteRowNames',true)
        
        % 再把data写到excel
        xlswrite(['E:\analysis\results\epoch\results\' savename{k} '\' file_name{i} '.xlsx'],data);
    end
    
end



