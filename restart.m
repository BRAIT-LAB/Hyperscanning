clear;
clc;
savename = {'resting','pol','poh','ppl','pph','no7','in7','now','inw','mol','moh','mpl','mph'};
for h =1:1
    disp(["------------->当前是：" savename{h} "<-----------------"]);
    %figure(h+2);
    path = ['E:\analysis\results\epoch\' savename{h} '\*.set'];
    namelist = dir(path);%dir()函数用于获得指定文件夹下的所有子文件夹和文件，并存放在一种文件结构体数组中。
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
    locfile = 'E:\analysis\code-2019-01-19 (2)\code-2019-01-19\BioSemi32.loc';
    comp_no = (1:1:32);
    Fbp1 = 8; Fbp2 = 13;   % 带通滤波
    eloc = readlocs('E:\BBI\EEGData\bbi.locs');
    eloc_label = {eloc.labels};
    for i = 1:2:1 % # of forms

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


        EEG_M = pop_eegfilt(EEG_M,Fbp1,Fbp2,[],[],1);


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

        EEG_W = pop_eegfilt(EEG_W,Fbp1,Fbp2,[],[],1);

        %pop_eegplot( EEG_M, 1, 1, 1);
        %pop_eegplot(EEG_W,1,1,2);

%         %PSD 分析
%         Fs = 250;
%         Fm = 50;
%         deltaf =0.5;
% 
%         [Svv_M,F_M,Ns_M,PSD_M] = xspectrum(EEG_M.data,Fs,Fm,deltaf);%PSD=输入脑电图数据的估计功率谱密度
%         [Svv_W,F_W,Ns_W,PSD_W] = xspectrum(EEG_W.start,Fs,Fm,deltaf);
% 
%         PSD_M = mean(PSD_M,1);
%         PSD_W = mean(PSD_W,1);
% 
%         subplot(5,5,(i+1)/2);
%         plot(F_M,10*log10(PSD_M),F_W,10*log10(PSD_W));
%         title([savename{h}]);
%         xlabel('F/Hz');ylabel('PSD/dB');
%         xlim([0 50]);
%         legend('M','W');




         %PLV
%          meth = 'Trad';
%          fs = 250;
%          f0 = 20;
%          point = min(length(EEG_M.data),length(EEG_W.data));
%          for ch=1:32
%              result(h).PLV((i+1)/2,ch) = PLV_RawSig(meth, EEG_M.data(ch,1:point), fs, f0, [], [], [], EEG_W.data(ch,1:point));
%          end
              
    end
    
%          save(['E:\analysis\results\epoch\PLV\' savename{2} ] ,'-struct', 'result_2');
%          boxplot(result(2).PLV,'Labels',{'Fp1','F3','F7','FC5','FC1','C3','T7','TP9','CP5','CP1','Pz','P3','P7','O1','O2','P4','P8','TP10','CP6','CP2','Cz','C4','T8','FC6', 'FC2','F4','F8','Fp2','Fz','Oz','FT9','FT10'});
% xlabel('Channel');ylabel('PLV');
% title('Phase locked value of EEG data in pol state');
%  save(['E:\analysis\results\epoch\PLV\' savename{h}'],'-struct','result(3).PLV');  
%               for i = 1:2:19
%                   for ch = 1:32
%               result_2(1).PLV((i+1)/2,ch) = result(3).PLV((i+1)/2,ch)
%                   end
%               end
         
end




