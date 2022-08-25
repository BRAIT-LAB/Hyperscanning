% clear all;
% clc
% path = 'E:\analysis\MtoW\exp4_music\mpl';
% addpath(path);
savename = {'resting','pol','poh','ppl','pph','no7','in7','now','inw','mol','moh','mpl','mph'};
% for h = 1:13
%     mkdir(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} ]);
% end
for h =4:4
    path = ['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\' savename{h} '\*.set'];
    namelist = dir(path);%dir()�������ڻ��ָ���ļ����µ��������ļ��к��ļ����������һ���ļ��ṹ�������С�
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
    Fbp1 = 8; Fbp2 = 13;   % ��ͨ�˲�
    eloc = readlocs('E:\BBI\EEGData\bbi.locs');
    eloc_label = {eloc.labels};
    %% ���׻���
    count = 1;
    for i = 1:2:68 % # of forms

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
                % ȥ��50Hz��Ƶ����
                EEG_M = pop_eegfilt(EEG_M,49,51,[],1,1);
    %             % 0-75Hz�˲�
    %             EEG_M = pop_eegfilt(EEG_M,0,75,[],[],1);
                % �زο�
                EEG_M = pop_reref(EEG_M,[]);
                % ICA ȥ��α��
                EEG_M = pop_runica(EEG_M,'icatype','runica');
                [kurt,rej_no] = rejkurt(EEG_M.icaactivations,1.64,[],1);
                rej_comp = comp_no(rej_no);
                EEG_M = pop_subcomp(EEG_M,rej_comp);


                EEG_M = pop_eegfilt(EEG_M,Fbp1,Fbp2,[],[],1);
                % Ů������
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
                % ȥ��50Hz��Ƶ����
                EEG_W = pop_eegfilt(EEG_W,49,51,[],1,1);
                % �زο�
                EEG_W = pop_reref(EEG_W,[]);
    %             % 0-75Hz�˲�
    %             EEG_W = pop_eegfilt(EEG_W,0,75,[],[],1);
                % ICA ȥ��α��
                EEG_W = pop_runica(EEG_W,'icatype','runica');
                [kurt,rej_no] = rejkurt(EEG_W.icaactivations,1.64,[],1);
                rej_comp = comp_no(rej_no);
                EEG_W = pop_subcomp(EEG_W,rej_comp);
                % ��ͨ�˲�
    %             EEG_W_delta = pop_eegfilt(EEG_W,1,3,[],[],1);
    %             EEG_W_theta = pop_eegfilt(EEG_W,4,7,[],[],1);
    %             EEG_W_alpha = pop_eegfilt(EEG_W,8,12,[],[],1);
    %             EEG_W_beta = pop_eegfilt(EEG_W,13,30,[],[],1);
    %             EEG_W_gamma = pop_eegfilt(EEG_W,31,49,[],[],1);
                EEG_W = pop_eegfilt(EEG_W,Fbp1,Fbp2,[],[],1);

                Locs_M = {EEG_M.chanlocs.labels};
                Locs_W = {EEG_W.chanlocs.labels};
                [~,idx_M,idx_W] = intersect(Locs_M,Locs_W,'stable');
                DATA_M = double(EEG_M.data');
                DATA_W = double(EEG_W.data');
                DATA_M = DATA_M(:,idx_M);
                DATA_W = DATA_W(:,idx_W);

        s1 = size(DATA_M,1);
        s2 = size(DATA_W,1);
        min_s = min(s1,s2);
        DATA_M = DATA_M(s1-min_s+1:s1,:); % mistake
        DATA_W = DATA_W(s2-min_s+1:s2,:); % mistake

        for k = 1:32
            ch1 = DATA_M(:,k)';
            ch2 = DATA_W(:,k)';
            % ���㽻�����ܶ�
            [c12,f] = cpsd(ch1,ch2,[],[],500,250);

            % ���㹦�����ܶ�
            [p1,f] = pwelch(ch1,[],[],500,250);
            [p2,f] = pwelch(ch2,[],[],500,250);

            idx_delta = (f>=1 & f<= 3);
            idx_theta = (f>=4 & f<= 7);
            idx_alpha = (f>=8 & f<= 13);
            idx_beta = (f>=14 & f<= 30);
            idx_gamma = (f>=31 & f<= 45);

            m_c12_delta = c12(idx_delta);
            m_p1_delta = p1(idx_delta);
            m_p2_delta = p2(idx_delta);

            m_c12_theta = c12(idx_theta);
            m_p1_theta = p1(idx_theta);
            m_p2_theta = p2(idx_theta);

            m_c12_alpha = c12(idx_alpha);
            m_p1_alpha = p1(idx_alpha);
            m_p2_alpha = p2(idx_alpha);

            m_c12_beta = c12(idx_beta);
            m_p1_beta = p1(idx_beta);
            m_p2_beta = p2(idx_beta);

            m_c12_gamma = c12(idx_gamma);
            m_p1_gamma = p1(idx_gamma);
            m_p2_gamma = p2(idx_gamma);

            isc12_delta = abs(m_c12_delta).^2 ./ (m_p1_delta .* m_p2_delta);
            isc12_theta = abs(m_c12_theta).^2 ./ (m_p1_theta .* m_p2_theta);
            isc12_alpha = abs(m_c12_alpha).^2 ./ (m_p1_alpha .* m_p2_alpha);
            isc12_beta = abs(m_c12_beta).^2 ./ (m_p1_beta .* m_p2_beta);
            isc12_gamma = abs(m_c12_gamma).^2 ./ (m_p1_gamma .* m_p2_gamma);


            ch_isc12_delta(i,k) = mean(isc12_delta);
            ch_isc12_theta(i,k) = mean(isc12_theta);
            ch_isc12_alpha(i,k) = mean(isc12_alpha);
            ch_isc12_beta(i,k) = mean(isc12_beta);
            ch_isc12_gamma(i,k) = mean(isc12_gamma);
        end
        if count<10
            eval( ['isc_coh_alpha.sub0', num2str(count), ['=squeeze(ch_isc12_alpha(' num2str(i) ',:));']] );
            eval( ['isc_coh_beta.sub0', num2str(count), ['=squeeze(ch_isc12_beta(' num2str(i) ',:));']] );
            eval( ['isc_coh_delta.sub0', num2str(count), ['=squeeze(ch_isc12_delta(' num2str(i) ',:));']] );
            eval( ['isc_coh_gamma.sub0', num2str(count), ['=squeeze(ch_isc12_gamma(' num2str(i) ',:));']] );
            eval( ['isc_coh_theta.sub0', num2str(count), ['=squeeze(ch_isc12_theta(' num2str(i) ',:));']] );
            count = count + 1;
        else
            eval( ['isc_coh_alpha.sub', num2str(count), ['=squeeze(ch_isc12_alpha(' num2str(i) ',:));']] );
            eval( ['isc_coh_beta.sub', num2str(count), ['=squeeze(ch_isc12_beta(' num2str(i) ',:));']] );
            eval( ['isc_coh_delta.sub', num2str(count), ['=squeeze(ch_isc12_delta(' num2str(i) ',:));']] );
            eval( ['isc_coh_gamma.sub', num2str(count), ['=squeeze(ch_isc12_gamma(' num2str(i) ',:));']] );
            eval( ['isc_coh_theta.sub', num2str(count), ['=squeeze(ch_isc12_theta(' num2str(i) ',:));']] );
            count = count + 1;
        end
    end

    save(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\isc_coh_alpha'],'-struct','isc_coh_alpha');
    save(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\isc_coh_beta'],'-struct','isc_coh_beta');
    save(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\isc_coh_delta'],'-struct','isc_coh_delta');
    save(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\isc_coh_gamma'],'-struct','isc_coh_gamma');
    save(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\isc_coh_theta'],'-struct','isc_coh_theta');

h = 4;
    name = {'isc_coh_alpha','isc_coh_beta','isc_coh_delta','isc_coh_gamma','isc_coh_theta'};
    for i= 1:length(name)
         file_path= ['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Depart_set\ISC_coh_mat\' savename{h} '\' name{i}  '.mat'];
         load(file_path)

        % ��ȡ�ֶ���������ȡһ���ж��ٸ�sub
        eval(['n_len = length(fieldnames(' name{i} '));']);

        % ȫ��д�뵽data
        data=[];
        for j=1:n_len
%              eval( ['data(', num2str(j), [',:)=' name{i} '.form1.sub0' num2str(j) ';']] );
           if j<10
            eval( ['data(', num2str(j), [',:)=' name{i} '.sub0' num2str(j) ';']] );
           else
               eval( ['data(', num2str(j), [',:)=' name{i} '.sub' num2str(j) ';']] );
           end
        end

        % �ٰ�dataд��excel
        xlswrite(['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\correct\Result\' savename{h} '\' name{i} '.xlsx'],data);
    end
end
