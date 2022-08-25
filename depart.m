 clear all;
clc
path = 'F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\';
addpath(path);
%�����������������ļ���
savename = {'resting','pol','poh','ppl','pph','no7','in7','now','inw','mol','moh','mpl','mph'};
%����ÿ�ζ�Ӧ���ļ���
% for h = 1:13
%     mkdir(['E:\analysis\results\epoch\correct_set\' savename{h} ]);
% end

% savename2 = {};
for i =  1:34
file_path = ['F:\�������\�ٶ�����\��������\qEEG - BBI - BBI\WFWX\H' num2str(i) 'W' num2str(i)];
    file_name_1 = ['H' num2str(i) '.vhdr']; file_name_2 = [ 'W' num2str(i) '.vhdr'];
%    eval(['[file' num2str(i) '_H,~] = pop_loadbv(file_path,file_name_1);']);
    %��ȡ����Ů�Եķֶ���Ϣ������
    for j =1:2
        %eeglab�����Ե������ļ�
       eval(['[file' num2str(i) '_H,~] = pop_loadbv(file_path,file_name_' num2str(j) ');']);
       eval(['[ALLEEG,~] = pop_loadbv(file_path,file_name_' num2str(j) ');']);
        
        %��������Ϊ250
        ALLEEG = pop_resample(ALLEEG,250);
         %��ȡevent�ļ����latency
        if j == 1
             eventname = ['E:\analysis\results\event_excle\file' num2str(i) '_H' num2str(i) '.xlsx'];
        else
             eventname = ['E:\analysis\results\event_excle\file' num2str(i) '_W' num2str(i) '.xlsx'];
        end
        
         %��Ϣ̬
         num_resting = xlsread(eventname,1,'A3:A4');
         begin_resting = num_resting(1,1);
         stop_resting = num_resting(2,1);
         %ͼƬ����
         %�����ͻ���
         num_pol = xlsread(eventname,1,'A5:A6');
         begin_pol = num_pol(1,1);
         stop_pol = num_pol(2,1);
         %�����߻���
          num_poh = xlsread(eventname,1,'A7:A8');
         begin_poh = num_poh(1,1);
         stop_poh = num_poh(2,1);
         %�����ͻ���
           num_ppl = xlsread(eventname,1,'A9:A10');
         begin_ppl = num_ppl(1,1);
         stop_ppl = num_ppl(2,1);
          %�����߻���
           num_pph = xlsread(eventname,1,'A11:A12');
         begin_pph = num_pph(1,1);
         stop_pph = num_pph(2,1);
           %���ɰ�
           %�޻���
          num_no7 = xlsread(eventname,1,'A13:A14');
         begin_no7 = num_no7(1,1);
         stop_no7 = num_no7(2,1);
         %�л���
           num_in7 = xlsread(eventname,1,'A15:A16');
         begin_in7 = num_in7(1,1);
         stop_in7 = num_in7(2,1); 
         %����
           %�޻���
          num_now = xlsread(eventname,1,'A17:A18');
         begin_now = num_now(1,1);
         stop_now = num_now(2,1);
            %�л���
          num_inw = xlsread(eventname,1,'A19:A20');
         begin_inw = num_inw(1,1);
         stop_inw = num_inw(2,1);
         %��������
          %�����ͻ���
         num_mol = xlsread(eventname,1,'A21:A22');
         begin_mol = num_mol(1,1);
         stop_mol = num_mol(2,1);
           %�����߻���
         num_moh = xlsread(eventname,1,'A23:A24');
         begin_moh = num_moh(1,1);
         stop_moh = num_moh(2,1);
           %�����ͻ���
         num_mpl = xlsread(eventname,1,'A25:A26');
         begin_mpl = num_mpl(1,1);
         stop_mpl = num_mpl(2,1);
           %�����߻���
         num_mph = xlsread(eventname,1,'A27:A28');
         begin_mph = num_mph(1,1);
         stop_mph = num_mph(2,1);
         
         if i<10 %����resting01-resting39
             if j==1 %j=1ʱfile_name_1Ϊ���ԣ�j=2ʱfile_name_2ΪŮ��
                 %ѭ���ֶβ�����set
                 for k =1:13 %�Ե����ݷ�13��
                     eval(['begin = begin_' savename{k} ';']);
                     eval(['stop = stop_' savename{k} ';']);
                     part = pop_select(ALLEEG,'point',[begin stop]);
                     part = pop_saveset(part,['E:\analysis\results\epoch\' savename{k} '\' savename{k} '0' num2str(i) '_H']);
                 end
                 
             else
                 for k = 1:13
                     eval(['begin = begin_' savename{k} ';']);
                     eval(['stop = stop_' savename{k} ';']);
                     part = pop_select(ALLEEG,'point',[begin stop]);
                     part = pop_saveset(part,['E:\analysis\results\epoch\' savename{k} '\' savename{k} '0' num2str(i) '_W']);
                 end  %for_k
             end %if_j
         else
             if j==1 %j=1ʱfile_name_1Ϊ���ԣ�j=2ʱfile_name_2ΪŮ��
                 %ѭ���ֶβ�����set
                 for k =1:13 %�Ե����ݷ�13��
                     eval(['begin = begin_' savename{k} ';']);
                     eval(['stop = stop_' savename{k} ';']);
                     part = pop_select(ALLEEG,'point',[begin stop]);
                     %            savename{k} = pop_select(ALLEEG,'point',['begin_' savename{k} 'stop_' savename{k} ]);
                     part = pop_saveset(part,['E:\analysis\results\epoch\' savename{k} '\' savename{k}  num2str(i) '_H']);
                 end
                 
             else
                 for k = 1:13
                     eval(['begin = begin_' savename{k} ';']);
                     eval(['stop = stop_' savename{k} ';']);
                     part = pop_select(ALLEEG,'point',[begin stop]);
                     part = pop_saveset(part,['E:\analysis\results\epoch\' savename{k} '\' savename{k}  num2str(i) '_W']);
                 end  %for_k
             end %if_j
             
         end%if_i<10
    end
end



