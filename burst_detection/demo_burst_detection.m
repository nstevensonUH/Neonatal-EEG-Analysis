% Demo_burst_detection

% SELECT FILE
filename = 'eeg_burst_eg.edf';
% DEFINE MONTAGE OF INTEREST
str = cell(18,2); 
str{1,1} = 'Fp2'; str{1,2} = 'F4';  
str{2,1} = 'F4'; str{2,2} = 'C4';    
str{3,1} = 'C4'; str{3,2} = 'P4';    
str{4,1} = 'P4'; str{4,2} = 'O2';   
str{5,1} = 'Fp1'; str{5,2} = 'F3';  
str{6,1} = 'F3'; str{6,2} = 'C3';    
str{7,1} = 'C3'; str{7,2} = 'P3';    
str{8,1} = 'P3'; str{8,2} = 'O1';    
str{9,1} = 'Fp2'; str{9,2} = 'F8';    
str{10,1} = 'F8'; str{10,2} = 'T4';    
str{11,1} = 'T4'; str{11,2} = 'T6';    
str{12,1} = 'T6'; str{12,2} = 'O2';   
str{13,1} = 'Fp1';  str{13,2} ='F7';  
str{14,1} = 'F7'; str{14,2} = 'T3';     
str{15,1} = 'T3'; str{15,2} = 'T5';  
str{16,1} = 'T5'; str{16,2} = 'O1'; 
str{17,1} = 'Fz'; str{17,2} = 'Cz';   
str{18,1} = 'Cz';  str{18,2} ='Pz';   

% READ FILE WITH RELEVANT MONTAGE
[eeg_data, label, fs] = read_edf_new(filename, str);
t1 = 0:1/fs(1):length(eeg_data{1})/fs(1)-1/fs(1);
% PERFORM BURST DETCTION PER CHANNEL
for ii = 1:length(eeg_data)
    [t2, ba(ii,:)] = detector_per_channel_palmu(eeg_data{ii}, fs(ii));
end

% PLOT RESULTS
figure; hold on; 
for ii = 1:length(eeg_data)
        plot(t1, eeg_data{ii}+100*ii); plot(t2, 50*ba(ii,:)+100*ii)
end
axis([0 max(t1) 0 1950])
set(gca, 'Ytick', [1:length(eeg_data)]*100, 'Yticklabel', label)
xlabel('time (s)')
