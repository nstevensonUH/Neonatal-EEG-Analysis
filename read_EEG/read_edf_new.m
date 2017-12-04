function [dat, label, fs] = read_edf_new(filename, varargin)
% [data, hdr, label] = read_edf_new(filename);
%
% This functions reads in a EDF file as per the format outlined in
%  http://www.edfplus.info/specs/edf.html
%
% INPUT: filename - EDF file name
%               varargin - specfic montage required (if there is not input then the recorded montage is used)
%
% OUTPUT: dat - a cell array containing the data in the file (int16 format)
%                   label - a cell array containing the labels for each
%                   channel in dat (ASCII format)
%                   fs - the sampling frequency
%
%
% Nathan Stevenson
% Univsirty of Helsinki, Finland
% 13/10/2017
%

fid = fopen(filename, 'r');
hdr = cell(1);
hdr{1} = fread(fid, 256, 'char');         % CONTAINS PATIENT INFORMATION, RECORDING INFORMATION
len_s = str2num(char(hdr{1}(235:244))');        % START DATE AND TIME and a RESERVED
rec_dur = str2num(char(hdr{1}(244:252))');
ns = char(hdr{1}(253:256))';
ns = str2num(ns);
hdr{2} = fread(fid, ns*16, 'char');    % LABEL channel label, temp or HR
hdr{3} = fread(fid, ns*80,'char');     % TRANSDUCER TYPE
hdr{4} = fread(fid, ns*8,'char');       % PHYSICAL DIMENSION, voltage - temperature
hdr{5} = fread(fid, ns*8,'char');       % PHYSICAL MIN
hdr{6} = fread(fid, ns*8,'char');       % PHYSICAL MAX
hdr{7} = fread(fid, ns*8,'char');       % DIGITAL MIN
hdr{8} = fread(fid, ns*8,'char');       % DIGITAL MAX
label = cell(1);
for jj=1:ns;
rf2 = jj*8; rf1 = rf2-7; 
label{jj} = char(hdr{2}(16*(jj-1)+1:16*(jj-1)+16));
phy_lo(jj) = str2num(char(hdr{5}(rf1:rf2))');
phy_hi(jj) = str2num(char(hdr{6}(rf1:rf2))');  
dig_lo(jj) = str2num(char(hdr{7}(rf1:rf2))');
dig_hi(jj) = str2num(char(hdr{8}(rf1:rf2))');
end
scle = (phy_hi-phy_lo)./(dig_hi-dig_lo);
offs = (phy_hi+phy_lo)/2;                                   

hdr{9} = fread(fid, ns*80,'char');                    % PRE FILTERING
hdr{10} = fread(fid, ns*8, 'char');                 % SAMPLING NO rec
nsamp = str2num(char(hdr{10})');
hdr{11} = fread(fid, ns*32,'char');     % RESERVED    
fs = nsamp/rec_dur;

% Build the empty data matrix of size INT16 not double
dat = cell(1, ns);
for jj = 1:length(nsamp);
 %   len_s.*nsamp(jj)
       dat{jj} = double(int16(zeros(1,len_s*nsamp(jj)))).*scle(jj)-offs(jj);        
end

% Load data into dat array from EDF file: there are length(nsamp) channels
% and the size of each channel will be len_s * nsamp(ii)
for ii = 1:len_s;
    for jj = 1:length(nsamp);
        r1 = nsamp(jj)*(ii-1)+1; r2 = ii*nsamp(jj);
        dat{jj}(r1:r2) = double(fread(fid, nsamp(jj), 'short')').*scle(jj)-offs(jj);    
    end
end

if length(varargin)==1 
    str = varargin{1};
    
    data_mont = cell(1,length(str));
    for jj = 1:length(str)
    ref1 = zeros(1,length(dat));
    ref2 = zeros(1,length(dat));
    for ii = 1:length(dat);
        ref1(ii) = length(findstr(label{ii}', str{jj,1})); 
        ref2(ii) = length(findstr(label{ii}', str{jj,2}));
    end
    
    qq1 = find(ref1==1,1);
    qq2 = find(ref2==1,1);
    if length(dat{qq1})~=length(dat{qq2})
        data_mont{jj} = zeros(1,length(dat{1}));
    else
    data_mont{jj} = dat{qq1}-dat{qq2}; 
    end
    end
       
    dat = data_mont;

    for jj = 1:length(str)
            label{jj} = [str{jj,1} '-' str{jj,2}];
    end
    
end


