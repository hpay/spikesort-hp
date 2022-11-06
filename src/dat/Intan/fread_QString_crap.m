function a = fread_QString(fid)
a = '';
length = fread(fid, 1, 'uint32');
if length == hex2num('ffffffff')
 return;
end
length = length / 2; % convert length from bytes to 16-bit Unicode words
for i = 1:length
 a(i) = fread(fid, 1, 'uint16');
end
return