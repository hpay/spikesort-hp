% Function to output the list of contact order for use by ProbeInterface
%
% e.g. for H6 probe:
% mapping_to_device = [ 53,55,61,63,62,60,58,56,54,52,50,46,48,47,43,45,
% 57,59,42,40,38,36,34,44,32,30,28,26,24,22,23,27,51,49,0,2,4,6,8,7,10,12,
% 14,16,18,20,21,19,41,39,37,35,33,31,29,25,17,9,11,15,13,5,3,1 ]
% probe.set_device_channel_indices(mapping_to_device)
function printchannelmap
P = readtable('D:\hannah\Dropbox\alab\Analysis\SILICON PROBE MAP H6.xlsx');

P = sortrows(sortrows(P,'ShankPosition'),'ShankLetter')
% P = sortrows(P,'MolexChannel')

output  = [];
for ii = 1:height(P)
    output = [output num2str(P.IntanChannel(ii)) ','];    
end
output