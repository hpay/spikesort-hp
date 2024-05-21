for ii = 1:length(Tmissing)
    
    % 1. Copy all files (not folders and excluding .avi or .mp4 files) from fromdir to todir
    fromdir1 = fullfile('Z:\Hannah\ephys\project2', Tmissing{ii});
    todir1 = fullfile('D:\data', Tmissing{ii});
    
    % Get list of files in fromdir1
    files1 = dir(fullfile(fromdir1, '*.*'));
    files1 = files1(~[files1.isdir]);  % Filter out directories

    % Copy each file to todir1, excluding .avi or .mp4 files
    for jj = 1:numel(files1)
        % Exclude .avi and .mp4 files
        [~, ~, ext] = fileparts(files1(jj).name);
        if ~strcmpi(ext, '.avi') && ~strcmpi(ext, '.mp4')
            disp(files1(jj).name)
            copyfile(fullfile(fromdir1, files1(jj).name), todir1, 'f');
        end
    end

    % 2. Copy all files (not folders) from fromdir to todir
    fromdir2 = fullfile('Z:\Hannah\ephys\project2', Tmissing{ii}, 'kilosort2_output');
    todir2 = fullfile('D:\data', Tmissing{ii}, 'kilosort2_output');

    % Get list of files in fromdir2
    files2 = dir(fullfile(fromdir2, '*.*'));
    files2 = files2(~[files2.isdir]);  % Filter out directories

    % Copy each file to todir2
    for kk = 1:numel(files2)
        copyfile(fullfile(fromdir2, files2(kk).name), todir2, 'f');
    end
end