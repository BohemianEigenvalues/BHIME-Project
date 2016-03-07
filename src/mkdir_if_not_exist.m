function mkdir_if_not_exist(dirpath)
    if dirpath(end) ~= filesep, dirpath = [dirpath, filesep]; end
    if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end
end 