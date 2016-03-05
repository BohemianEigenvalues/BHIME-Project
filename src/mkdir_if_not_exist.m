function mkdir_if_not_exist(dirpath)
    if dirpath(end) ~= '/', dirpath = [dirpath '/']; end
    if (exist(dirpath, 'dir') == 0), mkdir(dirpath); end
end 