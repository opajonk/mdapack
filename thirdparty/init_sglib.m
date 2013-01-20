[sglib_path] = fullfile(fileparts(mfilename('fullpath')), 'sglib');

if isdir(sglib_path)
    run(fullfile(sglib_path, 'startup.m'));
end