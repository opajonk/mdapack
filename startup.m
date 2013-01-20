this_file = mfilename('fullpath');

[mdapack_path] = fileparts(this_file);

% Set paths
addpath( mdapack_path );                       % main directory

% initialize 3rd party toolboxes by executing all M-files in the "thirdparty" directory

s = dir(fullfile(mdapack_path, 'thirdparty'));
for i=3:length(s)
    if ~isa(s(i), 'dir')
        fname = fullfile(mdapack_path, 'thirdparty', s(i).name);
        [~, ~, ext] = fileparts(fname);
        if strcmp(ext, '.m') == 1
            run(fname);
        end
    end
end


% cleanup
clear this_file mdapack_path s fname ext i;