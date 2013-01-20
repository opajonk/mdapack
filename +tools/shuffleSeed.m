function seed = shuffleSeed
% this method is copied from MATLAB R2011b, since older
% versions do not support 'Seed'='shuffle'

% Create a seed based on 1/100ths of a second, this repeats itself
% about every 497 days.

% Wait until the time changes enough to guarantee a unique seed for each call.
seed0 = mod(floor(now*8640000),2^31-1); % traditionally this was sum(100*clock)
for i = 1:100
    seed = mod(floor(now*8640000),2^31-1);
    if seed ~= seed0, break; end
    pause(.01); % smallest recommended interval
end
end