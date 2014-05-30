function [s1 s2] = asymetries(signal)
symbols_1 = [1 2 3; 3 1 2; 2 3 1]; % gowing up
symbols_2 = [3 2 1; 2 1 3; 1 3 2]; % gowing down

[unused x] = sort(cat(1,...
    signal(1:(end-2)),...
    signal(2:(end-1)),...
    signal(3:(end))));
x= x';
s1 = length(strmatch(symbols_1(1,:),x)) + length(strmatch(symbols_1(2,:),x));
s2 = length(strmatch(symbols_2(1,:),x)) + length(strmatch(symbols_2(2,:),x));
end