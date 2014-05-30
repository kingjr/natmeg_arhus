function out = science_notation(in,precision)
if nargin == 1, precision = 3; end

s = size(in);
out = reshape(in,[],1);
for ii = 1:numel(in)
    if abs(log10(in{ii})) > precision
        if in{ii}>0
            out{ii,1} = '10';
            out{ii,2} = num2str(ceil(log10(in{ii})));
        elseif in{ii}==0
            out{ii,2} = '0';
        else
            out{ii,1} = '-10';
            out{ii,2} = num2str(ceil(log10(abs(in{ii}))));
        end
    else
        out{ii,1} = num2str(round(in{ii} *10^precision)/10^precision,precision);
    end
end
out = reshape(permute(reshape(out,cat(2,s,2)),[1 length(s)+1 2:length(s)]),[s(1) s(2)*2 s(3:end)]);