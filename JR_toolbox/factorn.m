function x = factorn(n)
if round(sqrt(n))^2 >= n
    x = [round(sqrt(n))^2 round(sqrt(n))^2];
elseif ceil(sqrt(n))*floor(sqrt(n)) >= n
  x = [floor(sqrt(n)) ceil(sqrt(n))];
else
    x = [ceil(sqrt(n)) ceil(sqrt(n))];
end
