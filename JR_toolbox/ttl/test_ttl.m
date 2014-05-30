%/media/DONNEES/JR KING bureau/Pro/Toolbox/LinuxParportServer$ sudo ./parallelPortServer
addpath('/media/DONNEES/JR KING bureau/Pro/Toolbox/LinuxParportServer/');
ParportTTL('Open', 'localhost');   
value=32;
while 1,
     value=value+1;
    ParportTTL('Set', value,100);
    disp(value);pause(1);
end
 ParportTTL('Close', 'localhost');