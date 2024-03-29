clc;
clear all;
%Precision = 'float';
Precision = 'double';
fidp = fopen ('output/nT.dat', 'r', 'l');
if (fidp == -1)
    disp('File "nT.dat" not found');
    return;
end
datap = fread (fidp, 1, 'int');
fclose (fidp);
n = datap(1);
Size = [n];
fid = fopen ('output/dT.dat', 'r', 'l');
if (fid == -1)
    disp('File "dT.dat" not found');
    return;
end
U = fread (fid, Size, Precision);
fclose (fid);
semilogy(U);