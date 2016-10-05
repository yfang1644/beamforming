CHANNEL = 63;
TAPS = 67;
fl = fopen("filter.dat", "r");
val = fread(fl, [CHANNEL, TAPS], 'single');
fclose(fl);
j = ones(1, CHANNEL);
h0 = j*val;

Fs = 16000;
cs = 34000;
Dx = 35.0; dx = 35.0/8;
Dy = 25.0; dy = 25.0/6;

freq = 100:100:8000;
freq = 3000;
w = 2*pi*(freq / Fs);

plane = zeros(35, 35);
for theta = -85:5:85
    for phi = -85:5:85
        output = w*0;
        for ch=0:CHANNEL-1
            jch = floor(ch/9);
            ich = ch - jch*9;

            t = (dx*ich*sin(pi*theta/180) + dy*jch*sin(pi*phi/180))/cs;
            dz = t*Fs;

            hi = val(ch+1,:);
            pz = exp(i*(0:-1:(1-TAPS))'*w);
            dl = exp(-i*dz*w);
            hm = dl.*(hi*pz);
            output = output + hm;
        end;
        plane(theta/5+18, phi/5+18) = abs(output);
    end;
end;
surface(-85:5:85, -85:5:85, plane);
