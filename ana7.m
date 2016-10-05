CHANNEL = 7;
TAPS = 67;
alpha = [270 210 150 90 30 -30 0];
fl = fopen("filter.dat", "rb");
j = fread(fl, TAPS*CHANNEL, 'float32');
fclose(fl);
val=[];
for i = 0:CHANNEL-1
    val = [val j(i*TAPS+1:(i+1)*TAPS)];
end;
val = val';
j = ones(1, CHANNEL);
h0 = j*val;

Fs = 16000;
cs = 34000;
radius = 4.25;

freq = 100:100:8000;
freq = 500;
w = 2*pi*(freq / Fs);

plane = [];
for theta = 0:5:360
    output = w*0;
    pz = exp(i*(0:-1:(1-TAPS))'*w);
    for ch=1:CHANNEL-1
        dz = Fs*radius/cs*cos((theta - alpha(ch))*pi/180.0);
        hi = val(ch,:);
        dl = exp(-i*dz*w);
        hm = dl.*(hi*pz);
        output = output + hm;
    end;
    hi = val(CHANNEL,:);
    hm = hi*pz;
    output = output + hm;
    plane = [plane abs(output)];
end;

%plot(plane);
polar((0:5:360)*pi/180.0, plane/360);
