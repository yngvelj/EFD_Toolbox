function powerSpectra(data,x,d,U,fig)
% PIV data:
y = data(1).y;
yl = data(1).yl; %number of samples in y
nsamples = length(data);%number of samples (sample time/samplefreq)

% Tolerance is half of grid size
tolx = (abs(data(1).x(end))-abs(data(1).x(end-1)))/2;

%Find x indexes: of line
xind = find(abs(data(1).x-x) < tolx);

%initialize x,y,z vectors for powespectra plotting
xvec = zeros(1,nsamples);
yvec = zeros(1,nsamples);

linespecX = zeros(yl,round(nsamples/2));
linespecY = zeros(yl,round(nsamples/2));

for j = 1:yl
    for i = 1:nsamples
        xvec(i) = data(i).velXArray(j,xind);
        yvec(i) = data(i).velYArray(j,xind);
    end
    %Create X velocity spectrum
    Y = fft(xvec);           
    P2 = abs(Y/nsamples);   %Two sided spectrum
    P1 = P2(1:ceil(nsamples/2));%One sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);
    
    linespecX(j,:) = P1;
    
    %Create Y velocity spectrum
    Y = fft(yvec);           
    P2 = abs(Y/nsamples);
    P1 = P2(1:ceil(nsamples/2));
    P1(2:end-1) = 2*P1(2:end-1);
    
    linespecY(j,:) = P1;
end


Fs = 15;%sampling freq
f = Fs*(0:(nsamples/2))/nsamples;
%make sure its same length
if length(f) == length(linespecX(1,:))+1
    f = f(1:end-1);
end
figure(fig)
plot(f,P1); 

%X-spectrum
figure(fig+1)
hold on;
title("Spectrum of X - velocity");
xlabel('$\frac{f D}{U}$','interpreter','latex');
ylabel('$\frac{y}{d_{twine}}$','interpreter','latex');
surf(f.*(d*1e-3)./U,y./d,linespecX);

%Y-spectrum
figure(fig+2);
hold on;
title("Spectrum of Y - velocity");
xlabel('$\frac{f D}{U}$','interpreter','latex');
ylabel('$\frac{y}{d_{twine}}$','interpreter','latex');
surf(f,y./d,linespecY);

end