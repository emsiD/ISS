% ======================
% ISS
% 2018/2019
% Projekt
% xdubec00
% ======================

% Task #1
[Signal,Fs] = audioread('xdubec00.wav');% nacitaj audio vstup
%riadkovy vektor
Signal=Signal';
% vzorkovacia frekvencia
Fs;
% dlzka vo vzorkoch
SignalL=length(Signal);
% dlzka v sekundach
T=SignalL/Fs;
%pocet binarnych symbolov
B=SignalL/16;


% Task #2
figure(2)
i = 8;
index = 1;
while(i <= SignalL)
  if(Signal(i) > 0)
    binarySymbols(index, 1) = 1;
  elseif (Signal(i) < 0)
    binarySymbols(index, 1) = 0;
  end
  i = i + 16;
  index=index+1;
end
plot((1:1:320)/16, Signal(1:320));
hold on
difSignal = binarySymbols(1:end);
binarySymbols = binarySymbols(1:20);
binSignal = binarySymbols;
stem((8:16:320)/16, binarySymbols, 'or');
xlabel ('t [ms]');
ylabel ('s[n], symbols');
ylim([-1,1]);
hold off
difTxt = textread('xdubec00.txt');
difference = xor(difSignal, difTxt);

% Task #3
figure(3)
B=[0.0192 -0.0185 -0.0185 0.0192];
A=[1.0000 -2.8870 2.7997 -0.9113];
zplane(B,A) 
if abs(roots(A)) < 1
    fprintf('Filter je stabilný.\n');
else
    fprintf('Filter nie je stabilný.\n');
end




% Task #4
figure(4)
plot((0 : 255) / 256 * Fs / 2, abs(freqz(B, A, 256)));
xlabel('f [Hz]');
ylabel('|H(f)|');
legend('Modul kmitoètovej charakteristiky');
fprintf('Dolná priepus.\n');
fprintf('Hranièná frekvencia leží na intervale 500 Hz\nUrèil som to pomocou priblíženia grafu.\n');
grid;






% Task #5
figure(5)
FSignal=filter(B, A, Signal);
plot(FSignal(1:320));
hold on
plot(Signal(1:320));
hold off
%graficky som zistil, ze treba zrealizovat posun do lava o 15




% Task #6
figure(6)
SFSignal=FSignal(15:end);
%SFSignal je posunuty FSignal
i = 8;
index = 1;
SFSignalL = length(SFSignal);
while(i <= SFSignalL)
  if(SFSignal(i) > 0)
    binarySymbols(index, 1) = 1;
  elseif (SFSignal(i) < 0)
    binarySymbols(index, 1) = 0;
  end
  i = i + 16;
  index=index+1;
end
plot((1:1:320)/16, Signal(1:320));
hold on
difSignal = binarySymbols(1:end);
binarySymbols = binarySymbols(1:20);
stem((8:16:320)/16, binarySymbols, 'or');
xlabel ('t');
ylabel ('s[n], symbols');
ylim([-1,1]);
plot((1:1:320)/16, FSignal(1:320));
plot((1:1:320)/16, SFSignal(1:320));

hold off


% Task #7
fprintf('\n7)\n');
errs = 0;
for i = 1:length(binarySymbols)
    errs = errs + xor(binSignal(i), binarySymbols(i));
end
chybovost = errs / 2000 * 100;
fprintf('Pocet chyb posunuteho signalu je: %4f.\n', errs); %pocet chyb
fprintf('Chybovost posunuteho signalu je: %3.2f %%.\n', chybovost);



% Task #8
figure(8)
Spektrum=fft(Signal);
plot(abs(Spektrum(1:Fs/2))); 
xlabel('f (Hz)'); 
ylabel('Magnitude');
hold on;
filtered_spectrum=fft(FSignal);
plot(abs(filtered_spectrum(1:Fs/2)));
hold off;


% Task #9
figure(9)
histogram=hist(Signal, 50);
plot(histogram/SignalL);
%xlim([-1 1]); neviem upravit x osu







% Task #10
figure(10)
k = (-50 : 50);
R = xcorr(Signal,'biased');
R = R(k + SignalL);
plot(k, R);
xlabel('k');
legend('Autokorelaèný koeficient');
grid;



% Task #11
%R(1), R(16), R[0]
fprintf('R[0] je %f.\n', R(51));
fprintf('R[1] je %f.\n', R(52));
fprintf('R[16] je %f.\n', R(67));

% Task #12
figure(12)
L = 50;
x = linspace(min(Signal), max(Signal), 50);
h = zeros(L, L);
[~, ind1] = min(abs(repmat(Signal(:)', L, 1) - repmat(x(:), 1, SignalL)));
ind2 = ind1(1 + 1 : SignalL);
for i = 1 : SignalL - 1,
	d1 = ind1(i);
	d2 = ind2(i);
	h(d1, d2) = h(d1, d2) + 1;
end
surf = (x(2) - x(1)) ^ 2;
p = h / (SignalL -1) / surf;
imagesc(x, x, p);
axis xy;
colorbar;
xlabel('x2');
ylabel('x1');



% Task #13
fprintf('\n13)\n');
check = sum(sum(p)) * surf;
fprintf('Overenie, èi se jedná o správnu združenú funkciu hustoty rozdelenia pravdepodobnosti, 2D integrál by mal by rovný 1 a je rovný %f.\n', check);


% Task #14
fprintf('\n14)\n');
fprintf('Hodnota koeficientu R[1] je %f.\n', sum(sum(repmat(x(:), 1, L) .* repmat(x(:)', L, 1) .* p)) * surf);



