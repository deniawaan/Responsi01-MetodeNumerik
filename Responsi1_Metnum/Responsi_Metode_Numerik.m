clear;
clc;
close all;

%% KASUS 1 - DERET TAYLOR DAN GALAT
lambda = 12;
v = 36;
x = 10;
h0 = 5;
t = linspace(0, 5, 501);
z = t .* v ./ lambda;
hEksak = h0 .* (cos(z) + exp(x));
nList = [23 24 25];

hPendekatan = zeros(length(nList), length(t));
mae = zeros(size(nList));
mrpe = zeros(size(nList));

for idx = 1:length(nList)
    n = nList(idx);

    cosTaylor = zeros(size(z));
    sukuCos = ones(size(z));
    for k = 0:n-1
        if k == 0
            sukuCos = ones(size(z));
        else
            sukuCos = sukuCos .* (-z.^2) ./ ((2*k-1)*(2*k));
        end
        cosTaylor = cosTaylor + sukuCos;
    end

    expTaylor = 0;
    sukuExp = 1;
    for k = 0:n-1
        if k == 0
            sukuExp = 1;
        else
            sukuExp = sukuExp * x / k;
        end
        expTaylor = expTaylor + sukuExp;
    end

    hPendekatan(idx, :) = h0 .* (cosTaylor + expTaylor);
    galatMutlak = abs(hPendekatan(idx, :) - hEksak);
    mae(idx) = mean(galatMutlak);
    mrpe(idx) = mean(galatMutlak ./ abs(hEksak)) * 100;
end

figure('Name', 'Kasus 1 - Perbandingan Gelombang');
plot(t, hEksak, 'k-', 'LineWidth', 1.8);
hold on;
plot(t, hPendekatan(1, :), '--', 'LineWidth', 1.2);
plot(t, hPendekatan(2, :), '-.', 'LineWidth', 1.2);
plot(t, hPendekatan(3, :), ':', 'LineWidth', 1.8);
grid on;
xlabel('Waktu, t');
ylabel('Ketinggian, h');
title('Perbandingan Fungsi Eksak dan Deret Taylor');
legend('Eksak', 'n = 23', 'n = 24', 'n = 25', 'Location', 'best');

fprintf('\n=== KASUS 1: DERET TAYLOR DAN GALAT ===\n');
fprintf('Rentang t: 0 sampai 5; lambda = %g, v = %g, x = %g, h0 = %g\n', lambda, v, x, h0);
fprintf('Jumlah suku\tMAE\t\t\tGalat relatif rata-rata (%%)\n');
for idx = 1:length(nList)
    fprintf('%d\t\t%.8f\t%.8f\n', nList(idx), mae(idx), mrpe(idx));
end
fprintf('Toleransi galat relatif rata-rata = 1%%\n');

%% KASUS 2 - SPL DENGAN METODE GAUSS-SEIDEL
A = [4.00, 0.75, 1.00, 2.00;
     1.00, 4.25, 2.00, 0.90;
     2.50, 1.50, 5.00, 0.50;
     2.00, 1.50, 0.50, 6.00];
b = [149000; 173000; 161000; 251500];

xGS = zeros(4, 1);
tol = 1e-8;
maxIter = 1000;
konvergen = false;

for iter = 1:maxIter
    xLama = xGS;

    for i = 1:length(b)
        jumlahKiri = 0;
        jumlahKanan = 0;
        if i > 1
            jumlahKiri = A(i, 1:i-1) * xGS(1:i-1);
        end
        if i < length(b)
            jumlahKanan = A(i, i+1:end) * xLama(i+1:end);
        end
        xGS(i) = (b(i) - jumlahKiri - jumlahKanan) / A(i, i);
    end

    galatIterasi = norm(xGS - xLama, inf);
    if galatIterasi < tol
        konvergen = true;
        break;
    end
end

fprintf('\n=== KASUS 2: SPL METODE GAUSS-SEIDEL ===\n');
fprintf('Urutan variabel: beras, telur, tepung, ayam.\n');
if konvergen
    statusGS = 'Ya';
else
    statusGS = 'Tidak';
end
fprintf('Konvergen: %s; iterasi: %d; galat iterasi akhir: %.12g\n', statusGS, iter, galatIterasi);
fprintf('Harga beras  = Rp %.2f per kg\n', xGS(1));
fprintf('Harga telur  = Rp %.2f per kg\n', xGS(2));
fprintf('Harga tepung = Rp %.2f per kg\n', xGS(3));
fprintf('Harga ayam   = Rp %.2f per kg\n', xGS(4));
fprintf('Norma residu maksimum = %.12g\n', norm(A*xGS-b, inf));

%% KASUS 3 - PERSAMAAN NONLINEAR DENGAN METODE BISEKSI
f = @(r) (4/3) .* pi .* r.^3 - 300;
rGrafik = linspace(0, 6, 401);

figure('Name', 'Kasus 3 - Grafik Persamaan Nonlinear');
plot(rGrafik, f(rGrafik), 'LineWidth', 1.8);
hold on;
plot(rGrafik, zeros(size(rGrafik)), 'k--', 'LineWidth', 1.1);
grid on;
xlabel('Jari-jari, r (cm)');
ylabel('f(r)');
title('Grafik f(r) = (4/3) pi r^3 - 300');

batasBawah = 4;
batasAtas = 5;
tolBiseksi = 1e-10;
maxIterBiseksi = 200;
fa = f(batasBawah);
fb = f(batasAtas);

if fa * fb > 0
    error('Interval awal tidak mengurung akar. Pilih interval dengan tanda f(a) dan f(b) berbeda.');
end

konvergenBiseksi = false;
for iterBiseksi = 1:maxIterBiseksi
    rTengah = (batasBawah + batasAtas) / 2;
    ft = f(rTengah);

    if abs(ft) < tolBiseksi || (batasAtas - batasBawah)/2 < tolBiseksi
        konvergenBiseksi = true;
        break;
    end

    if fa * ft < 0
        batasAtas = rTengah;
        fb = ft;
    else
        batasBawah = rTengah;
        fa = ft;
    end
end

r = rTengah;
tinggiTabung = 2 * r;
volumeBola = (4/3) * pi * r^3;
volumeTabung = pi * r^2 * tinggiTabung;
volumeLuarBola = volumeTabung - volumeBola;

fprintf('\n=== KASUS 3: METODE BISEKSI ===\n');
fprintf('f(4) = %.8f; f(5) = %.8f\n', f(4), f(5));
if konvergenBiseksi
    statusBiseksi = 'Ya';
else
    statusBiseksi = 'Tidak';
end
fprintf('Konvergen: %s; iterasi: %d; residu f(r): %.12g\n', statusBiseksi, iterBiseksi, f(r));
fprintf('Jari-jari bola dan tabung = %.10f cm\n', r);
fprintf('Tinggi tabung = %.10f cm\n', tinggiTabung);
fprintf('Volume bola = %.10f cm^3\n', volumeBola);
fprintf('Volume tabung = %.10f cm^3\n', volumeTabung);
fprintf('Volume tabung di luar bola = %.10f cm^3\n', volumeLuarBola);
