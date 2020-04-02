clear all
close all
clc

%% landline
% landline = [0 2 0 round(9*rand(1, 8))];
% disp(landline)

landline = [0 2 0 8 3 9 1 6 2 2 6];

%% Part 1
time = 0:1/32768:5.25-1/32768;
y = [];

for i = 1:11
    if(landline(i) >= 7)
        first(i) = 852;
    elseif(landline(i) >= 4)
        first(i) = 770;
    elseif(landline(i) >= 1)
        first(i) = 697;
    else
        first(i) = 941;
    end
    
    if(landline(i) == 1 || landline(i) == 4 || landline(i) == 7)
        second(i) = 1209;
    elseif(landline(i) == 3 || landline(i) == 6 || landline(i) == 9)
        second(i) = 1477;
    else
        second(i) = 1336;
    end
    
    y = [y sin(2*pi*first(i)*time((i-1)*0.5*32768+1 : ((i-1)*2+1)*0.25*32768)) + sin(2*pi*second(i)*time((i-1)*0.5*32768+1 : ((i-1)*2+1)*0.25*32768))];
    if (i ~= 11)
        y = [y 0*time(1 : 0.25*32768)];
    end
    
end

which_number = 4;
which_number = (which_number-1)*2;
figure()
plot(time(which_number*0.25*32768+1 : (which_number+3)*0.25*32768), y(which_number*0.25*32768+1 : (which_number+3)*0.25*32768))
title(sprintf('Time domain signal for the fourth and fifth digits of the landline number %.0f%.0f%.0f %.0f%.0f%.0f%.0f%.0f%.0f%.0f%.0f', landline(1), landline(2), landline(3), landline(4), landline(5), landline(6), landline(7), landline(8), landline(9), landline(10), landline(11)))
ylabel('Signal amplitude (Au)')
xlabel('Time (s)')
xlim([which_number*0.25 (which_number+3)*0.25])

%% Part 2
[s,f] = spectrogram(y, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis');

figure()
plot(time, y)
title(sprintf('Time domain signal for the landline number %.0f%.0f%.0f %.0f%.0f%.0f%.0f%.0f%.0f%.0f%.0f', landline(1), landline(2), landline(3), landline(4), landline(5), landline(6), landline(7), landline(8), landline(9), landline(10), landline(11)))
ylabel('Signal amplitude (Au)')
xlabel('Time (s)')
xlim([0 5.25])

figure()
spectrogram(y, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis')
title('Spectrogram of the dial tone signal')
ylim([0 2])

figure(), hold on
for i = 1:2:21
    plot(f, abs(s(:, i)))
    xlim([600 1600])
end
title('Frequency spectrum of the dial tone signal')
ylabel('Amplitude (Au')
xlabel('frequency (Hz)')

%% Noise
noise_low = randn(1, 172032);
noise_medium = 10 * randn(1, 172032);
noise_high = 25 * randn(1, 172032);

%% Part 4
show_signals = 0;
show_spectrogram = 1;
show_frequency = 1;

y_low = y + noise_low;
y_medium = y + noise_medium;
y_high = y + noise_high;

s_low = spectrogram(y_low, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis');
s_medium = spectrogram(y_medium, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis');
s_high = spectrogram(y_high, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis');

if (show_signals == 1) % plot the noisy signals
    figure()
    plot(time, y_low)
    title(sprintf('Time domain signal for %.0f%.0f%.0f %.0f%.0f%.0f%.0f%.0f%.0f%.0f%.0f with added noise', landline(1), landline(2), landline(3), landline(4), landline(5), landline(6), landline(7), landline(8), landline(9), landline(10), landline(11)))
    ylabel('Signal amplitude (Au)')
    xlabel('Time (s)')
    xlim([0.23 0.26])

    figure()
    plot(time, y_medium)
    title(sprintf('Time domain signal for %.0f%.0f%.0f %.0f%.0f%.0f%.0f%.0f%.0f%.0f%.0f with added noise', landline(1), landline(2), landline(3), landline(4), landline(5), landline(6), landline(7), landline(8), landline(9), landline(10), landline(11)))
    ylabel('Signal amplitude (Au)')
    xlabel('Time (s)')
    xlim([0.23 0.26])

    figure()
    plot(time, y_high)
    title(sprintf('Time domain signal for %.0f%.0f%.0f %.0f%.0f%.0f%.0f%.0f%.0f%.0f%.0f with added noise', landline(1), landline(2), landline(3), landline(4), landline(5), landline(6), landline(7), landline(8), landline(9), landline(10), landline(11)))
    ylabel('Signal amplitude (Au)')
    xlabel('Time (s)')
    xlim([0.23 0.26])
end

if (show_spectrogram == 1)
    figure()
    spectrogram(y_low, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis')
    title('Spectrogram of the dial tone signal with added noise')
    ylim([0 2])

    figure()
    spectrogram(y_medium, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis')
    title('Spectrogram of the dial tone signal with added noise')
    ylim([0 2])

    figure()
    spectrogram(y_high, hann(0.25*32768), 0, 0.25*32768, 32768, 'yaxis')
    title('Spectrogram of the dial tone signal with added noise')
    ylim([0 2])
end

if (show_frequency == 1)
    figure(), hold on
    for i = 1:2:21
        plot(f, abs(s_low(:, i)))
        xlim([600 1600])
    end
    title('Frequency spectrum of the dial tone signal with added noise')
    ylabel('Amplitude (Au)')
    xlabel('frequency (Hz)')

    figure(), hold on
    for i = 1:2:21
        plot(f, abs(s_medium(:, i)))
        xlim([600 1600])
    end
    title('Frequency spectrum of the dial tone signal with added noise')
    ylabel('Amplitude (Au)')
    xlabel('frequency (Hz)')

    figure(), hold on
    for i = 1:2:21
        plot(f, abs(s_high(:, i)))
        xlim([600 1600])
    end
    title('Frequency spectrum of the dial tone signal with added noise')
    ylabel('Amplitude (Au)')
    xlabel('frequency (Hz)')
end

