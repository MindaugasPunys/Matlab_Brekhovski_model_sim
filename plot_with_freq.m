function [freqArray, ampArray] = plot_with_freq(freq, amp)

changeIndex = find(diff(sign(freq)) < 0, 1, 'first') + 1;
freqArray = [freq(changeIndex:end), freq(1:changeIndex-1)];
ampArray = [amp(changeIndex:end), amp(1:changeIndex-1)];

end