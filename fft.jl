using WAV, Plots, Noise

function fft(a)
    y1 = Any[]; y2 = Any[]
    n = length(a)
    if n ==1 return a end
    wn(n) = exp(-2*π*im/n)
    y_even = fft(a[1:2:end])
    y_odd = fft(a[2:2:end])
    w = 1
    for k in 1:Int(n/2)
        push!(y1, y_even[k] + w*y_odd[k])
        push!(y2, y_even[k] - w*y_odd[k])
        w = w*wn(n)
    end
    return vcat(y1,y2)
end

function ifft(a)
    a = conj.(a)
    a = fft(a)
    a = conj.(a)
    a = a./length(a)
end

fs = 2^13 # Frequency
t = 0.0:1/fs:prevfloat(1.0) # prevfloat(1.0) = 0.99999999999999
f = 1e3
y = sin.(2pi * f * t) * 0.1 # 1kHz sine tone

noisy = add_gauss(y);
## write clean and noisy to file
wavwrite(y, "clean.wav", Fs=fs)
wavwrite(noisy, "noisy.wav", Fs=fs)

fhat = fft(y.+0.0im)

# recall that the amplitude of the signal is encoded as the magnitude (use abs)
# of the complex numbers in fhat, whilst the phase is as the angle.
x = LinRange(0.0, 1.0, fs)
p_original = plot(x, y, label = "clean")
p_noisy = plot(x, noisy, label = "added noise")
p_bins = plot(x, abs.(fhat), label = "bins") # plot the amplitude

## 0 for elements below 100 amplitude and 1 otherwise
indices = abs.(fhat).>100;

## multiplying fhat by indices zeroes out the elements below 100 amplitude
fhat = fhat.*indices;
## invert to recover the signal
fixed = real(ifft(fhat))
p_fixed = plot(x, fixed, label = "fixed")

p = plot(p_original,p_noisy,p_bins, p_fixed , layout = (2,2))
savefig(p, "plot.png")
display(p)

# write fixed audio to file
wavwrite(fixed, "fixed.wav", Fs=fs)
yclean, fs = wavread("fixed.wav")

