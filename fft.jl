using WAV, Plots, Noise

function fft(f)
    N = length(f)
    if N == 1 return f end
    N_half = Int(N/2)
    w = exp(-2π*im/N)
    G = fft(f[1:2:end])
    H = fft(f[2:2:end])
    F = zeros(Complex{Float64}, N)
    wk = 1;
    for k in 1:N_half
	wHk = wk*H[k]
	F[k] = G[k] + wHk
	F[k + N_half] = G[k] - wHk
	wk = wk*w
    end
    return F
end

function ifft(a)
    a = conj.(a)
    a = fft(a)
    a = conj.(a)
    a = a./length(a)
end

N = 2^13
fs = 2^13 # sampling frequency
#fs = 44100

t = 0.0:1/fs:prevfloat(1.0) # time (prevfloat(1.0) = 0.99999999999999)
#t = 0:1/fs:(N - 1)/fs
f = 1e3
#f = 440.0
y = sin.(2pi * f * t) * 0.1 # 1kHz sine tone

noisy = add_gauss(y);
## write clean and noisy to file
wavwrite(y, "clean.wav", Fs=fs)
wavwrite(noisy, "noisy.wav", Fs=fs)

## compute DFT
fhat = fft(y.+0.0im)

# plot
p_original = plot(y, label = "clean", xlabel = "sampel index", ylabel="s[i]")
p_noisy = plot(noisy, label = "added noise", xlabel = "sampel index", ylabel="s[i]")
p_bins = plot(abs.(fhat)[1:4000], label = "fft", xlabel = "freq index", ylabel="mag",
              line = :stem, marker = :o, color=:black)

## 0 for elements below 100 amplitude and 1 otherwise
indices = abs.(fhat).>10;

## multiplying fhat by indices zeroes out the elements below 100 amplitude
fhat = fhat.*indices;
## invert to recover the signal
fixed = real(ifft(fhat))
p_fixed = plot(fixed, label = "fixed", xlabel = "sample index", ylabel = "s[i]")

p = plot(p_original,p_noisy,p_bins, p_fixed , layout = (2,2))
savefig(p, "plot.png")
display(p)

# write fixed audio to file
wavwrite(fixed, "fixed.wav", Fs=fs)
#yclean, fs = wavread("fixed.wav")

function fftshift(x)
    n = length(x)
    m = Int(n/2)
    return [x[m+1:end]; x[1:m]]
end
