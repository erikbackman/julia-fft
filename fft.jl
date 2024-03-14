using WAV, Plots, Noise

function fft(a)
    y1 = Any[]; y2 = Any[]
    n = length(a)
    if n == 1 return a end
    wn(n) = exp(-2*π*im/n) # twiddle factor
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

# in place impl
function fft2(p)
    n = length(p)
    if n == 1 return p end
    w = exp(-2*pi*im/n)
    y_even = fft2(p[1:2:end])
    y_odd = fft2(p[2:2:end])
    y = zeros(Complex{Float64}, n)
    wk = 1;
    for k in 1:Int(n/2)
        y[k] = y_even[k] + wk*y_odd[k]
        y[k + Int(n/2)] = y_even[k] - wk*y_odd[k]
        wk = wk*w
    end
    return y
end

function ifft(a)
    a = conj.(a)
    a = fft(a)
    a = conj.(a)
    a = a./length(a)
end

fs = 2^13 # sampling frequency
t = 0.0:1/fs:prevfloat(1.0) # time (prevfloat(1.0) = 0.99999999999999)
f = 1e3
y = sin.(2pi * f * t) * 0.1 # 1kHz sine tone

noisy = add_gauss(y);
## write clean and noisy to file
#wavwrite(y, "clean.wav", Fs=fs)
#wavwrite(noisy, "noisy.wav", Fs=fs)

#fhat = fft(y.+0.0im)
fhat = fft2(y.+0.0im)

# recall that the amplitude of the signal is encoded as the magnitude (use abs)
# of the complex numbers in fhat, whilst the phase as the angle.
x = LinRange(0.0, 1.0, fs)
p_original = plot(y, label = "clean", xlabel = "sampel index", ylabel="s[i]")
p_noisy = plot(noisy, label = "added noise", xlabel = "sampel index", ylabel="s[i]")
p_bins = plot(abs.(fhat)[1:4000], label = "fft", xlabel = "freq index", ylabel="mag") # plot the amplitude

## 0 for elements below 100 amplitude and 1 otherwise
indices = abs.(fhat).>20;

## multiplying fhat by indices zeroes out the elements below 100 amplitude
fhat = fhat.*indices;
## invert to recover the signal
fixed = real(ifft(fhat))
p_fixed = plot(fixed, label = "fixed", xlabel = "sample index", ylabel = "s[i]")

p = plot(p_original,p_noisy,p_bins, p_fixed , layout = (2,2))
savefig(p, "plot.png")
display(p)

# write fixed audio to file
#wavwrite(fixed, "fixed.wav", Fs=fs)
#yclean, fs = wavread("fixed.wav")

function fftshift(x)
    n = length(x)
    m = Int(n/2)
    return [x[m+1:end]; x[1:m]]
end
