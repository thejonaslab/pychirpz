import numpy as np
import numba
import math 

def dtft(x, omegas):
    """
    evaluate the DTFT at the indicated points omega for the signal x
    
    x is assumed to run from [-n/2, n/2-1]
    """
    N = len(x)
    ns = np.arange(N)
    W = np.zeros((len(omegas), N), dtype=np.complex128)
    for wi, w in enumerate(omegas):
        W[wi, :] = np.exp(-1.0j * w * ns)
        
    return np.dot(W, x)

@numba.jit(nopython=True)
def nextpow2(n):
    """
    Return the smallest power of two greater than or equal to n.
    """
    return int(math.ceil(math.log(n)/math.log(2)))

# now try ourselves a chirp-z transform
@numba.jit
def chirpz(x, M, A, W):
    """
    chirp z transform per Rabiner derivation pp1256
    x is our (complex) signal of length N
    

    
    """
    N = len(x)
    L = 2**(nextpow2(N + M -1))  # or nearest power of two
    yn = np.zeros(L, dtype=np.complex128)
    for n in range(N):
        yn_scale =  A**(-n) * W**((n**2.0)/2.0)
        yn[n] = x[n] * yn_scale
    Yr = np.fft.fft(yn)
    
    vn = np.zeros(L, dtype=np.complex128)
    for n in range(M):
        vn[n] = W**((-n**2.0)/2.0)
        
    for n in range(L-N+1, L):
        vn[n] = W**(-((L-n)**2.0)/2.0)
        
    Vr = np.fft.fft(vn)
    
    Gr = Yr * Vr
    
    gk = np.fft.ifft(Gr)
    #gk = np.convolve(yn, vn)
    
    Xk = np.zeros(M, dtype=np.complex128)
    for k in range(M):
        g_scale = W**((k**2.0)/2.0) 
        Xk[k] = g_scale * gk[k]
        
    return Xk

@numba.jit
def chirpz2d(x, M, A, W):
    N = len(x)

    out = np.zeros((N, M), dtype=np.complex128)
    for i in range(N):
        out[i] = chirpz(x[i], M, A, W)
    out2d = np.zeros((M, M), dtype=np.complex128)

    for i in range(M):
        out2d[i] = chirpz(out[:, i], M, A, W)
    return out2d

@numba.jit
def fchirpz2d(x, M, A, W):
    """
    chirp z transform per Rabiner derivation pp1256
    x is our (complex) signal of length N
    assume x is square, output M will be square, dims are the same on all sides
    
    
    """
    N = len(x)
    L = 2**(nextpow2(N + M -1))  # or nearest power of two
    yn = np.zeros((L, L), dtype=np.complex128)
    ns = np.arange(N)
    ms = np.arange(M)
    
    yn_scale =  A**(-ns) * W**((ns**2.0)/2.0)
        
    yn[:N, :N] = x * np.outer(yn_scale, yn_scale)

    Yr = np.fft.fft2(yn)
    
    vn = np.zeros(L, dtype=np.complex128)
    for n in range(M-1):
        vn[n] = W**((-n**2.0)/2.0)
        
    for n in range(L-N+1, L):
        vn[n] = W**(-((L-n)**2.0)/2.0)
        
    Vr = np.fft.fft2(np.outer(vn, vn))
    
    Gr = Yr * Vr
    
    gk = np.fft.ifft2(Gr)
    
   
    Xk = W**((ms**2.0)/2.0) 
        
    return gk[:M, :M] * np.outer(Xk, Xk)

def zoom_fft(x, theta_start, step_size, M):
    """
    "zoomed" version of the fft, produces M step_sized samples
    around the unit circle starting at theta_start
    
    """
    A = np.exp(1j * theta_start)
    W = np.exp(-1j * step_size)
    
    return chirpz(x, M, A, W)
