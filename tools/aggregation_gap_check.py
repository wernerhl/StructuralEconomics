import numpy as np
rng=np.random.default_rng(0)
def gap(W):
    ev=np.sort(np.abs(np.linalg.eigvals(W)))[::-1]; return 1-ev[1]
def stat(W):
    w,v=np.linalg.eig(W.T); i=np.argmin(abs(w-1)); p=np.real(v[:,i]); return p/p.sum()
# (a) balanced closed T: D^-1 T vs T D^-1 spectra
# build balanced T: circulation = sum of random cycles
n=6; T=np.zeros((n,n))
for _ in range(30):
    k=rng.integers(2,n+1); c=rng.permutation(n)[:k]; w=rng.random()
    for a,b in zip(c,np.roll(c,-1)): T[a,b]+=w
y=T.sum(1); assert np.allclose(y,T.sum(0))
R=T/y[:,None]; C=(T/y[None,:]).T
print("a) gap receipt",gap(R)," gap expend",gap(C))
# (c) stationary-weighted collapse, non-reversible: count gap decreases
cnt=0;tot=0;worst=0
for trial in range(20000):
    n=6; K=rng.random((n,n))**3
    W=K/K.sum(1,keepdims=True)
    blocks=[[0,1,2],[3,4,5]]
    p=stat(W)
    Wb=np.zeros((2,2))
    for a,B in enumerate(blocks):
        for b,Bp in enumerate(blocks):
            Wb[a,b]=sum(p[i]*W[i,j] for i in B for j in Bp)/p[B].sum()
    g,gb=gap(W),gap(Wb); tot+=1
    if gb<g-1e-9: cnt+=1; worst=max(worst,g-gb)
print("c) nonrev collapse lowered gap in",cnt,"of",tot,"max drop",worst)
# reversible check
cnt=0
for trial in range(5000):
    K=rng.random((6,6)); K=K+K.T; W=K/K.sum(1,keepdims=True); p=stat(W)
    Wb=np.array([[sum(p[i]*W[i,j] for i in B for j in Bp)/p[B].sum() for Bp in blocks] for B in blocks])
    if gap(Wb)<gap(W)-1e-9: cnt+=1
print("c) reversible collapse lowered gap in",cnt,"of 5000 (note: |lambda2| gap, not 1-lambda2)")
