"""Dissimulation des blocs manquants dans le domaine temporel (répétition de la forme d'onde précédente
avec fondus), au lieu de moyenner des spectres TDAC."""
import numpy as np
def conceal(audio, missing, hop=512, span=768, fade=96):
    """audio : (ch, n) synthétisé avec des spectres nuls pour les blocs manquants ; missing : indices de blocs."""
    out=audio.copy(); n=audio.shape[1]
    bad=np.zeros(n,bool)
    for t in missing: bad[t*hop:min(n,t*hop+span)]=True
    # zones contiguës à reconstruire
    i=0
    while i<n:
        if not bad[i]: i+=1; continue
        j=i
        while j<n and bad[j]: j+=1
        L=j-i; src0=i-L-fade
        if src0>=0:
            seg=out[:,src0:src0+L+2*fade].copy()           # forme d'onde juste avant le trou
        else:
            seg=np.zeros((out.shape[0],L+2*fade))
        # fondu d'entrée et de sortie avec le signal existant
        w=np.ones(L+2*fade); r=np.linspace(0,1,fade)
        a=max(0,i-fade); b=min(n,j+fade)
        piece=seg[:,:b-a]
        win=np.ones(b-a)
        k=i-a
        if k>0: win[:k]=np.linspace(0,1,k)
        k2=b-j
        if k2>0: win[-k2:]=np.linspace(1,0,k2)
        out[:,a:b]=out[:,a:b]*(1-win)+piece*win
        i=j
    return out

def conceal2(audio, missing, hop=512, span=768, fade=128):
    """Comme conceal, mais en interpolant des deux côtés : la forme d'onde qui précède le trou s'estompe
    pendant que celle qui le suit (décalée vers l'arrière) apparaît. Plus lisse quand le trou est court."""
    out=audio.copy(); n=audio.shape[1]
    bad=np.zeros(n,bool)
    for t in missing: bad[t*hop:min(n,t*hop+span)]=True
    i=0
    while i<n:
        if not bad[i]: i+=1; continue
        j=i
        while j<n and bad[j]: j+=1
        L=j-i; a=max(0,i-fade); b=min(n,j+fade); M=b-a
        prev=out[:,max(0,a-L):max(0,a-L)+M] if a-L>=0 else np.zeros((out.shape[0],M))
        nxt=out[:,b-M+L:b+L] if b+L<=n else None
        if nxt is None or nxt.shape[1]!=M: nxt=prev
        if prev.shape[1]!=M: prev=nxt
        x=np.linspace(0,1,M); fill=prev*(1-x)+nxt*x            # passage progressif de l'avant à l'après
        win=np.ones(M)
        k=i-a
        if k>0: win[:k]=np.linspace(0,1,k)
        k2=b-j
        if k2>0: win[-k2:]=np.linspace(1,0,k2)
        out[:,a:b]=out[:,a:b]*(1-win)+fill*win
        i=j
    return out


def _ar_coeffs(x, p):
    """Coefficients de prédiction linéaire (méthode d'autocorrélation, Levinson-Durbin)."""
    x = x - x.mean(); n = len(x)
    r = np.correlate(x, x, 'full')[n - 1:n + p] if n < 4096 else np.array([np.dot(x[:n - k], x[k:]) for k in range(p + 1)])
    r = r.astype(np.float64); r[0] *= 1.0 + 1e-6
    a = np.zeros(p + 1); a[0] = 1.0; e = r[0]
    if e <= 0: return np.zeros(p)
    for k in range(1, p + 1):
        acc = r[k] + np.dot(a[1:k], r[k - 1:0:-1])
        ref = -acc / e
        a[1:k + 1] = a[1:k + 1] + ref * np.concatenate([a[k - 1:0:-1], [1.0]])
        e *= (1.0 - ref * ref)
        if e <= 0: break
    return -a[1:]                              # x[n] ~ somme c[k] x[n-1-k]

def _ar_extend(ctx, L, p):
    """Prolonge le contexte de L échantillons par prédiction linéaire."""
    p = min(p, len(ctx) // 3)
    if p < 4: return np.zeros(L)
    m = ctx.mean(); c = _ar_coeffs(ctx, p); buf = np.concatenate([ctx[-p:] - m, np.zeros(L)])
    for i in range(L):
        buf[p + i] = np.dot(c, buf[p + i - 1::-1][:p])
    return buf[p:] + m

def conceal3(audio, missing, hop=512, span=768, order=160, context=2048):
    """Masquage par prédiction linéaire : la zone touchée par un bloc manquant (768 échantillons) est recalculée
    en prolongeant le son qui précède vers l'avant et le son qui suit vers l'arrière, puis en fondant les deux.
    Contrairement à la copie de forme d'onde, la hauteur et le timbre sont continus à travers le trou."""
    out = audio.astype(np.float64).copy(); n = out.shape[1]
    bad = np.zeros(n, bool)
    for t in missing: bad[t * hop:min(n, t * hop + span)] = True
    i = 0
    while i < n:
        if not bad[i]: i += 1; continue
        j = i
        while j < n and bad[j]: j += 1
        L = j - i
        a0 = i
        while a0 > 0 and i - a0 < context and not bad[a0 - 1]: a0 -= 1
        b0 = j
        while b0 < n and b0 - j < context and not bad[b0]: b0 += 1
        w = 0.5 - 0.5 * np.cos(np.pi * (np.arange(L) + 0.5) / L)          # 0 -> 1
        for ch in range(out.shape[0]):
            pre = out[ch, a0:i]; post = out[ch, j:b0]
            fwd = _ar_extend(pre, L, order) if len(pre) >= 64 else None
            bwd = _ar_extend(post[::-1], L, order)[::-1] if len(post) >= 64 else None
            if fwd is None and bwd is None: out[ch, i:j] = 0
            elif bwd is None: out[ch, i:j] = fwd * (1 - w)
            elif fwd is None: out[ch, i:j] = bwd * w
            else: out[ch, i:j] = fwd * (1 - w) + bwd * w
        i = j
    return out


def conceal_tdac(spec_list, missing, synthesize, KC, KS, KBD, W0, order=320, passes=2):
    """Masquage qui respecte la transformée : le bloc manquant est d'abord estimé dans le temps par prédiction
    linéaire, puis ses coefficients MDCT / MDST sont recalculés à partir de cette estimation et remis dans la
    suite des blocs. La synthèse normale fait alors le reste : aux deux bords de la zone, les vrais blocs voisins
    retrouvent leur partenaire de recouvrement, et seule la partie centrale (256 échantillons) est inventée."""
    M = KC.shape[1]; N = KC.shape[0]
    specs = [np.zeros_like(spec_list[0]) if (s is None) else s.copy() for s in spec_list]
    miss = sorted(missing)
    for t in miss: specs[t] = np.zeros_like(spec_list[0] if spec_list[0] is not None else specs[0])
    alpha = 2.0 / M                                  # analyse = alpha x transposée de la synthèse (vérifié numériquement)
    est = conceal3(synthesize(specs), miss, order=order)
    ms = set(miss)
    for it in range(passes):
        if it > 0:
            # affinage : les bords de la zone sont maintenant en partie exacts ; on ne re-prédit que le centre
            # (256 échantillons), depuis des contextes qui le touchent
            for t in miss:
                if (t - 1) in ms or (t + 1) in ms: continue
                i0, j0 = t * 512 + 256, t * 512 + 512; L = j0 - i0
                w = 0.5 - 0.5 * np.cos(np.pi * (np.arange(L) + 0.5) / L)
                for ch in range(est.shape[0]):
                    pre = est[ch, max(0, i0 - 2048):i0]; post = est[ch, j0:j0 + 2048]
                    if len(pre) < 64 or len(post) < 64: continue
                    est[ch, i0:j0] = _ar_extend(pre, L, order) * (1 - w) + _ar_extend(post[::-1], L, order)[::-1] * w
        for t in miss:
            for ch in range(est.shape[0]):
                a = est[ch, t * 512:t * 512 + N]; b = est[ch, t * 512 + M:t * 512 + M + N]
                if len(a) < N or len(b) < N: continue
                specs[t][ch, 0] = alpha * (KC.T @ (a * KBD)); specs[t][ch, 1] = alpha * (KS.T @ (b * KBD))
        est = synthesize(specs)
    return est
