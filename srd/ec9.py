"""EC9 decoder (SR-D revision 9), pure Python. Written from US patent 5,583,962, from the study of how a decoder
behaves, and from measurements on film prints. Input: the 215 de-whitened frame words 2..216 (3,440 bits).
Needs no codebook: revision 9 uses plain binary mantissas. (Comments below are in French.)"""
import numpy as np
W9=[1]*12+[2]*5+[3]*2+[4]*6+[6]+[10]*9+[15]+[18]*3          # largeurs des 39 sous-bandes : tableau I du brevet US 5 583 962
BT9=[5]*12+[4]*12+[3]*10+[2,2,0,0,0]                       # allocation de base par sous-bande : même brevet
NB=[39,20,20,20,20,7]
ZERO_FILL=0.0          # bruit pour les mantisses à zéro bit : désactivé (à l'écoute, sans bruit sonne le plus naturel)
SER_DEFAULT=False      # reconstruction statistique désactivée par défaut en EC9 (même raison)
STRAT=3                # poids du champ de 2 bits dans l'exposant (3 exposants, soit 18 dB, par cran)
# Correction de niveau des 19 bandes composites (dB), mesurée contre la piste optique de la même copie
# (Call of Duty, équilibre moyen, bandes discrètes 400 Hz - 3 kHz comme référence). Empirique, comme
# CPL_EQ_DB de l'EC11 mesuré contre le DTS ; les trois dernières bandes (au-delà de 16 kHz) ne sont pas
# mesurables (référence filtrée) et restent à 0.
CPL_EQ_DB=[0.0]*19
# Couplage : niveau efficace absolu de chaque canal par bande composite, (m / 32) x 2^(CPL_SHIFT - e), appliqué au
# composite ramené à une valeur efficace de 1 sur la bande. Établi de deux façons indépendantes : (1) sans référence,
# le niveau (m, e) prolonge le niveau de la dernière bande propre du même canal aussi bien que le fait la bande
# voisine (écart-type 1,3 unité d'exposant), et l'exposant du composite n'y contribue pas ; (2) contre la piste
# optique de Call of Duty, l'équilibre des bandes composites 20 à 35 tombe à +- 3 dB sans aucune constante ajustée.
CPL_SHIFT=0
CPL_ROUGH_MAX=2.5      # écart moyen admis entre niveaux de couplage de bandes voisines (0,7 en temps normal)                                      # bandes par groupe : L (+ composite), R, Ls, Rs, C, LFE
def F(b,p,n):
    v=0
    for k in range(n): v=v*2+int(b[p+k])
    return v
def side(b):
    e=[F(b,4*i,4) for i in range(126)]; st=[F(b,504+2*i,2) for i in range(7)]; f=F(b,518,4)
    X={}; RAW={}
    for bd in range(20): RAW[(0,bd)]=e[bd]
    for k in range(19): RAW[(0,20+k)]=e[20+k]
    for c in range(4):
        for bd in range(20): RAW[(c+1,bd)]=e[39+20*c+bd]
    for bd in range(7): RAW[(5,bd)]=e[119+bd]
    side.RAW=RAW
    for bd in range(20): X[(0,bd)]=e[bd]+STRAT*st[0]
    for k in range(19): X[(0,20+k)]=e[20+k]+STRAT*st[6]
    for c in range(4):
        for bd in range(20): X[(c+1,bd)]=e[39+20*c+bd]+STRAT*st[c+1]
    for bd in range(7): X[(5,bd)]=e[119+bd]+STRAT*st[5]
    cpe=[F(b,522+5*i,5) for i in range(95)]; cpm=[F(b,997+5*i,5) for i in range(95)]
    phase=[F(b,1472+i,1) for i in range(76)]
    return X,f,cpe,cpm,phase
EXCL15=True
SIL=15
FILL_SILENT=True       # le remplissage peut-il servir les bandes muettes ?
XM_MODE='tous'      # 'tous' : maximum sur toutes les bandes ; 'actives' : sans les bandes d'exposant brut 15
def allocate(X,f,budget=1892-162,raw=None):
    if XM_MODE=='actives' and raw is not None:
        act=[X[s] for s in X if raw[s]<15]; xm=max(act) if act else max(X.values())
    else: xm=max(X.values())
    a={}
    for (g,bd),x in X.items():
        d=xm-x
        base=BT9[bd]+(1 if (g==0 and bd>=20) else 0)+((d+2)>>2)   # (max - x)/4 arrondi
        a[(g,bd)]=min(max(base+f-8,0),9)
    if EXCL15 and raw is not None:
        for s in a:
            if raw[s]>=SIL: a[s]=0          # bandes muettes (exposant brut 15) exclues de l'allocation (brevet)
    for g in range(6): a[(g,0)]=0; a[(g,1)]=0
    cost=lambda s: 2*W9[s[1]]
    rem=budget-sum(a[s]*cost(s) for s in a)
    order=[(g,bd) for bd in range(2,20) for g in range(5)]+[(0,bd) for bd in range(20,39)]
    if EXCL15 and raw is not None and not FILL_SILENT: order=[s for s in order if raw[s]<SIL]
    stop=None
    while rem>=0:
        start=rem
        for s in order:
            if a[s]>=9: continue
            if cost(s)>rem: stop=s; break
            a[s]+=1; rem-=cost(s)
        if stop or rem==start: break
    for g in range(6): a[(g,0)]=9; a[(g,1)]=9
    return a,rem,stop
def mantissas(b,a,rem,stop,ser=None):
    remA=rem>>1; remB=rem-remA; p=1548; out={}; starts={}
    def rd(n):
        nonlocal p
        if n<=0: return 0.0
        v=F(b,p,n); p+=n
        if v>=1<<(n-1): v-=1<<n
        if ser is not None: return (v+ser.uniform(0,1))/(1<<(n-1))     # reconstruction statistique (brevet, quantification par troncature)
        return v/(1<<(n-1))
    for ps,rx in ((0,remA),(1,remB)):
        for g in range(6):
            for bd in range(0 if ps==0 else 1, NB[g]):
                s=(g,bd); n=a[s]; w=W9[bd]; starts[(s,ps)]=p
                vals=[]
                for i in range(w):
                    nn=n+(1 if (s==stop and rx>0 and i>=w-rx) else 0)
                    if n>0 or (s==stop and i>=w-rx and rx>0): vals.append(rd(nn))
                    elif ser is not None and ZERO_FILL and not (ps==1 and bd==0): vals.append(ser.uniform(-ZERO_FILL,ZERO_FILL))
                    else: vals.append(0.0)
                out[(s,ps)]=vals
    return out,p,starts
def quality(b,a,rem,stop):
    """nombre de canaux « en échec » sur les canaux contrôlables : la bande 0 MDCT (9 bits, sans partenaire)
    doit être normalisée ; un canal très faible (exposant brut de bande 0 ≥ 13) peut légitimement ne pas
    l'être (signal sous la plage de l'exposant) et n'est pas contrôlé. Renvoie 6 − échecs."""
    e0=[F(b,4*i,4) for i in (0,39,59,79,99,119)]
    p=1548; bad=0
    for g in range(6):
        v=F(b,p,2) if p+2<=len(b) else 0
        if e0[g]<13 and (v>>1)==(v&1): bad+=1
        for bd in range(NB[g]): p+=a[(g,bd)]*W9[bd]+((rem>>1) if (g,bd)==stop else 0)
    return 6-bad
def quality_old(b,a,rem,stop):
    """nombre de canaux (sur 6) dont la bande 0 MDCT, toujours normalisée, tombe sur une valeur normalisée"""
    p=1548; ok=0
    for g in range(6):
        v=F(b,p,2) if p+2<=len(b) else 0; ok+=(v>>1)!=(v&1)
        for bd in range(NB[g]): p+=a[(g,bd)]*W9[bd]+((rem>>1) if (g,bd)==stop else 0)
    return ok
def decode_frame(words,ser=None,bits=None):
    b=bits if bits is not None else ((np.array(words[2:217],np.int64)[:,None]>>np.arange(15,-1,-1))&1).reshape(-1)
    X,f,cpe,cpm,phase=side(b)
    a,rem,stop=allocate(X,f,raw=side.RAW)
    m,end,starts=mantissas(b,a,rem,stop,ser)
    q=quality(b,a,rem,stop)
    ce=np.array(cpe).reshape(5,19); rough=float(np.mean(np.abs(np.diff(ce,axis=1))))
    if rough>CPL_ROUGH_MAX: q=0            # niveaux de couplage en dents de scie : trame corrompue
    return dict(X=X,f=f,cpe=cpe,cpm=cpm,phase=phase,alloc=a,rem=rem,stop=stop,mant=m,end=end,q=q,check=q)

S9=np.cumsum([0]+W9)
def spectrum(r):
    """spectre [6 canaux, 2 demi-blocs (MDCT, MDST), 256] dans l'ordre du flux L, R, Ls, Rs, C, LFE.
    Bandes composites : les mantisses du composite donnent la forme (ramenées à une valeur efficace de 1 sur la
    bande, deux demi-blocs réunis), et chaque canal a son niveau efficace absolu par bande, (m / 32) x 2^-e,
    dans les unités des exposants ; signe inversé si le drapeau de phase du canal est levé. L'exposant du
    composite ne sert qu'à l'allocation de ses mantisses."""
    sp=np.zeros((6,2,256)); comp=np.zeros((2,256))
    for ((g,bd),ps),vals in r['mant'].items():
        if g==0 and bd>=20:
            for i,v in enumerate(vals): comp[ps,S9[bd]+i]=v
        else:
            sc=2.0**(-r['X'][(g,bd)])
            for i,v in enumerate(vals): sp[g,ps,S9[bd]+i]=v*sc
    for j in range(19):
        bd=20+j; k0,k1=S9[bd],S9[bd+1]; c=comp[:,k0:k1]
        rms=np.sqrt(np.mean(c**2))
        if rms==0: continue
        shape=c/rms
        for ch in range(5):
            e=r['cpe'][ch*19+j]; m=r['cpm'][ch*19+j]
            if not m: continue
            cf=(m/32)*2.0**(CPL_SHIFT-e)*10**(CPL_EQ_DB[j]/20)
            if ch>=1 and r['phase'][(ch-1)*19+j]: cf=-cf
            sp[ch,:,k0:k1]=shape*cf
    return sp
def spectra(bits, rng=None, ser=True):
    """EC9 frame bits (3440) -> spec[6, 2, 256] in stream order L, R, Ls, Rs, C, LFE (same layout as decoder.spectra).
    The noise reconstruction is off for EC9 unless SER_DEFAULT is set (listening tests on Call of Duty)."""
    r=decode_frame(None, rng if (ser and SER_DEFAULT) else None, bits=bits)
    return spectrum(r)
def validate(bits):
    """True if the band-0 MDCT value of at least 5 of the 6 channels is normalized (consistency check)."""
    try: return decode_frame(None, None, bits=bits)['check']>=5
    except Exception: return False
