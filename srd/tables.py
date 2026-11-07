import os, numpy as np
_T = None
def load(path=None):
    global _T
    if _T is None:
        path = path or os.path.join(os.path.dirname(__file__), '..', 'tables.npz')
        if not os.path.exists(path):
            raise SystemExit(
                'tables.npz is missing.\n\n'
                'This repository does not ship Dolby firmware data. Generate it once from a dump of\n'
                'hardware you own (the 84147 flash on a Cat. No. 673A card, 128 KB):\n\n'
                '    python tools/extract_tables.py your_84147_dump.bin\n\n'
                'See NOTICE.md for why it works this way.')
        z = np.load(path); _T = {k: z[k] for k in z.files}
    return _T
def H(a):
    t = load(); return int(t['hq'][a - int(t['hq_base'])])
