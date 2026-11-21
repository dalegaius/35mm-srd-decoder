#!/usr/bin/env python3
"""One-click front end: pick a folder of SR-D images, get WAV files (EC9 and EC11).

What it runs depends on what the folder holds:
  - camera RAW photos (.dng): srd_prepare.py (reads the blocks) then srd_photos.py (corrects them, combines the
    reads of a block seen in several photos, decodes). This is the path that gives the most blocks;
  - a folder that already holds points.npz and rapport.csv (an earlier srd_prepare output): srd_photos.py only;
  - other images (JPEG, PNG, TIFF) or .npz block files from srd_extract_video.py: srd_decode.py.
Results go to <folder>/srd_decoded/.
"""
import os, re, sys, time, threading, subprocess, queue
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = {'.jpg', '.jpeg', '.png', '.tif', '.tiff', '.bmp'}


def python_exe():
    """a console python even when the GUI itself runs under pythonw (so the child scripts can print)"""
    exe = sys.executable
    if os.path.basename(exe).lower() == 'pythonw.exe':
        cand = os.path.join(os.path.dirname(exe), 'python.exe')
        if os.path.exists(cand): return cand
    return exe


def plan(folder, name, workers, fast):
    """-> (list of (title, command), output prefix, error message or None)"""
    out = os.path.join(folder, 'srd_decoded'); prefix = os.path.join(out, name)
    files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
    ext = lambda e: [f for f in files if os.path.splitext(f)[1].lower() in e]
    py = [python_exe(), '-u']
    photos = [py + [os.path.join(HERE, 'srd_photos.py')], ['-o', prefix, '--workers', str(workers)] + (['--fast'] if fast else [])]
    if 'points.npz' in files and 'rapport.csv' in files:
        return [('Correcting and decoding the blocks', photos[0] + [folder] + photos[1])], prefix, None
    if ext({'.dng'}):
        steps = []
        if os.path.exists(os.path.join(out, 'points.npz')) and os.path.exists(os.path.join(out, 'rapport.csv')):
            steps.append(('Photos already read (srd_decoded/points.npz found): delete that file to read them again', None))
        else:
            steps.append((f'Reading the blocks of {len(ext({".dng"}))} photos (about 8 s per photo and per process)',
                          py + [os.path.join(HERE, 'srd_prepare.py'), folder, '-o', out, '--no-crops', '--workers', str(max(1, min(workers, (os.cpu_count() or 2) // 2, 16)))]))
        steps.append(('Correcting and decoding the blocks', photos[0] + [out] + photos[1]))
        return steps, prefix, None
    inputs = [os.path.join(folder, f) for f in sorted(ext({'.npz'}))] or ([folder] if ext(IMG) else [])
    if not inputs:
        return [], prefix, 'No DNG photos, images, or block files (.npz) in this folder.'
    return [('Reading and decoding', py + [os.path.join(HERE, 'srd_decode.py')] + inputs + ['-o', prefix, '--workers', str(workers)])], prefix, None


PROGRESS = re.compile(r'(\d+)\s*/\s*(\d+)')


def fmt(sec):
    sec = int(sec); return f'{sec // 60} min {sec % 60:02d} s' if sec >= 60 else f'{sec} s'


def run_steps(steps, log, status=lambda *a: None):
    """runs the commands one after the other, sending their output to log(). status(step title, start time,
    done, total) is called whenever a line of output carries a 'done/total' count. -> True if all succeeded"""
    env = dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8')
    flags = getattr(subprocess, 'CREATE_NO_WINDOW', 0); k = 0; real = sum(1 for _, c in steps if c)
    for title, cmd in steps:
        log('== ' + title)
        if cmd is None: continue
        k += 1; t0 = tp = time.time(); last_total = 0; label = f'Step {k} of {real}: {title.split(" (")[0]}'; status(label, t0, 0, 0, tp)
        p = subprocess.Popen(cmd, cwd=HERE, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                             encoding='utf-8', errors='replace', creationflags=flags)
        for line in p.stdout:
            line = line.rstrip(); log(line)
            m = PROGRESS.search(line)
            if m and line.lstrip().startswith(('[', m.group(1))):      # "[12/300] ..." or "500/8000 blocs ..."
                if int(m.group(2)) != last_total: last_total = int(m.group(2))      # a new count: time it from the previous phase change
                status(label, t0, int(m.group(1)), int(m.group(2)), tp)
            elif line and not line.startswith(' '):
                tp = time.time(); last_total = 0; status(label, t0, 0, 0, tp)       # a new phase inside the step
        if p.wait() != 0:
            log(f'This step stopped with an error (code {p.returncode}).'); return False
        log(f'   step finished in {fmt(time.time() - t0)}')
    return True


def open_folder(path):
    if sys.platform.startswith('win'): os.startfile(path)
    elif sys.platform == 'darwin': subprocess.Popen(['open', path])
    else: subprocess.Popen(['xdg-open', path])


def status_text(label, t0, done, total, tp=None, now=None):
    """elapsed time of the step; the time left is estimated for the current phase only (tp = its start)"""
    now = now or time.time(); el = now - t0; ph = now - (tp if tp is not None else t0); txt = f'{label}   |   elapsed {fmt(el)}'
    if total and done:
        txt += f'   |   {done}/{total} ({100 * done // total} %)'
        if 0 < done < total and ph > 3: txt += f'   |   about {fmt(ph * (total - done) / done)} left in this phase'
    return txt


class App:
    def __init__(self, root):
        self.root = root; root.title('SR-D decoder'); root.geometry('820x560')
        tk.Label(root, text='Dolby SR-D decoder', font=('Helvetica', 16, 'bold')).pack(pady=(14, 2))
        tk.Label(root, text='Choose a folder of photos of the SR-D blocks (DNG), of images, or of block files (.npz).\n'
                            'Revisions 9 (EC9) and 11 (EC11) are both decoded.').pack()
        row = tk.Frame(root); row.pack(pady=8)
        self.btn = tk.Button(row, text='Choose a folder and decode...', font=('Helvetica', 13), command=self.pick)
        self.btn.pack(side='left', padx=8)
        tk.Label(row, text='Processes:').pack(side='left')
        self.workers = tk.IntVar(value=max(1, (os.cpu_count() or 2) - 2))
        tk.Spinbox(row, from_=1, to=64, width=4, textvariable=self.workers).pack(side='left', padx=(2, 10))
        self.fast = tk.BooleanVar(value=False)
        tk.Checkbutton(row, text='Fast (a few blocks fewer)', variable=self.fast).pack(side='left')
        self.status = tk.Label(root, text='', font=('Helvetica', 11, 'bold'), anchor='w'); self.status.pack(fill='x', padx=12)
        self.log = scrolledtext.ScrolledText(root, height=20, state='disabled'); self.log.pack(fill='both', expand=True, padx=12, pady=(0, 12))
        self.q = queue.Queue(); self.prog = None; root.after(200, self.poll)
        if not os.path.exists(os.path.join(HERE, 'codebooks.npz')):
            self.write('Note: codebooks.npz is not in the decoder folder. EC9 prints decode without it; EC11 prints need it '
                       '(see NOTICE.md).')

    def write(self, msg):
        self.log.configure(state='normal'); self.log.insert('end', msg + '\n'); self.log.see('end'); self.log.configure(state='disabled')

    def poll(self):
        while not self.q.empty():
            kind, payload = self.q.get()
            if kind == 'log': self.write(payload)
            elif kind == 'status': self.prog = payload
            elif kind == 'done':
                self.prog = None; self.status.configure(text=f'Finished in {fmt(time.time() - self.t_start)}.')
                self.btn.configure(state='normal'); self.write('Done.'); open_folder(payload)
            elif kind == 'error':
                self.prog = None; self.status.configure(text='Stopped.')
                self.btn.configure(state='normal'); self.write('Error: ' + payload); messagebox.showerror('SR-D decoder', payload)
        if self.prog: self.status.configure(text=status_text(*self.prog))
        self.root.after(200, self.poll)

    def pick(self):
        folder = filedialog.askdirectory(title='Folder of SR-D images')
        if not folder: return
        folder = os.path.normpath(folder); name = os.path.basename(folder) or 'srd'
        steps, prefix, err = plan(folder, name, int(self.workers.get()), bool(self.fast.get()))
        self.write(f'Folder: {folder}')
        if err: self.write(err); return
        os.makedirs(os.path.dirname(prefix), exist_ok=True)
        self.btn.configure(state='disabled'); self.t_start = time.time()
        threading.Thread(target=self.work, args=(steps, os.path.dirname(prefix)), daemon=True).start()

    def work(self, steps, out):
        try:
            ok = run_steps(steps, lambda m: self.q.put(('log', m)), lambda *s: self.q.put(('status', s)))
            self.q.put(('done', out) if ok else ('error', 'Decoding did not finish. The messages above say why.'))
        except Exception as e:
            self.q.put(('error', repr(e)))


if __name__ == '__main__':
    root = tk.Tk(); App(root); root.mainloop()
