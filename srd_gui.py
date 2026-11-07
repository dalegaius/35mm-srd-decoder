#!/usr/bin/env python3
"""One-click front end: pick a folder of SR-D images, get WAV files next to them."""
import os, sys, threading, subprocess, queue, multiprocessing
import tkinter as tk
from tkinter import filedialog, scrolledtext, messagebox

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def open_folder(path):
    if sys.platform.startswith('win'): os.startfile(path)
    elif sys.platform == 'darwin': subprocess.Popen(['open', path])
    else: subprocess.Popen(['xdg-open', path])

class App:
    def __init__(self, root):
        self.root = root; root.title('SR-D decoder'); root.geometry('720x440')
        tk.Label(root, text='Dolby SR-D decoder', font=('Helvetica', 16, 'bold')).pack(pady=(14, 2))
        tk.Label(root, text='Choose a folder of photos or scans of the SR-D blocks (JPEG, PNG, TIFF or camera RAW).').pack()
        self.btn = tk.Button(root, text='Choose a folder and decode...', font=('Helvetica', 13), command=self.pick)
        self.btn.pack(pady=10)
        self.log = scrolledtext.ScrolledText(root, height=16, state='disabled'); self.log.pack(fill='both', expand=True, padx=12, pady=(0, 12))
        self.q = queue.Queue(); root.after(200, self.poll)

    def write(self, msg):
        self.log.configure(state='normal'); self.log.insert('end', msg + '\n'); self.log.see('end'); self.log.configure(state='disabled')

    def poll(self):
        while not self.q.empty():
            kind, payload = self.q.get()
            if kind == 'log': self.write(payload)
            elif kind == 'done':
                self.btn.configure(state='normal')
                if payload:
                    self.write('Done.'); open_folder(os.path.dirname(payload[0]))
            elif kind == 'error':
                self.btn.configure(state='normal'); self.write('Error: ' + payload); messagebox.showerror('SR-D decoder', payload)
        self.root.after(200, self.poll)

    def pick(self):
        folder = filedialog.askdirectory(title='Folder of SR-D images')
        if not folder: return
        self.btn.configure(state='disabled'); self.write(f'Folder: {folder}')
        out = os.path.join(folder, 'srd_decoded'); os.makedirs(out, exist_ok=True)
        threading.Thread(target=self.work, args=(folder, os.path.join(out, os.path.basename(folder.rstrip('/\\')) or 'srd')), daemon=True).start()

    def work(self, folder, out):
        try:
            import srd_decode
            paths = srd_decode.run([folder], out, log=lambda m: self.q.put(('log', m)))
            self.q.put(('done', paths))
        except SystemExit as e:
            self.q.put(('error', str(e)))
        except Exception as e:
            self.q.put(('error', repr(e)))

if __name__ == '__main__':
    multiprocessing.freeze_support()
    root = tk.Tk(); App(root); root.mainloop()
