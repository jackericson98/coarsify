import tkinter as tk
from tkinter import filedialog
from System.system import System


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    scheme = input("Choose coarsify scheme: (1) Average Distance, (2) Encapsulate Residues, (3) CG Martini >>>   ")
    root = tk.Tk()
    root.withdraw()
    file = filedialog.askopenfilename()
    sys = System(file=file, scheme=scheme)
