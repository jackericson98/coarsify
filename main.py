import tkinter as tk
from tkinter import filedialog
from System.system import System


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    root = tk.Tk()
    root.withdraw()
    file = filedialog.askopenfilename()
    sys = System(file=file)