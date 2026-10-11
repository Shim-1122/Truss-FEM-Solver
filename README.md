# Truss FEM Solver

![Tests](https://github.com/Shim-1122/Truss-FEM-Solver/actions/workflows/tests.yml/badge.svg)

Linear-elastic 2D truss solver using the direct stiffness method. Computes nodal displacements, member forces and stresses, and support reactions, including for statically indeterminate trusses. Trusses are defined in JSON; results are printed and plotted with members colored by axial force.

![Pratt Bridge](figures/pratt_bridge.png)

## Validation

| Case | Analytical | Solver |
|---|---|---|
| Pratt bridge reactions | 250 kN / 2 = 125 kN | 125.00 kN |
| Bottom chord 2-3 (method of sections, moments about node 8) | (125·8 − 50·4) / 4 = 200 kN | 200.00 kN |
| Top chord 8-9 (method of sections, moments about node 3) | −(125·12 − 50·4 − 50·8) / 4 = −225 kN | −225.00 kN |
| Center vertical 3-9 | zero-force member | 0.00 kN |
| Three-bar indeterminate truss, vertical bar | P / (1 + 2cos³α) | match, rel. error < 1e-6 |
| Single bar elongation | PL / EA | match, rel. error < 1e-6 |
| Unsupported truss stiffness matrix | 3 rigid-body modes | 3 zero eigenvalues |

## Indeterminate Example

![X-Braced Cantilever](figures/cantilever_x_braced.png)

Wall-mounted cantilever, 4 panels of 1 x 1 m, X-braced, 10 kN tip load. Degree of indeterminacy is 4 (m + r − 2j = 20 + 4 − 20), so equilibrium alone cannot solve it. The stiffness method adds compatibility of displacements at the joints.

| Output | Value |
|---|---|
| Horizontal reactions | ±40 kN (couple = 10 kN x 4 m over 1 m depth) |
| Vertical reactions | 5 kN at each support |
| Tip deflection | 4.86 mm |
| Peak stress | 70 MPa |

## Method

For a member with axial stiffness k = EA/L and direction cosines (c, s):

```
g = [-c, -s, c, s]        elongation = g · u_e        k_e = k · g gᵀ  (4x4)
```

where u_e = [u_i, v_i, u_j, v_j] are the member end displacements.

1. **Assemble:** all element matrices are built in one vectorized array of shape (n_elements, 4, 4) and passed to a sparse COO constructor, which sums entries at shared (row, col) positions.
2. **Apply supports:** supported DOFs have zero displacement, leaving K_ff · u_f = F_f.
3. **Solve** with sparse LU. A singular K_ff (mechanism) raises an error.
4. **Recover:** N = k · (g · u_e), σ = N / A, reactions = K·u − F at supported DOFs.

The degree of static indeterminacy (m + r − 2j) is reported with each solution.

## Getting Started (Mac or Windows, No Experience Needed)

You do not need to know Python. Each step says exactly what to click or type; copying and pasting is fine.

**1. Install Python (once per computer)**

- **Mac:** open Terminal (press Cmd+Space, type `Terminal`, press Return). Type `python3 --version` and press Return.
  - If it shows Python 3.9 or newer, go to step 2.
  - If a box offers to install the "command line developer tools", click **Install**, wait for it to finish, and check again.
  - Otherwise download the macOS installer from https://www.python.org/downloads/, open it, and click through.
- **Windows:** download the installer from https://www.python.org/downloads/ and open it. On the first screen tick **Add python.exe to PATH**, then click **Install Now**.

**2. Download this project**

Click the green **Code** button at the top of this page, then **Download ZIP**. Unzip it: on a Mac double-click the ZIP in Downloads; on Windows right-click it, choose **Extract All**, then **Extract**. You now have a folder named `Truss-FEM-Solver-main`.

**3. Start it**

- **Windows:** open the `Truss-FEM-Solver-main` folder and double-click **start.bat**.
- **Mac:** in Terminal type `cd ` (with a space after it), drag the `Truss-FEM-Solver-main` folder from Finder into the Terminal window, and press Return. Then type `python3 start.py` and press Return.

The first time, it sets itself up: it creates a private Python environment inside the project folder and downloads the libraries it needs. This takes a few minutes and needs internet. Then a menu appears. Type a number and press Enter:

| Option | What it does | What you get |
| --- | --- | --- |
| 1 | Solve the example Pratt bridge | Prints every member force, the reactions and the largest deflection, saves `figures/pratt_bridge.png` and opens it |
| 2 | Solve the example X-braced cantilever | Same for the statically indeterminate cantilever |
| 3 | Solve your own truss file | Asks for the path to a `.json` file written as described in Input Format below (drag the file into the window to paste its path) |
| 4 | Run the automated tests | Checks the solver against hand-calculated solutions |

To use it again later, repeat step 3; setup is skipped and the menu appears in a few seconds. Type `q` and press Enter to quit.

**Troubleshooting**

| Problem | Fix |
| --- | --- |
| Mac: `command not found: python3` | Install Python (step 1). |
| Windows: `Python was not found`, or the Microsoft Store opens | Install Python from python.org (step 1) and tick **Add python.exe to PATH**. If you already installed it without that box, run the installer again, choose **Modify**, then **Next**, and tick **Add Python to environment variables**. |
| Windows: "Windows protected your PC" when opening start.bat | Click **More info**, then **Run anyway**. The file only runs `py start.py`; you can open it in Notepad to check. |
| `This project needs Python 3.9 or newer` | Install the current Python from https://www.python.org/downloads/ |
| `FileNotFoundError` for your own truss | Check the path. Dragging the file into the window types the full path for you. |
| `Installing libraries failed` or setup stops with an error | Check your internet connection, delete the `.venv` folder inside the project folder, and start again. (On a Mac the folder is hidden: press Cmd+Shift+. in Finder to show it.) |
| A figure did not open by itself | Look in the project's `figures` folder. |

**For developers.** `start.py` only creates `.venv`, installs `requirements.txt` and runs the commands below, so you can also work by hand (Python 3.9 or newer):

```bash
# macOS / Linux
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

```powershell
# Windows (PowerShell or Command Prompt)
py -m venv .venv
.venv\Scripts\activate       # if PowerShell blocks this: Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
pip install -r requirements.txt
```

```bash
python solve.py examples/pratt_bridge.json
python -m pytest
```

`python start.py <shortcut>` runs one menu option directly (`python start.py --list` shows the shortcuts). CI runs the tests through `start.py` on macOS, Windows and Linux with Python 3.9 and 3.13.

### Command Line

Output (abridged):

```
12 nodes, 21 members, statically determinate

member    nodes  force (kN)  stress (MPa)  state
     2      2-3      200.00         40.00  tension
     7      8-9     -225.00        -45.00  compression
    14      3-9        0.00          0.00  zero-force
...
Max displacement: 9.714 mm at node 3
```

Options: `--scale` sets deformation magnification (default: auto), `--out` sets the figure path (default: `figures/<name>.png`).

## Input Format

```json
{
  "name": "Example Truss",
  "E": 200e9,
  "A": 0.005,
  "nodes":    [[0, 0], [4, 0], [2, 3]],
  "elements": [[0, 1], [1, 2], [2, 0]],
  "supports": {"0": "pin", "1": "roller"},
  "loads":    {"2": [0, -10000]}
}
```

SI units (m, N, Pa). `E` and `A` accept one value or one value per member. Support types: `"pin"` (x and y fixed), `"roller"` (y fixed), `"roller_x"` (x fixed), or `[fix_x, fix_y]`.

## Project Structure

```
start.py                    Beginner launcher: sets up .venv, shows the menu (start.bat on Windows)
truss/
  model.py     truss definition, JSON loading, indeterminacy
  solver.py    stiffness assembly, solve, force recovery
  plot.py      deformed shape colored by axial force
solve.py       command-line interface
examples/      Pratt bridge (determinate), X-braced cantilever (indeterminate)
tests/
```

## Tests

- Single bar: δ = PL/EA
- Member force is invariant to rotating the whole problem (5 orientations)
- Two-bar hanger and Pratt bridge forces match statics
- Three-bar indeterminate truss matches the closed-form compatibility solution
- Indeterminate cantilever satisfies ΣF = 0 and ΣM = 0
- Unsupported stiffness matrix is symmetric with exactly 3 rigid-body modes
- Mechanisms raise a clear error

## License

MIT
