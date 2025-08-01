# Buzcode Guide for New Users

Welcome to the Buzsaki Lab's buzcode repository. This guide expands on the basic instructions in `readme.txt` and is intended for newcomers to MATLAB and this codebase. After completing this guide you should feel comfortable navigating and extending the repository.

## 1. Initial Setup

1. **Clone or download** the repository and open MATLAB in the repository root.
2. **Compile the MEX files** by running:

   ```matlab
   compileBuzcode
   ```

   This command recompiles C/MEX sources and pulls required external packages into your MATLAB path. You must run it whenever the repository or your platform changes.

3. **Explore the tutorials** in the `tutorials/` folder. They contain example data and scripts demonstrating typical workflows.

4. **Consult the wiki** for additional documentation, including [data formatting standards](https://github.com/buzsakilab/buzcode/wiki/Data-Formatting-Standards).

## 2. Repository Structure

- `analysis/` – core analysis functions (e.g., ripple detection, place fields).
- `detectors/` – event detectors organized by output type.
- `preprocessing/` – tools for spike sorting, LFP extraction, and behavioral tracking.
- `io/` – utilities to load and save data from disk.
- `utilities/` – general helper functions (filtering, plotting, etc.).
- `database/` – scripts for interfacing with the Buzsaki Lab database.
- `tutorials/` – example scripts and datasets for learning workflows.
- `externalPackages/` – third-party libraries used by buzcode.

Use these directories to locate functions or contribute new code in the appropriate place.

## 3. Typical Workflow

1. **Load data** using I/O functions such as `bz_GetLFP`, `bz_GetSpikes`, and `bz_LoadBehavior`.
2. **Run preprocessing** (if needed) to prepare your data. For example, `preprocessing/clustering` scripts handle spike sorting.
3. **Apply detectors or analyses** from the `analysis/` directory. Many functions output structures with standardized field names (e.g., `.events` or `.states`).
4. **Visualize results** using tools in `visualization/` or plotting examples in the tutorials.
5. **Save results** in the recommended format (`.events.mat`, `.states.mat`, etc.) following the data formatting standards.

## 4. Contributing and Getting Help

- If you encounter bugs or missing features, open an issue on GitHub.
- Fork the repository and submit pull requests for improvements.
- Check the `docs/CONTRIBUTING.md` for guidelines on reporting issues and contributing code.
- Engage with other users via the issue tracker to share tips and solutions.

## 5. Learning MATLAB Basics

If you are new to MATLAB, start with these resources:

- MATLAB's official [Getting Started](https://www.mathworks.com/help/matlab/getting-started-with-matlab.html) guide.
- The tutorials in `tutorials/` demonstrate practical examples of loading data, running analysis functions, and visualizing output.

## 6. Tips for Mastery

- Study function headers and comments – they often contain usage examples.
- Explore the `tutorials/exampleDataStructs` folder to see expected data structures.
- Use MATLAB's `help` and `doc` commands to read documentation from within the environment (e.g., `help bz_GetLFP`).
- Experiment with small datasets to understand each function's inputs and outputs before processing large files.

Following these steps will help you become proficient with the buzcode repository. Continue to explore the codebase and contribute improvements as you learn more.

