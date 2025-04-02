# ResGA

**ResGA** is a C++/Python-based genetic algorithm framework for fitting resonance parameters in hadronic transport models, developed alongside the study:

📄 **[Systematic optimization of resonance parameters in a transport approach](https://arxiv.org/abs/2503.05504)**  
_Carl B. Rosenkvist and Hannah Elfner, 2025_

The algorithm is designed to optimize properties such as masses, widths, and branching ratios by comparing transport model output (e.g., from SMASH) to experimental data on exclusive cross-sections.

It supports modular simulation wrapping and flexible scoring strategies.

---

## ⚠️ Requirements and Customization

This repository is **not expected to run out of the box**, since it requires:

- **Experimental data** for target reactions
- **Simulation output** (e.g., SMASH runs) to be evaluated for fitness

The `emulator.py` script is a user-defined interface between the GA and your model. You can **replace or adapt `emulator.py`** to fit your needs — whether you're working with different observables, models, or channels.

Future updates may include minimal examples or templates to help you get started.
