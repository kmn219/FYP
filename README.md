# Observability of Modes in Power Systems Through Signal Injection

This repository contains the report, presentation, and code used for a final year project submission as part of the MEng degree in Electrical and Electronic Engineering at Imperial College London.

## Abstract

Due to rising pressures of climate change, modern power grids are relying more and more on renewable methods of energy generation. Many of these are inverter-based
resources (IBRs), such as wind turbines (e.g. doubly fed induction generators) and solar panels, that interface with the grid using power electronics. This results
in grid stability dynamics being increasingly influenced by controller profiles which are often not fully reflected in simulation models, leading to unexpected
oscillatory modes showing up under normal operation.

This project explores the methods of online impedance measurement for stability analysis through small-signal injections. The method is first contextualised based
on literature. It is then validated in simulation on a basic single phase model, showing that measurements can be taken only at one node to get a full picture of
the system with no noise present. The impacts of noise are discussed. Finally, a signal injection device topology based on small- signal current and voltage injection
is proposed, exploring the possible implementation using available technologies.

## Repository Structure

The code folder contains the MATLAB scripts and Simulink models used for this project, including the [Simplus Grid Tool](https://github.com/Future-Power-Networks/Simplus-Grid-Tool) developed as part of a PhD thesis by Yue Zhu in the EEE Department at Imperial College London. You can find the interim report, outlining the intial project goals in the project management folder.
