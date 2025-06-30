# ase\_asro\_calculator

`ase_asro_calculator` - some routines for calculating Warren-Cowley atomic short-range order (ASRO) parameters for simulation cells described by an [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase//index.html) Atoms object.

## Citation

There is no publication associated with these (very lightweight) scripts. However, if they were useful for your research, you may wish to consider citing my paper, [Phys. Rev. Mater **7**, 103801 (2023)](https://doi.org/10.1103/PhysRevMaterials.7.013801), where a definition of the Warren-Cowley parameters in the multicomponent setting is given (Eq. 12), along with appropriate citations to the original literature where they were first defined. And, of course, if you are using ASE, you should include a citation to the [relevant paper for that](https://doi.org/10.1088/1361-648X/aa680e), too.

## Generating special quasirandom structures (SQSs)

In the `structures` directory, I have a few disordered SQS supercells generated using `icet` ([https://icet.materialsmodeling.org/index.html](https://icet.materialsmodeling.org/index.html)). They have a [nice tutorial on the generation of SQS supercells](https://icet.materialsmodeling.org/advanced_topics/sqs_generation.html). If you use this, you should give them a citation, too! (Details on their website.)

Copyright (C) C. D. Woodgate 2025. Released under the GNU Lesser General Public License, version 3.
