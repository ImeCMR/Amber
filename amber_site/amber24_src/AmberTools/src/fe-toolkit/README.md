One can read the FE-ToolKit documentation online at: https://rutgerslbsr.gitlab.io/fe-toolkit

This repository is a collection of programs useful for free energy analysis.
The programs are stored in individual tarballs. Each program has its own
license, which can be viewed by extracting the contents of the program
distribution tarball.

The programs included in this distribution include:


ndfes     - a program for evaluating free energy profiles from umbrella
            window simulations using either vFEP or MBAR analysis.

edgembar  - a program for performing networkwide free energy analysis
            of alchemical free energy transformation graphs typically
	    constructed to compute relative binding (or solvation)
	    free energies.

The programs rely on the presence of external software not included.
Specifically, these programs require an implementation of BLAS, LAPACK,
and NLopt (https://github.com/stevengj/nlopt). The configuration and
installation is performed with cmake (version 3.12) and pip (a recent
version).

If you are using Fedora, you can install the missing dependencies with:
```
   sudo dnf install cmake.x86_64 openblas-serial.x86_64 NLopt.x86_64 \
                    NLopt-devel.x86_64 python3-numpy.x86_64 \
		    python3-scipy.x86_64 python3-matplotlib.x86_64
```

Alternatively, if the BLAS/LAPACK and NLopt libraries are missing,
the cmake build system will download them from github and compile them.
To install cmake and update your version of pip, you can run:
- USERBASE=\$(python3 -m site --user-base)
- USERSITE=\$(python3 -m site --user-site)
- export PATH="\${USERBASE}/bin:\${PATH}"
- export PYTHONPATH="\${USERSITE}:\${PYTHONPATH}"
- python3 -m pip install pip --upgrade --user
- python3 -m pip install cmake --upgrade --user
You should consider adding the above export commands to your
${HOME}/.bashrc and then "source ~/.bashrc".
If, for whatever reason, pip is unavailable on your system, you
can install it using the directions here:
https://pip.pypa.io/en/stable/installation

To install fe-toolkit,
- cd build
- bash ./run_cmake.sh
- make install VERBOSE=1 -j4
- cd ../local
- export PATH="\${PWD}/bin:\${PATH}"
- export PYTHONPATH="\${PWD}/lib/python3.XX/site-packages:\${PYTHONPATH}"
where python3.XX should be replaced by the appropriate python version.


We request that if you use this software in a publication, to please reference
as appropriate:

[1] Alchemical Binding Free Energy Calculations in AMBER20: Advances and Best 
Practices for Drug Discovery
Tai-Sung Lee, Bryce K. Allen, Timothy J. Giese, Zhenyu Guo, Pengfei Li, 
Charles Lin, T. Dwight McGee, David A. Pearlman, Brian K. Radak, Yujun Tao, 
Hsu-Chun Tsai, Huafeng Xu, Woody Sherman, Darrin M. York
J. Chem. Inf. Model. (2020) 60, 5595-5623
DOI: 10.1021/acs.jcim.0c00613

[2] Variational Method for Networkwide Analysis of Relative Ligand Binding Free 
Energies with Loop Closure and Experimental Constraints
Timothy J. Giese, Darrin M. York
J. Chem. Theory Comput. (2021)
DOI: 10.1021/acs.jctc.0c01219

[3] Extension of the Variational Free Energy Profile and Multistate Bennett 
Acceptance Ratio Methods for High-Dimensional Potential of Mean Force Profile 
Analysis
Timothy J. Giese, Şölen Ekesan, Darrin M. York
J. Phys. Chem. A (2021) 125, 4216-4232
DOI: 10.1021/acs.jpca.1c00736

[4] Multireference Generalization of the Weighted Thermodynamic
Perturbation Method
Timothy J. Giese, Jinzhe Zeng, and Darrin M. York
J. Phys. Chem. A (2022) 126, 8519-8533
DOI: 10.1021/acs.jpca.2c06201

