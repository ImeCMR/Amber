# Conda recipe for AmberTools

This is the `conda` recipe used to build `ambertools` for the Conda Forge ecosystem.

The files included in this directory weere copied over from [AmberTools Feedstock PR#30](https://github.com/conda-forge/ambertools-feedstock/pull/30) on 2020.04.27, corresponding to commit `a51b9ac`. More up-to-date features might be available on that repository in the future.

Limitations at this point:

* This only builds the `serial` version for Linux (`amd64` and `ppc64le`) and MacOS.
* MacOS package does not include GUI tools.
* By default, unit tests are not run for every build. You can change this by editing `conda_build_config.yaml` so `unit_tests` is set to `run` instead of `skip`. Running the full test suite will add all the test files to the final package (~650MB). Unit tests require `csh` to be installed system-wide.

`conda_build_config.yaml` governs the pinned versions for each critical package in the Conda Forge ecosystem. If you want to create 100% compatible CF packages, make sure to update its contents after the first commented block (leave `unit_tests: skip` there!). You can also clear all the pinnings for maximum flexibility, if needed.

# How to build

You might need to edit `meta.yaml` so it points to the latest source release. The `sha256` hash must match as well. Check the comments in that file for more instructions. Then, install `conda build` if needed. It must be in the base `conda` environment!

```
conda activate base
conda install conda-build
```

Then, make sure you have `conda-forge` as part of your channels. Assumming you are running this from the repository root:

```
conda build -c conda-forge recipe_at/
```

# Contact

Check the [AmberTools Feedstock](https://github.com/conda-forge/ambertools-feedstock) for help strictly related to this recipe.