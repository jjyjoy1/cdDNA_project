1. Introduction

1.1 Overview
Environment Modules allow you to dynamically modify your user environment using modulefiles. Modules let you load and unload different software versions without cluttering your system’s global PATH or managing complex environment variables manually.
EasyBuild is a software build and installation framework that allows you to automate the process of building, installing, and managing (scientific) software on HPC systems or local servers. It automatically generates corresponding modulefiles, ensuring your users can load the installed software easily.
Using both Environment Modules and EasyBuild together provides an efficient, scalable, and reproducible approach for deploying software in HPC or other shared computing environments.

2. Environment Modules

2.1 What Are Environment Modules?
Environment Modules is a tool that lets you modify environment variables (e.g., PATH, LD_LIBRARY_PATH) in a modular way. Modules are typically managed through simple commands: module load <module-name> and module unload <module-name>.

Environment Modules:

Make software setups easier by organizing different versions.
Avoid version conflicts and environment pollution.
Provide a systematic approach to environment management.

2.2 Installing Environment Modules
Note: On many HPC systems, environment modules might already be installed. If your system already has module commands available, you can skip this step.
If you are have sudo user permission, you can set environment modules path on your personal directory

Download the Environment Modules source from the official website:
https://github.com/cea-hpc/modules
Extract and navigate to the source directory:
tar -xzf modules-<version>.tar.gz
cd modules-<version>
Configure the build (prefix where you want Modules to be installed):
./configure --prefix=/opt/modules
Compile and install:
make
make install
Configure shell startup scripts (e.g., .bashrc, .zshrc) to initialize modules. Typically, you add something like:
# Initialize Environment Modules
source /opt/modules/init/bash
Verify installation:
module --version
If you see the version number, your Environment Modules installation was successful.

2.3 Using Environment Modules
List available modules:
module avail
Load a module:
module load <module-name>
Unload a module:
module unload <module-name>
Show currently loaded modules:
module list
Show what a module does (e.g., which environment variables it sets):
module show <module-name>

3. Introduction to EasyBuild

3.1 What Is EasyBuild?
EasyBuild is an open-source framework that simplifies building and installing (scientific) software on HPC or local servers. It automates:

Fetching source code
Handling dependencies
Applying patches
Configuring
Compiling
Installing software
Generating modulefiles
With EasyBuild, you reduce the manual steps required to build software from source. It also helps you maintain a consistent environment across different systems.

###You can also use easybuild to build easybuild environment module, and and load easybuild when you need it. 

3.2 EasyBuild Key Components
EasyBuild Framework – The core scripts that coordinate building and installing software.
EasyBlocks – Python modules that encapsulate the building logic for different software packages.
EasyConfigs – Files that specify how a particular piece of software is to be built. (E.g., versions, compiler toolchains, dependencies.)
3.3 Installing EasyBuild
Prerequisite: Python 3 is recommended for installing EasyBuild.
Install from pip (recommended):
pip install easybuild
Or, if you need a user-based install:

pip install --user easybuild
Alternatively, clone from GitHub:
git clone https://github.com/easybuilders/easybuild-framework.git
cd easybuild-framework
pip install .
Verify installation:
eb --version
You should see the current version of EasyBuild printed.

4. Building & Installing Software with EasyBuild

4.1 Common Workflow
Locate or create an EasyConfig file for the software you want to build (e.g., GROMACS-2021.3-foss-2021a.eb).
You can often find community EasyConfigs in the EasyBuild repository.
Or write a custom EasyConfig if none exists for your software.
Check dependencies. Ensure that:
The compiler toolchain modules (like foss, intel, GCC, etc.) are available.
Any prerequisite modules or library versions are loaded or can be located.
Run EasyBuild to install:
eb GROMACS-2021.3-foss-2021a.eb --installpath=/opt/software --sourcepath=/opt/sources
--installpath defines where EasyBuild places the installed software.
--sourcepath is where EasyBuild looks for the source tarballs (or will store them if it needs to download).
Verify installation & generated module:
Module will typically be generated in (e.g.) /opt/software/modules/all/GROMACS/2021.3-foss-2021a.
Load the module and check that the software runs.

4.2 Using Custom EasyConfig Files
If you have a custom or modified EasyConfig file, place it in a location recognized by EasyBuild. For example:

export EASYBUILD_PREFIX=/opt/easybuild
export EASYBUILD_CONFIGFILES=/path/to/myconfig/easybuild.cfg
export EASYBUILD_BUILDPATH=$EASYBUILD_PREFIX/build
export EASYBUILD_INSTALLPATH=$EASYBUILD_PREFIX/software
export EASYBUILD_SOURCEPATH=$EASYBUILD_PREFIX/sources

eb /path/to/GROMACS-2021.3-foss-2021a.eb


5. Integrating Environment Modules & EasyBuild

One of the strengths of EasyBuild is that it automatically generates modulefiles. By default, EasyBuild uses the Environment Modules naming scheme. This means that after a successful build, you can do:

module use /opt/software/modules/all
module load GROMACS/2021.3-foss-2021a
Your environment now has all the paths set for that specific GROMACS installation.

5.1 Example Directory Structures
A typical directory structure for an EasyBuild-managed environment:

/opt/easybuild
├── sources        # Source tarballs
├── software       # Installed software
└── modules
    └── all        # Modulefiles generated by EasyBuild
Ensure that your module use command points to /opt/easybuild/modules/all so that you can discover the newly built modules:

module use /opt/easybuild/modules/all

6. Best Practices

Centralize your EasyBuild configuration:
Use a system-level config file (often /etc/easybuild.d/ or $HOME/.config/easybuild/) to define the standard paths, compilers, etc.
Keep track of your toolchains:
EasyBuild toolchains (e.g., foss, intel, gompi) define your compiler and MPI stack. Be consistent about which toolchain you use to avoid compatibility issues.
Version control your EasyConfig files:
If you create custom EasyConfig files, store them in a version-controlled repository (e.g., Git) to track changes and share with colleagues.
Avoid environment pollution:
Make sure to module purge or carefully unload modules before building software with EasyBuild, so that you build in a clean environment.
Document:
Each time you build new software, document the build environment, EasyConfig changes, and any custom patches. This ensures reproducibility.
Automate:
Consider using a CI/CD or HPC job scheduling approach to automate large sets of builds with EasyBuild.

7. Example: Full Walk-Through

Below is an example illustrating the entire workflow of setting up a module environment, installing software with EasyBuild, and using the resulting module.

Load your base environment (if Modules is already installed):
# HPC system usually has Modules by default
module load modules
Install EasyBuild if needed:
pip install --user easybuild
Set up EasyBuild environment variables (in ~/.bashrc or session):
export EASYBUILD_PREFIX=$HOME/.local/easybuild
export EASYBUILD_BUILDPATH=$EASYBUILD_PREFIX/build
export EASYBUILD_INSTALLPATH=$EASYBUILD_PREFIX/software
export EASYBUILD_SOURCEPATH=$EASYBUILD_PREFIX/sources
Fetch or create an EasyConfig file:
cd $EASYBUILD_PREFIX
wget https://raw.githubusercontent.com/easybuilders/easybuild-easyconfigs/develop/easybuild/easyconfigs/h/HDF5/HDF5-1.10.7-gompi-2021a.eb
Build & install:
eb HDF5-1.10.7-gompi-2021a.eb
Wait for the compilation and installation to complete.
Use the newly generated module:
module use $EASYBUILD_PREFIX/modules/all
module load HDF5/1.10.7-gompi-2021a
Now your environment is configured to use HDF5 with that specific toolchain.
Verify:
h5cc -showconfig
You should see details matching the built version.

8. Troubleshooting & Tips

Environment Conflicts: If builds fail, check that you are not loading conflicting modules or environment variables from system-level configurations.
Permission Issues: Make sure you have write permissions to your chosen install and modulefile paths.
Check Logs: EasyBuild logs each build in the $EASYBUILD_BUILDPATH directory. Consult the logs if you encounter problems.
Module Naming Issues: Ensure your module purge and module use <path> calls do not conflict with system default modules.

9. References & Further Reading

Environment Modules:
https://modules.readthedocs.io/
EasyBuild:
https://easybuild.io/
https://github.com/easybuilders/easybuild
Conclusion

By combining Environment Modules and EasyBuild, you can maintain a clean, modular, and reproducible software environment on both HPC clusters and local servers. Environment Modules simplifies how end users load and unload different software versions, and EasyBuild streamlines the building process and automatically creates the corresponding modulefiles.

This workflow:

Improves consistency across development, testing, and production environments.
Ensures reproducibility by capturing all build information in EasyConfig files.
Simplifies user experience by letting them just module load <software>.
Feel free to adapt the above structure and examples to suit your organization’s specific paths, HPC scheduler environment, or configuration management tools.


