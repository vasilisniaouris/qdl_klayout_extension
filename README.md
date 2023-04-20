# Description 
The quantum defect laboratory (QDL) KLayout Extension provides additional functionalities complementary to the KLayout Python API.

# Installation
To install the package you can either
1.  directly install it from GitHub you can use the following command:
    ~~~shell
    pip install git+https://github.com/vasilisniaouris/qdl_klayout_extension.git
    ~~~
    This will install the latest version of the package from the master branch. 

    If you want to install a specific release, you can use the tag name instead of master, like this:
    ~~~shell
    pip install git+https://github.com/vasilisniaouris/qdl_klayout_extension.git@v0.1.0
    ~~~
    Replace v0.1.0 with the tag name of the release you want to install.

2. Or, Alternatively, you can clone the repository:

    ~~~shell
    git clone https://github.com/vasilisniaouris/qdl_klayout_extension.git
    ~~~
    
    Navigate to the cloned directory:
    ~~~shell
    cd qdl_klayout_extension
    ~~~
    
    And, finally, use pip to install the package:
    ~~~shell
    pip install .
    ~~~

And you are all done. The qdl_klayout_extension should be available to you as a python module.

# Dependencies
The QDL KLayout Extension requires the following dependencies:

~~~
klayout>=0.28.6
matplotlib>=3.7.1
multipledispatch>=0.6.0
numpy>=1.24.2
pint>=0.20.1
~~~

# License
The QDL KLayout Extension is released under the GNU GPL v3 license. See LICENSE for more information.
Find a copy of the GNU General Public License [here](https://www.gnu.org/licenses/gpl-3.0.html).

# Copyright
Copyright (C) 2023, Vasilis Niaouris

# Version Control Log
[Found here](./version_control_log.md).

# Credits
The QDL KLayout Extension was created by Vasilis Niaouris under the Quantum Defect Laboratory at the University of Washington.  

This module was inspired and supported by code written by Alan Logan and Nicholas Yama.

# Contact
If you have any questions or comments, please feel free to contact me at vasilisniaouris*gmail.com.
