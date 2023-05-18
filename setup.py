from glob import glob
from os.path import join as pjoin

from setuptools import setup, find_packages

PKG = 'demerit'
DESCRIPTION = "MeerKAT Demerit score tools"

setup(name=PKG,
      description=DESCRIPTION,
      long_description=DESCRIPTION,
      url='https://github.com/ska-sa/katsdpcontim',
      classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "Topic :: Scientific/Engineering :: Astronomy",
      ],
      author='Tom Mauch',
      author_email='tmauch@sarao.ac.za',
      python_requires=">=3.0",
      install_requires=['astropy', 'numpy'],
      scripts=glob(pjoin('scripts', '*.py')),
      packages=find_packages(),
      package_data={PKG: [pjoin('static', '*')]},
)


