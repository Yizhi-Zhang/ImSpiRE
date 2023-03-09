from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')


setup(
    name='imspire',  # Required

    version='1.1.2',  # Required

    description='ImSpiRE is a python script (python 3.8+) for spatial resolution enhancenment by solving the entropic regularized fused Gromov-Wasserstein transport (FGW) problem for in situ capturing (ISC) datasets.',  # Optional

    long_description=long_description,  # Optional

    long_description_content_type='text/markdown',  # Optional (see note above)

    url='https://github.com/TongjiZhanglab/ImSpiRE',  # Optional

    author='Yuwei Hua',  # Optional

    author_email='ywhua@tongji.edu.com',  # Optional

    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',

        # Pick your license as you wish
        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 3.8',
    ],

    package_dir={'': 'src'},
    packages=[''],
    scripts=['bin/ImSpiRE'],

    python_requires='>=3.8',

    install_requires=[
        'pandas==1.4.1',
        'numpy==1.23.4',
        'matplotlib==3.1.3',
        'scanpy==1.9.0',
        'scipy==1.9.3',
        'scikit-learn==1.1.3',
        'numba==0.53.0',
        'backgroundremover==0.1.9',
        # 'cellprofiler==4.2.1',
        'scikit-image==0.19.3',
        'anndata==0.8.0',
        'tifffile',
        'POT',
        'opencv-python'
    ],

)
