from setuptools import setup

with open('README.md', 'r') as fp:
    long_description = fp.read()

setup(
    name='ndtest',
    version='0.1',
    description='Multi-dimensional statistical tests with python, including the 2D Kolmogorov–Smirnov test and energy distance statistics.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/syrte/ndtest/',
    keywords=['Fortran', 'Numpy'],
    author='Zhaozhou Li',
    author_email='lizz.astro@gmail.com',
    py_modules=['ndtest'],
    install_requires=['numpy', 'scipy'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
)
