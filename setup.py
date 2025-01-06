from setuptools import setup, find_packages

setup(
    name='RASP',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'scanpy',
        'squidpy',
        'matplotlib',
        'scipy',
        'scikit-learn',
        'rpy2',
        'mclustpy',
        'igraph',
        'glasbey',
    ],
    author='Ian K Gingerich',
    author_email='ian.gingerich.gr@darmouth.edu',
    description='A computationally efficient method for dimensionality reduction of high-resolution spatial transcriptomics data',
    url='https://github.com/yourusername/my_rasp_package',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
