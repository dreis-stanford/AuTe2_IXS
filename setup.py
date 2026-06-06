from setuptools import setup, find_packages

setup(
    name="aute2-ixs",
    version="0.1.0",
    description="Phonon dispersion and IXS cross-section calculator",
    author="Your Name",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "scipy>=1.7.0",
        "matplotlib>=3.3.0",
    ],
    python_requires=">=3.8",
    entry_points={
        'console_scripts': [
            'aute2-ixs=analyze_q:main',
            'si-ixs=analyze_si:main',
        ],
    },
)
