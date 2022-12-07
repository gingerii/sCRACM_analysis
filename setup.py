import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sCRACM_analysis",
    version="0.1.1",  # Don't forget to update in sCRACM_analysis/__init__.py
    author="Ian Gingerich",
    author_email="gingerichik@gmail.com",
    description="Package for importing, wrangling, analyzing, and "
                "visualization of brain cell data (optogenetic mapping data). in development",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gingerii/sCRACM_analysis",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
        "Liscense::OSI Approved::MIT Liscense",
    ],
    install_requires=['plotly>= 3.5.0',
                      'nested-lookup>=0.2.18',
                      'pandas>=0.23.4',
                      'numpy>=1.15.3',
                      'scipy>=1.2.0',
                      'matplotlib>=3.1.2',
                      'Pillow>=5.4.1',
                      'seaborn>=0.9.0',
                      'pyqrcode>=1.2.1',
                      'tables>=3.5.2',
                      'xarray>=0.20.2',
                      ]
)

