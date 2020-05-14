import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    # Project details
    name="sensei2-shaoheng",
    version="0.0.1",
    description="sample size estimation for cell abundance in single-cell sample",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KChen-lab/sensei",

    # Author details
    author="Shaoheng Liang",
    author_email="sliang3@mdanderson.org",

    # Dependencies
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "scipy>=1.4.0", # Betabinomial is only supported after scipy 1.4.0
        "numpy",
        "matplotlib",
        "pandas",
        "plotly"
    ]
)