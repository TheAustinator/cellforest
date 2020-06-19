import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cellforest",
    version="0.0.1",
    author="Austin McKay",
    author_email="austinmckay303@gmail.com",
    description="An interactive single-cell bioinformatics workflow manager",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TheAustinator/cellforest",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)