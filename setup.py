from setuptools import find_packages, setup

from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()


def parse_requirements(requirements_path):
    with open(Path(__file__).parent / requirements_path) as f:
        return f.read().splitlines()


requirements = parse_requirements("requirements.txt")
test_requirements = parse_requirements("requirements/requirements-test.txt")
dev_requirements = parse_requirements("requirements/requirements-dev.txt")
print("REQUIREMENTS: ", requirements)

setup(
    name="cellforest",
    version="0.0.1",
    author="Austin McKay",
    author_email="austinmckay303@gmail.com",
    description="An interactive single-cell bioinformatics workflow manager",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TheAustinator/cellforest",
    packages=find_packages(exclude=("tests",)),
    # TODO: this works in "all", but not in `install_reqs` -- fix
    install_requires=requirements,
    extras_require={"all": requirements, "test": test_requirements, "dev": dev_requirements,},
    # data_files=[("cellforest/config/", ["cellforest/config/process_schema.yaml"])],
    package_data={"cellforest": ["config/process_schema.yaml"]},
    include_package_data=True,
    manifest="MANIFEST.in",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: AGPL 3.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
)
