import pathlib
from setuptools import setup
from setuptools import find_packages

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='gpi_core',
    version='2.4.0',
    description="Graphical Programming Interface",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/gpilab/framework",
    author="AbdulRahman Alfayad",
    author_email="alfayad.abdulrahman@mayo.edu",
    license="GNU",
    packages=find_packages(),
    include_package_data=True,
    python_requires='>=3.7',
)