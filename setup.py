from setuptools import setup, find_packages
import pypandoc


try:
    description=pypandoc.convert('README.md', 'rst')
except (IOError, ImportError):
    description=open('README.md').read()

    
# This call to setup() does all the work
setup(
    name="Heavy-Neutrino-Limits",
    version="0.0.2",
    description="A package with a compilation of all publicly available experimental limits on heavy neutrinos.",
    long_description = description,
    long_description_content_type="text/markdown",
    url="https://github.com/mhostert/Heavy-Neutrino-Limits",
    author="Matheus Hostert",
    author_email="mhostert@perimeterinstitute.ca",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    packages=find_packages(),
    include_package_data=True,
    install_requires=["numpy", "matplotlib", "seaborn", "pandas","hepunits"]
)