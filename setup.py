from setuptools import find_packages, setup

REQUIRED_PACKAGES = [
    "numpy>=1.8.1",
    "dipy>=1.1.1",
    "nipype>=1.4.0",
    "scipy>=1.4.0",
    "nibabel",
    "boto3",
]

setup(
    name="hcp_connectomes",
    packages=find_packages(),
    version="0.1.0",
    description="HCP DWI data preprocessing",
    author="j1c",
    license="Apache License 2.0",
    install_requires=REQUIRED_PACKAGES,
)
