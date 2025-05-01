from setuptools import setup, find_packages

setup(
    name="scCNA",  # Name of the package
    version="0.1",  # Package version
    description="A package for detecting and annotating Copy Number Aberrations in single-cell RNA-seq data",
    author="Your Name",  # Replace with your name
    author_email="your.email@example.com",  # Replace with your email
    url="https://github.com/yujiafeng8888/scCNA",  # Replace with your project URL
    packages=find_packages(),  # Automatically find packages in the project
    install_requires=[
        "pandas",
        "numpy",
        "scanpy",
        "seaborn",
        "matplotlib",
        "anndata",
        "typing-extensions",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",  # Minimum Python version required
    include_package_data=True,
    zip_safe=False,
)
