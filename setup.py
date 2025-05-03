from setuptools import setup, find_packages

setup(
    name="scCNA",
    version="0.1",
    description="A package for detecting and annotating Copy Number Aberrations in single-cell RNA-seq data",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yujiafeng8888/scCNA",
    packages=find_packages(where="src"),
    package_dir={"": "src"},  
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
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    include_package_data=True,
    zip_safe=False,
)
