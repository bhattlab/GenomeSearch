from setuptools import setup, find_packages

setup(
    name="genomesearch",
    version='0.0.8_dev',
    description='A command line tool to quickly search for closely related microbial genomes using a marker-gene based approach.',
    url='https://github.com/bhattlab/GenomeSearch',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'wget==3.2',
        'click==7.0',
        'biopython==1.76'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'genomesearch = genomesearch.main:cli',
        ],
}
)
