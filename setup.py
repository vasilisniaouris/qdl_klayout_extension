import setuptools
import datetime

from configparser import ConfigParser

# get long description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# read package metadata from setup.cfg
config = ConfigParser()
config.read('setup.cfg')
metadata = dict(config.items('metadata'))
options = dict(config.items('options'))

# make the copyright string by using the current year
_year = str(datetime.datetime.today().year)
copyright_years = _year if _year == metadata['first_released_year'] else f"{metadata['first_released_year']} - {_year}"
_copyright = f"Copyright (C) {copyright_years}, {metadata['author']}"

setuptools.setup(
    name=metadata['name'],
    version=metadata['version'],
    author=metadata['author'],
    author_email=metadata['autor_email'],
    license=metadata['licence'],
    url=metadata['url'],
    copyright=_copyright,
    description=metadata['description'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
    package_dir=dict(config.items('options.package_dir')),
    package_data={key: value.split('\n') for key, value in config.items('options.package_data')},
    include_package_data=bool(options['include_package_data']),
    install_requires=options['install_requires'].split('\n'),
    classifiers=metadata['classifiers'].split('\n'),
    python_requires=metadata['python_requires'],
)
