import setuptools
import datetime

from configparser import ConfigParser

# get long description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# read package metadata from setup.cfg
_config = ConfigParser()
_config.read('setup.cfg')
_metadata = dict(_config.items('metadata'))
_options = dict(_config.items('options'))

# make the copyright string by using the current year
_year = str(datetime.datetime.today().year)
copyright_years = _year if _year == _metadata['first_released_year'] else f"{_metadata['first_released_year']} - {_year}"
_copyright = f"Copyright (C) {copyright_years}, {_metadata['author']}"

setuptools.setup(
    # The rest are defined in the setup.cfg file!
    copyright=_copyright,
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=setuptools.find_packages(),
)
