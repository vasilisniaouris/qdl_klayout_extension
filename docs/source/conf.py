# Check out this example! https://github.com/eeholmes/rtd-github-pages

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import tomli

from urllib.parse import quote


with open("../../pyproject.toml", "rb") as f:
    _config = tomli.load(f)

_project_metadata = _config['project']
_docs_metadata = _config['docs']

# set paths
_basedir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../', _project_metadata['name']))
sys.path.insert(0, _basedir)

# -- Project information -----------------------------------------------------

project = _docs_metadata["fancy-name"]
copyright = _docs_metadata["copyright"]
author = _project_metadata['authors'][0]["name"]

# The full version, including alpha/beta/rc tags
release = _project_metadata['version']

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_search.extension',  # pip install readthedocs-sphinx-search
    # 'sphinx.ext.autodoc',  # sphinx extension for automatic documentation from files. Could not figure it out...
    # 'sphinx.ext.napoleon',  # autodoc for numpy anf google styles
    'autoapi.extension',  # pip install sphinx-autoapi easier than autodoc
    'myst_parser',  # for md support. pip install -U myst-parser
    'sphinx.ext.githubpages',  # to publish on GitHub Pages
    # 'sphinx.ext.linkcode',  # to add source to github page!
    'sphinx.ext.viewcode',  # to see source code on the API reference
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ['.rst', '.md']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# -- Options for automated documentation -------------------------------------

# autoapi extension settings
autoapi_type = 'python'
autoapi_dirs = [f"../../{_project_metadata['name']}"]
autoapi_add_toctree_entry = True
autoapi_python_class_content = 'both'  # both class doc and __init__ doc
autoapi_member_order = 'groupwise'

# napoleon extension settings
# napoleon_google_docstring = True
# napoleon_numpy_docstring = True
# napoleon_include_init_with_doc = False
# napoleon_include_private_with_doc = False
# napoleon_include_special_with_doc = True
# napoleon_use_admonition_for_examples = False
# napoleon_use_admonition_for_notes = False
# napoleon_use_admonition_for_references = False
# napoleon_use_ivar = False
# napoleon_use_param = True
# napoleon_use_rtype = True
# napoleon_preprocess_types = False
# napoleon_type_aliases = None
# napoleon_attr_annotations = True

# -- Options for MyST parser -------------------------------------------------

# myst_heading_anchors = 3

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

# default theme
# html_theme = 'alabaster'

# Basic Read the Docs theme. Install with pip install sphinx-rtd-theme
# Currently search is broken: e.g. https://stackoverflow.com/questions/52474177/read-the-docs-search-broken
# html_theme = 'sphinx_rtd_theme'

# Cute theme, not great for API. Install with pip install git+https://github.com/bashtage/sphinx-material.git
#  for more info: https://bashtage.github.io/sphinx-material/index.html
# html_theme = 'sphinx_material'

# Another cute theme, not greate for API. Install with pip install sphinx-book-theme
# html_theme = 'sphinx_book_theme'

# Good alternative for read the docs, nice for API. Install with pip install furo
html_theme = 'furo'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_show_sourcelink = False  # the documentation source files, not the code itself.

# for rtd theme
# html_theme_options = {
#     'style_external_links': True,
#     'style_nav_header_background': '#4B0082',
#     # Toc options
#     'sticky_navigation': True,
#     'navigation_depth': -1}

# for material theme
# html_theme_options = {
#     'nav_title': project,
#     # 'theme_color': '4D1A7F',
#     'color_primary': 'indigo',
#     'color_accent': 'red',
#     'globaltoc_depth': 2,
# }

# for books theme
# html_theme_options = {
#     "repository_url": metadata['url'],
#     "use_repository_button": True,
# }

# for furo  # https://pradyunsg.me/furo/customisation/
html_theme_options = {
    "sidebar_hide_name": False,
    "source_repository": _project_metadata['urls']["repository"],
    "source_branch": "main",
    "source_directory": "docs/",
 }

# for sphinx specifically
html_theme_options.update({
    'navigation_with_keys': True,
})

# html_sidebars = {
#     "**": ["logo-text.html", "globaltoc.html", "localtoc.html", "searchbox.html"]
# }


def linkcode_resolve(domain, info):
    # https://stackoverflow.com/questions/48298560/how-to-add-link-to-source-code-in-sphinx
    if domain != 'py' or not info['module']:
        return None
    filename = quote(info['module'].replace('.', '/'))
    anchor = "#:~:text=" + quote(info["fullname"].split(".")[-1]) if "fullname" in info else ""

    # github
    result = f"{_project_metadata['urls']['repository']}/blob/master/%s.py%s" % (filename, anchor)
    return result
