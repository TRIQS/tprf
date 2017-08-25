# -*- coding: utf-8 -*-
#
# TRIQS documentation build configuration file

import sys
sys.path.insert(0, "@TRIQS_SPHINXEXT_PATH@/numpydoc")

extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.mathjax',
              'sphinx.ext.intersphinx',
              'sphinx.ext.doctest',
              'sphinx.ext.todo',
              'sphinx.ext.viewcode',
              'sphinx.ext.autosummary',
              'plot_directive',
              'numpydoc']

source_suffix = '.rst'

project = u'TPRF: Two-Particle Response Function tool box for TRIQS'
copyright = u'2017, H. U.R. Strand'
version = '0.0.0'
release = '0.0.0'

math_number_all = True
#mathjax_path = "@TRIQS_MATHJAX_PATH@/MathJax.js?config=default"
templates_path = ['_templates']

#html_theme = 'triqs'
#html_theme_path = ['@TRIQS_THEMES_PATH@']
html_show_sphinx = False
html_context = {'header_title': 'tprf',
                'header_subtitle': 'two-particle response function tools based on the <a class="triqs" style="font-size: 12px" href="http://ipht.cea.fr/triqs">TRIQS</a> library',
                'header_links': [['Install', 'install'],
                                 ['Documentation', 'documentation'],
                                 ['Issues', 'issues'],
                                 ['About TPRF', 'about']]}
html_static_path = ['_static']
html_sidebars = {'index': ['sideb.html', 'searchbox.html']}

htmlhelp_basename = 'TRIQSTPRFdoc'

intersphinx_mapping = {'python': ('http://docs.python.org/2.7', None), 'triqslibs': ('@TRIQS_INVENTORY_DIR@', None)}
