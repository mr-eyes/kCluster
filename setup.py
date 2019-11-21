"""

   ▄█   ▄█▄    ▄████████    ▄███████▄  ▄█  ████████▄     ▄████████    ▄████████
  ███ ▄███▀   ███    ███   ███    ███ ███  ███   ▀███   ███    ███   ███    ███
  ███▐██▀     ███    █▀    ███    ███ ███▌ ███    ███   ███    █▀    ███    ███
 ▄█████▀      ███          ███    ███ ███▌ ███    ███  ▄███▄▄▄      ▄███▄▄▄▄██▀
▀▀█████▄    ▀███████████ ▀█████████▀  ███▌ ███    ███ ▀▀███▀▀▀     ▀▀███▀▀▀▀▀
  ███▐██▄            ███   ███        ███  ███    ███   ███    █▄  ▀███████████
  ███ ▀███▄    ▄█    ███   ███        ███  ███   ▄███   ███    ███   ███    ███
  ███   ▀█▀  ▄████████▀   ▄████▀      █▀   ████████▀    ██████████   ███    ███
  ▀                                                                  ███    ███

"""

from setuptools import setup, find_packages
import sys
from os import path
from src.version import __version__  # pylint: disable=relative-beyond-top-level

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()



setup(
    name='kSpider',
    version=__version__,
    description = "kSpider sequence clustering software.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/mr-eyes/kSpider',
    author='Mohamed Abuelanin (Nile University), Tamer Mansour (UC Davis)',
    author_email='mabuelanin@gmail.com, drtamermansour@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='pairwise-similarity clustering',
    python_requires='>=3.6',
    install_requires=[
        'Click',
        'kProcessor',
    ],
    packages=find_packages(),
    include_package_data=True,
    entry_points='''
        [console_scripts]
        kSpider=src.kSpider:cli
    ''',
    project_urls={
        'Bug Reports': 'https://github.com/mr-eyes/kSpider/issues',
        'Say Thanks!': 'https://saythanks.io/to/mr-eyes',
        'Source': 'https://github.com/mr-eyes/kSpider',
        },
)