import setuptools
#from distutils.core import setup
#long_description = ""
#with open("README.txt", "r") as fh:
 #   long_description = fh.read()

setuptools.setup(
    name='bioformats_io',
    version='1.0',
    description='Bioformats-to-python reader/writer',    
  #  long_description=long_description,
license='GPLv3',
    author='Olivia Mariani',
    author_email='olivia.mariani@idiap.ch',
    packages=['bioformats_io'],
    python_requires='>=3.6.0',
    classifiers=[
        "Programming Language :: Python :: 3",],
    install_requires=["python-dateutil==2.8.0","numpy>=1.14.0", "argparse>=1.4.0", "javabridge>=1.0.15", "python-bioformats>=1.3.2",
                      
]
     )

