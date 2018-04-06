from setuptools import setup

setup(name='omnitool',
    version='0.1',
    description='A series of classes for photometric and asteroseismic analysis',
    url='https://github.com/ojhall94/omnitool',
    author='Oliver James Hall',
    author_email='ojh251@student.bham.ac.uk',
    license='MIT',
    packages=['omnitool']
    install_requires=[
        'numpy',
        'astropy',
        'pandas',
        'dustmaps',
        'pystellibs',
        'synphot',
        'scipy'
        ]
    zip_safe=False)
