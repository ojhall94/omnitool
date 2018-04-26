from setuptools import setup, find_packages

setup(name='omnitool',
    version='0.3',
    description='A series of classes for photometric and asteroseismic analysis',
    url='https://github.com/ojhall94/omnitool',
    author='Oliver James Hall',
    author_email='ojh251@student.bham.ac.uk',
    license='MIT',
    packages=find_packages(),
    package_data = {'omnitool':['data/*']},
    include_package_data = True,
    install_requires=[
        'numpy',
        'astropy',
        'pandas',
        'dustmaps',
        'pystellibs',
        'synphot',
        'scipy',
        'tqdm'
        ],
    zip_safe=False)
