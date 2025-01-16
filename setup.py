import setuptools

setuptools.setup(
    name='trachoma',
    version='1.0.1dev',
    url='https://www.artrabbit.com/',
    maintainer='ArtRabbit',
    maintainer_email='support@artrabbit.com',
    description='Trachoma transmission model',
    long_description='Set MDA times, initialise infection, and simulate transmission / MDA interventions.',
    packages=setuptools.find_packages(),
    python_requires='>=3.8',
    install_requires=['numpy', 'pandas', 'joblib', 'google-cloud-storage', 'matplotlib', 'openpyxl', 'pytest'],
    include_package_data=True
)
