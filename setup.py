import setuptools

setuptools.setup(
    name='trachoma',
    version='0.0.2',
    url='https://www.artrabbit.com/',
    maintainer='ArtRabbit',
    maintainer_email='support@artrabbit.com',
    description='Trachoma transmission model',
    long_description='Set MDA times, initialise infection, and simulate transmission / MDA interventions.',
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=['numpy', 'pandas', 'joblib', 'google-cloud-storage', 'matplotlib', 'openpyxl'],
    include_package_data=True
)
