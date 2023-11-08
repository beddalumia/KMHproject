import setuptools

setuptools.setup(
    name='strawberrypy_tmp',
    #authors='Roberta Favata, Nicolas Baù and Antimo Marrazzo',
    version='0.2.0',
    license='MIT License',
    author_email='favata.roberta@gmail.com',
    description='Python package for calculation of topological invariants through single-point formulas and local markers',
    packages=setuptools.find_packages(),
    install_requires=['pythtb==1.8.0',
                      'tbmodels==1.4.3',
                      'opt-einsum==3.3.0',],
)
