from distutils.core import setup

setup(
    name='SSOPipeline',
    author='Max Mahlke',
    author_email='max.mahlke@cab.inta-csic.es',
    description='SSO Recovery Pipeline',
    url='https://github.com/maxmahlke/SSO_Pipeline.git',
    version='0.9dev',
    packages=['ssopipeline',],
    license='GNU Affero General Public License',
    long_description=open('README.md').read(),
)